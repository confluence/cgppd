#include "Simulation.h"

Simulation::Simulation() : waitingThreads(0), exchanges(0), tests(0),  totalExchanges(0), totalTests(0), exchangeFraction(0), accumulativeExchangeFraction(0), offset(0), steps(0), thread_created(false)
{
    REMCRng = gsl_rng_alloc(gsl_rng_mt19937);
    _300kReplica = &replica[0];

    // thread is created in init, because parameters.threads is needed
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_cond_init (&waitingThreadCond, NULL);
    pthread_cond_init (&waitingReplicaExchangeCond, NULL);
    pthread_mutex_init(&waitingCounterMutex, NULL);
    pthread_mutex_init(&writeFileMutex, NULL);
#ifdef VERBOSE_THREAD_LOGGING
    pthread_mutex_init(&logMutex, NULL);
#endif

    // TODO: move these to an init method on aminoAcidData
    LOG(ALWAYS, "Loading amino acid data: %s\n", AMINOACIDDATASOURCE);
    aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);

    LOG(ALWAYS, "Loading pair lookup table: %s\n", LJPDSOURCE);
    aminoAcidData.loadLJPotentialData(LJPDSOURCE);

#ifdef displayLJpotentials
    printPotentials(&aminoAcidData);
#endif

#if USING_CUDA
    cout << "Initialising CUDA. "<< endl ;
    cuInit(0);
    cudaInfo();
#endif // USING_CUDA

#if INCLUDE_TIMERS
    cutCreateTimer(&RELoopTimer);
    cutCreateTimer(&MCLoopTimer);
#endif
}

Simulation::~Simulation()
{
    gsl_rng_free(REMCRng);
    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&waitingCounterMutex);
    pthread_mutex_destroy(&writeFileMutex);
#ifdef VERBOSE_THREAD_LOGGING
    pthread_mutex_destroy(&logMutex);
#endif
    pthread_cond_destroy(&waitingThreadCond);
    pthread_cond_destroy(&waitingReplicaExchangeCond);
    if (thread_created)
    {
        delete [] thread;
        delete [] data;
    }

#if INCLUDE_TIMERS
    cutDeleteTimer(RELoopTimer);
    cutDeleteTimer(MCLoopTimer);
#endif
}

void Simulation::init(int argc, char **argv, int pid)
{
    // Get the parameters

    parameters.pid = pid;

    getArgs(argc, argv);
    loadArgsFromFile();
    getArgs(argc, argv); // second pass to override any variables if doing performance tests

    // sanity check and calculation of some secondary parameters
    check_and_modify_parameters();

    // Create the initial replica, and compare CPU and GPU potential

    // TODO: remove magic number; make initial array size a constant
    initialReplica.init_first_replica(parameters.mdata, aminoAcidData, parameters.bound, 30);
    cout << "Loaded: " << initialReplica.residueCount << " residues in " << initialReplica.moleculeCount << " molecules:\n";

#if USING_CUDA
        // set box size
    if (parameters.auto_blockdim)
    {
        (initialReplica.residueCount < 1024) ? parameters.cuda_blockSize = 32 : parameters.cuda_blockSize = 64;
    }

    initialReplica.setBlockSize(parameters.cuda_blockSize);
#endif

    for (int i=0; i < initialReplica.moleculeCount; i++)
    {
        printf("%2d: %3d %12.3f A^3 %s\n",i , initialReplica.molecules[i].residueCount, initialReplica.molecules[i].volume,initialReplica.molecules[i].filename);
    }
    initialReplica.countNonCrowdingResidues();
    printf("counted : %3d complex residues\n",initialReplica.nonCrowderResidues);

#if INCLUDE_TIMERS
    initialReplica.initTimers();
#endif

    float p = initialReplica.E();
    float pnc = initialReplica.E(&initialReplica.molecules[0],&initialReplica.molecules[1]);
    // TODO: important: for reasons which are temporarily unclear to me, we need to set this initial potential and copy it to the child replicas, otherwise the GL display doesn't move. o_O
    initialReplica.potential = p;

    cout << initialReplica.countpairs() << " residue interaction pairs." << endl;
    printf("CPU initial energy value: %12.8f (%12.8f) kcal/mol\n", p, pnc);

#if USING_CUDA
    CUDA_setBoxDimension(parameters.bound);

#if CUDA_STREAMS
    cudaStreamCreate(&initialReplica.cudaStream);
#endif

    initialReplica.ReplicaDataToDevice();
    cudaMalloc((void**)&initialReplica.device_LJPotentials, LJArraySize);
    copyLJPotentialDataToDevice(initialReplica.device_LJPotentials, &aminoAcidData);

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    bindLJTexture(initialReplica.device_LJPotentials);
#endif

    double GPUE = initialReplica.EonDevice();
    double GPUE_NC = initialReplica.EonDeviceNC();
    initialReplica.potential = GPUE;

    printf("GPU initial energy value: %12.8f (%12.8f) kcal/mol\n", GPUE, GPUE_NC);
    printf("Absolute error:           %12.8f (%12.8f) kcal/mol\n", abs(GPUE-p), abs(GPUE_NC-pnc));
    printf("Relative error:           %12.8f (%12.8f) kcal/mol\n", abs(p-GPUE)/abs(p), abs(pnc-GPUE_NC)/abs(pnc));

#if CUDA_STREAMS
    cudaStreamSynchronize(initialReplica.cudaStream);
    cudaStreamDestroy(initialReplica.cudaStream);
#endif

    if (abs(p-GPUE)/abs(p)>0.01)
    {
        cout << "!\n! WARNING: Inconsistency between GPU and CPU result, check simulation \n!" << endl;
    }
    initialReplica.FreeDevice();
    cutilCheckMsg("Error freeing device");
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    unbindLJTexture();
    cutilCheckMsg("Error freeing texture");
#endif

    cudaFree(initialReplica.device_LJPotentials);
#endif

    if (initialReplica.nonCrowderCount < initialReplica.moleculeCount)
    {
#if REPULSIVE_CROWDING
        cout << "-!- Crowding modelled using: u(r) = (6/r)^(12)." << endl;
#else
        cout << "-!- Crowding is modelled using the full potential calculations." << endl;
#endif
    }

    // now set up all the replicas

    geometricTemperature = pow(double(parameters.temperatureMax/parameters.temperatureMin),double(1.0/double(parameters.replicas-1)));
    geometricTranslate = pow(double(MAX_TRANSLATION/MIN_TRANSLATION),double(1.0/double(parameters.replicas-1)));
    geometricRotation = pow(double(MAX_ROTATION/MIN_ROTATION),double(1.0/double(parameters.replicas-1)));

    for (size_t i = 0; i < parameters.replicas; i++)
    {
        replica[i].init_child_replica(initialReplica, int(i), geometricTemperature, geometricRotation, geometricTranslate, parameters);
        printf ("Replica %d %.3f %.3f %.3f\n", int(i), replica[i].temperature, replica[i].translateStep, replica[i].rotateStep);
    }

    // Temperature map stuff

    // create a lookup map for replica temperatures to allow for sampling temperatures and in place exchanges
    // samples are by temperature not replica index, hence temperatures get shuffled about in replica[],
    // so use replicaTemperaturesmap[] to track where a specific temperature is

    // build a map of lookups for temperature lookups while sampling
    for (size_t i = 0; i < parameters.replicas; i++)
    {
        // lookup of where replica with temperature x is
        position_of_temperature[replica[i].temperature] = i;
        // lookup of which replica has temperature x
        temperature[i] = replica[i].temperature;
    }

    // note which replica is the 300K replica
    // TODO: verify that there really is always a 300k replica
    _300kReplica = &replica[position_of_temperature[300.0f]];

    // Thread stuff

    thread = new pthread_t[parameters.threads];
    data = new SimulationData[parameters.threads];
    thread_created = true;

    // Set up parameters to pass to threads
    for (int i = 0; i < parameters.threads; i++)
    {
        data[i].replica = replica;
        data[i].replicaCount = parameters.replicas;
        data[i].index = i;
        data[i].threads = parameters.threads;
        data[i].streams = parameters.streams;
        data[i].MCsteps = parameters.MCsteps;
        data[i].REsteps = parameters.REsteps;
        data[i].MC_steps_per_RE = parameters.MCsteps/parameters.REsteps;
        data[i].sampleFrequency = parameters.sampleFrequency;
        data[i].sampleStartsAfter = parameters.sampleStartsAfter;
        data[i].bound = parameters.bound;
        data[i].streams = parameters.streams;
        data[i].waitingThreadCount = &waitingThreads;
        data[i].fractionBound = fractionBoundFile;
        data[i].boundConformations = boundConformationsFile;

        data[i].waitingCounterMutex = &waitingCounterMutex;
        data[i].writeFileMutex = &writeFileMutex;
#ifdef VERBOSE_THREAD_LOGGING
        data[i].logMutex = &logMutex;
#endif
        data[i].waitingThreadCond = &waitingThreadCond;
        data[i].waitingReplicaExchangeCond = &waitingReplicaExchangeCond;

        // assign gpus in rotation per thread, t0 = gpu0, t1 = gpu1 etc
        // % #gpus so they share if threads > gpus
        // will perform best if threads:gpus = 1:1
        data[i].GPUID = i % parameters.gpus;
        cout << "Assign thread " << i << " to GPU " << data[i].GPUID << endl;

        data[i].max_replicas_per_thread = parameters.max_replicas_per_thread;

        // we have already checked that the replicas can be distributed like this
        data[i].replicas_in_this_thread = parameters.max_replicas_per_thread;
        if (i == parameters.threads - 2) // last thread
        {
            data[i].replicas_in_this_thread -= parameters.replicas - parameters.max_replicas_per_thread * parameters.threads;
        }

    }

    // File stuff

    int make_dirs = system("mkdir -p output output_pdb checkpoints");

    cout << "Last CUDA error: " << cudaGetErrorString(cudaGetLastError()) << endl;
    cout.flush();
}

void Simulation::run()
{
    // we can't use pthreads and CUDA at the moment, but we are going to use streams

    LOG(ALWAYS, "Beginning simulation\n");

    // TODO: add stuff for resuming here

    if (initialReplica.moleculeCount == 0)  // make sure something is loaded
    {
        cout << "-!- No molecules loaded. -!-" << endl;
        cout << "-!- Aborting run. -!- " << endl;
        return;
    }

    cout << "Output files will be prefixed by " << parameters.prefix << "_" << parameters.pid << endl;
    writeFileIndex();
    initSamplingFiles();

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(RELoopTimer) );
#endif

    printf ("--- Starting simulation ---.\n");

    pthread_mutex_lock(&waitingCounterMutex);

    //TODO: assign streams as a function of the number of GPUs and threads

    //spawn N threads to do the Monte-Carlo mutations
    for (int i = 0; i < parameters.threads; i++)
    {
        //start the MC loops
        pthread_create(&thread[i], &attr, MCthreadableFunction, (void*)&data[i]);
    }

    while (++steps < parameters.REsteps)  // until enough steps taken
    {
#if INCLUDE_TIMERS
        CUT_SAFE_CALL( cutStartTimer(MCLoopTimer) );
#endif

        pthread_cond_wait(&waitingReplicaExchangeCond, &waitingCounterMutex);

#if INCLUDE_TIMERS
        CUT_SAFE_CALL( cutStopTimer(MCLoopTimer) );
#endif

        //if the sampling has started
        if (parameters.MCsteps/parameters.REsteps*(steps-1) >= parameters.sampleStartsAfter)
        {
            // do fraction bound and acceptance ratio

            fprintf(fractionBoundFile, "%9d: ", parameters.MCsteps / parameters.REsteps * steps);
            fprintf(acceptanceRatioFile, "%9d: ", parameters.MCsteps / parameters.REsteps * steps);

            for(map<int, float>::const_iterator iterator = temperature.begin(); iterator != temperature.end(); ++iterator)
            {
                int index = position_of_temperature[iterator->second];
                replica[index].fraction_bound(fractionBoundFile);
                replica[index].acceptance_ratio(acceptanceRatioFile);
            }

            fprintf(fractionBoundFile,"\n");
            fflush(fractionBoundFile);
            fprintf(acceptanceRatioFile,"\n");
            fflush(acceptanceRatioFile);
        }

        for (int ti = offset; ti < parameters.replicas - 1; ti += 2)
        {
            int tj = ti + 1;

            // ti and tj are labels of replicas which correspond to particular temperatures; ti and tj are actual indices within the replica array
            // when we exchange replicas we swap temperatures and labels, keeping almost everything else in place

            int i = position_of_temperature[temperature[ti]];
            int j = position_of_temperature[temperature[tj]];

            double delta = (1.0/replica[i].temperature - 1.0/replica[j].temperature) * (replica[i].potential - replica[j].potential) * (4184.0f/Rgas);

            if (gsl_rng_uniform(REMCRng) < min(1.0, exp(delta)))
            {
                replica[i].exchangeReplicas(replica[j]);

                position_of_temperature[replica[i].temperature] = i;
                position_of_temperature[replica[j].temperature] = j;

                exchanges++; // sampling
            }

            tests++;// tests performed samplng
        }

        offset = 1 - offset; // switch the swap neighbours in the replica exchange

        waitingThreads = 0;
        pthread_cond_broadcast(&waitingThreadCond);


        // update pointer to 300k replica
        _300kReplica = &replica[position_of_temperature[300.0f]];

#if GLVIS
        GLreplica = _300kReplica;
        GlutDisplay();
#endif

        LOG(ALWAYS, "Replica Exchange step %d of %d complete (Fraction bound @ 300K: %f)\n", steps, parameters.REsteps, _300kReplica->accumulativeFractionBound);

        exchange_frequency();
    }

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(MCLoopTimer) );
#endif

    pthread_cond_broadcast(&waitingThreadCond);
    pthread_mutex_unlock(&waitingCounterMutex);  // release the mutex so MC threads can continue.

    printf ("--- Replica Exchanges Complete.---\n");
    printf ("--- Waiting for threads to exit. ---\n");

    // join the threads that have finished
    for (int i=0; i<parameters.threads; i++)
        pthread_join(thread[i],NULL);

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(MCLoopTimer) );
#endif

    printf ("--- All threads complete.---\n");
    fprintf(fractionBoundFile,"%9d: ", parameters.MCsteps);
    fprintf(acceptanceRatioFile,"%9d: ", parameters.MCsteps);

    // final fraction bound and acceptance ratio
    for(map<int, float>::const_iterator iterator = temperature.begin(); iterator != temperature.end(); ++iterator)
    {
        int index = position_of_temperature[iterator->second];
        replica[index].fraction_bound(fractionBoundFile);
        replica[index].acceptance_ratio(acceptanceRatioFile);
    }

    fprintf(fractionBoundFile,"\n");
    fflush(fractionBoundFile);
    fprintf(acceptanceRatioFile,"\n");
    fflush(acceptanceRatioFile);

    printf ("--- Simulation finished.---.\n\n");


#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(RELoopTimer) );
    printf("Simulation Timers\n");
    printf("MC Loop:   Tot  %10.5f ms  Ave %10.5fms (%d steps, %d replicas, %d threads, %d streams)\n"      ,cutGetTimerValue(MCLoopTimer),cutGetTimerValue(MCLoopTimer)/float(parameters.REsteps),parameters.MCsteps,parameters.replicas,parameters.threads,parameters.streams);
    printf("Simulation:     %10.5f ms  (%d exchanges)\n"    ,cutGetTimerValue(RELoopTimer),parameters.REsteps);

#endif

    pthread_mutex_lock(&writeFileMutex);
    closeSamplingFiles();
    pthread_mutex_unlock(&writeFileMutex);

    LOG(ALWAYS, "Simulation done\n");

#if INCLUDE_TIMERS
    for (size_t i=0; i<parameters.replicas; i++)
    {
        cout << "Replica " << i << " timers" << endl;
        replica[i].printTimers();
    }
#endif

    return;
}
// END OF run

void Simulation::exchange_frequency()
{
    totalExchanges += exchanges;
    totalTests += tests;

    exchangeFraction = float(exchanges)/float(tests);
    accumulativeExchangeFraction = float(totalExchanges)/float(totalTests);

    fprintf(exchangeFrequencyFile, "%10d %6.4f %6.4f\n", steps, accumulativeExchangeFraction, exchangeFraction);

    tests = 0;
    exchanges = 0;
}

void setup_CUDA(int device_id, float box_dimension, float * device_LJ_potentials, AminoAcids * amino_acid_data)
{
    // initialise cuda for use in this thread
    cuInit(0);
    cutilCheckMsg("Failed to initialise CUDA runtime.");

    cudaSetDevice(device_id);
    cutilCheckMsg("Failed to pick device for the CUDA runtime.");

    CUDA_setBoxDimension(box_dimension);
    cutilCheckMsg("Failed to copy box dimensions to GPU.");

    // copy the LJpotentials to gpu memory in this thread context so it can access it
    cudaMalloc((void**)&device_LJ_potentials, LJArraySize);
    cutilCheckMsg("Failed to allocate contact potential memory on the GPU.");

    copyLJPotentialDataToDevice(device_LJ_potentials, amino_acid_data);
    cutilCheckMsg("Failed to copy contact potentials to device.");

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    bindLJTexture(device_LJ_potentials);
#endif
}

void teardown_CUDA(float * device_LJ_potentials)
{
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    unbindLJTexture();
#endif

    cudaFree(device_LJ_potentials);
}

#if CUDA_STREAMS
void setup_CUDA_streams(cudaStream_t * streams, int * streams_per_thread, int num_streams, int num_threads, int num_replicas_in_thread)
{
    // the stream/replica ration must be a whole number otherwise there will be lots of waste, ie dormant streams etc
    int replicas_per_stream = int(ceil(float(num_replicas_in_thread)/float(num_streams / num_threads)));
    *streams_per_thread  = replicasInThisThread/replicas_per_stream;

    for (int i = 0; i < *streams_per_thread; i++)
    {
        cudaStreamCreate(&streams[i]); // TODO: is this right?
    }
}

void teardown_CUDA_streams(cudaStream_t * streams, int streams_per_thread)
{
    for (int i = 0; i < streams_per_thread; i++)
    {
        cudaStreamDestroy(streams[i]);
    }
}
#endif

void *MCthreadableFunction(void *arg)
{

    SimulationData *data = (SimulationData *) arg;
    long threadIndex = data->index;
    Replica *replica = data->replica;
    // TODO: this is a temporary patch; we need to clean up the indices
    int tReplicas = data->max_replicas_per_thread;
    int replicasInThisThread = data->replicas_in_this_thread;

    int replica_offset = threadIndex*tReplicas;

    printf ("--- Monte-Carlo thread %d running. ---.\n", int(threadIndex + 1));

// TODO: why do we do this?
#if GLVIS
    GLreplica = &replica[0];
#endif

#ifdef VERBOSE_THREAD_LOGGING
    pthread_mutex_lock( data->logMutex );
    for (int tx = 0; tx < replicasInThisThread; tx++)
    {
        cout << "+ Thread " << threadIndex << " running replica " << tx + threadIndex * tReplicas << endl;
    }
    pthread_mutex_unlock( data->logMutex );
#endif

#if USING_CUDA
    float * deviceLJpTmp;
    setup_CUDA(data->GPUID, replica->boundingValue, deviceLJpTmp, &replica[threadIndex * tReplicas].aminoAcids);

    // create streams for each subsequent MCSearch
    // ensure the replica can find the lookup table
    // initialise data on the device and copy the initial batch

#if CUDA_STREAMS
    cudaStream_t streams[16];   // reserve 16 stream slots but create only as many as needed
    int streamsPerThread(0);
    setup_CUDA_streams(streams, &streamsPerThread);
#endif

    for (int tx = 0; tx < replicasInThisThread; tx++)
    {
        replica[tx + replica_offset].setup_CUDA(deviceLJpTmp);
#if CUDA_STREAMS
        replica[tx + replica_offset].setup_CUDA_streams(streams, tx % streamsPerThread);
#endif
    }
#endif // USING_CUDA

    // run the MC loop
#if !CUDA_MC && !CUDA_STREAMS   // if no MC on device
    int mcstep = 0;
    while (mcstep < data->MCsteps)
    {
        //cout << "Starting mc loop at " << mcstep << " steps" << endl;
        // do all MC steps before RE
        for (int tx = 0; tx < replicasInThisThread; tx++)
        {
            // to sample at the correct rate split this into another set of loops
            for (int s = 0; s < data->MC_steps_per_RE / data->sampleFrequency; s++)
            {
                replica[tx + replica_offset].MCSearch(data->sampleFrequency);

                if (mcstep + s * data->sampleFrequency >= data->sampleStartsAfter)
                {
                    // if E < -1.1844 kcal/mol then its bound
                    replica[tx + replica_offset].sample(data->boundConformations, mcstep + s*data->sampleFrequency, BOUND_ENERGY_VALUE, data->writeFileMutex);
                }
            }
        }
        mcstep += data->MC_steps_per_RE;
        //cout << "thread " << threadIndex << " waits for mutex " << endl;

        // do replica exchange
        pthread_mutex_lock(data->waitingCounterMutex);                 // lock the counter
        //cout << "thread " << threadIndex << " gets mutex " << endl;

        data->waitingThreadCount[0] = data->waitingThreadCount[0] + 1;
        if (data->waitingThreadCount[0] == data->threads)                 // if all threads in waiting state
        {
            pthread_cond_signal(data->waitingReplicaExchangeCond);
            //cout << "thread " << threadIndex << " signals parent" << endl;
        }
        if (mcstep < data->MCsteps)                       // wait if another MC loop must be done
        {
            //cout << "thread " << threadIndex << " waits after " << mcstep << endl;

            pthread_cond_wait(data->waitingThreadCond, data->waitingCounterMutex); // wait for replica exchange to finish.
            // NB! unlock the mutex locked upon the condition above
            // being met because this will unblock all threads such
            // that they will continue concurrently, if not unlocked
            // other threads will run sequentially
            //cout << "thread " << threadIndex << " releases mutex after wait" << endl;

            //  pthread_mutex_unlock(data->waitingCounterMutex);
        }
        pthread_mutex_unlock(data->waitingCounterMutex);


    } // continue MC
#endif


#if CUDA_STREAMS    // use streams
    // TODO: move this to Replica?
    int mcstep = 0;
    //int sampleIn = data->sampleStartsAfter+1;
    while (mcstep < data->MCsteps)
    {
        // to sample at the correct rate split this into another set of loops
        for (int mcx = 0; mcx < data->MC_steps_per_RE; mcx++) // at each mc step
        {
            for (int index = 0; index < replicasInThisThread; index += replicasPerStream)
            {
                for (int rps = 0; rps < replicasPerStream; rps++)
                {
                    // batch replicas such that no stream is shared per batch
                    replica[replica_offset + index + rps].MCSearchMutate();
                    replica[replica_offset + index + rps].MCSearchEvaluate();

                }
                for (int rps = 0; rps < replicasPerStream; rps++)
                {
                    replica[replica_offset + index + rps].MCSearchAcceptReject();
                }

                for (int rps = 0; rps < replicasPerStream; rps++)
                {
                    //sampleAsync
                    if (mcstep % data->sampleFrequency == 0 && mcstep >= data->sampleStartsAfter) // when enough steps are taken && sampleFrequency steps have passed
                    {
                        replica[replica_offset + index + rps].sample(data->boundConformations, mcstep + mcx * data->sampleFrequency, BOUND_ENERGY_VALUE, data->writeFileMutex);

                        //if (abs(replica[replica_offset+index+rps].temperature-300.0f)<1.0f)
                        //  cout << "Sampled: t=" << replica[replica_offset+index+rps].temperature << " mcstep=" << mcstep << endl;
                    }

                }
            }
            mcstep++;
        }

        // do replica exchange
        pthread_mutex_lock(data->waitingCounterMutex);                 // lock the counter
        data->waitingThreadCount[0] = data->waitingThreadCount[0] + 1;
        if (data->waitingThreadCount[0] == data->threads)
        {
            pthread_cond_signal(data->waitingReplicaExchangeCond);
        }
        if (mcstep < data->MCsteps)                                           // wait if another MC loop must be done
        {
            pthread_cond_wait(data->waitingThreadCond, data->waitingCounterMutex);     // wait for replica exchange to finish.
        }
        pthread_mutex_unlock(data->waitingCounterMutex);

    } // continue MC
#endif



#if USING_CUDA
    // sync all streams and free gpu memory
    for (int tx = 0; tx < replicasInThisThread; tx++)
    {
#if CUDA_STREAMS
        replica[tx + replica_offset].teardown_CUDA_streams();
#endif
        replica[tx + replica_offset].teardown_CUDA();
    }

#if CUDA_STREAMS
    teardown_CUDA_streams(streams, streamsPerThread);
#endif

    teardown_CUDA(deviceLJpTmp);
#endif
    printf (" >>> Monte-Carlo thread %d exited.\n",int(threadIndex));
    return 0;
}
// END OF *MCthreadableFunction

void Simulation::initSamplingFile(const char * name, FILE ** file_addr)
{
    char filename[256];
    FILE * file;

    sprintf(filename, "output/%s_%d_%s", parameters.prefix, parameters.pid, name);
    file = fopen (filename,"a+"); // attempt append a file of the same name
    if (!file)
    {
        file = fopen(filename, "wt"); // create if that failed, and open and append failed (eg: file permissions r--)
    }
    if (!file)
    {
        printf("Cannot open/create file: %s.\n", filename);
        return;
    }

    *file_addr = file;
}

void Simulation::closeSamplingFile(const char * name, FILE ** file_addr)
{
    fclose(*file_addr);
    LOG(ALWAYS, ">>> Closing sampling file: output/%s_%d_%s\n", parameters.prefix, parameters.pid, name);
}

void Simulation::initSamplingFiles()
{
    initSamplingFile("fractionBound", &fractionBoundFile);
    initSamplingFile("boundconformations", &boundConformationsFile);
    initSamplingFile("acceptance_ratios", &acceptanceRatioFile);
    initSamplingFile("exchange_freq", &exchangeFrequencyFile);

    fprintf(fractionBoundFile, "# Fraction Bound\n#Iteration InstantaneousAve CumulativeAve\n");
    // TODO: see what we're actually printing now
    fprintf(boundConformationsFile, "# iteration; complex free energy (replica free energy); temperature;\n# molecule: pdb file; 1 line per other molecule in the bound state\n");
    fprintf(acceptanceRatioFile, "# Iteration AcceptanceRatio\n");
    fprintf(exchangeFrequencyFile, "# Iteration continous_ratio, inst_ratio\n");

    fflush(fractionBoundFile);
    fflush(boundConformationsFile);
    fflush(acceptanceRatioFile);
    fflush(exchangeFrequencyFile);

    fprintf(fractionBoundFile, "Iteration:  ");
    fprintf(acceptanceRatioFile, "Iteration:  ");

    for(map<float, int>::const_iterator iterator = position_of_temperature.begin(); iterator != position_of_temperature.end(); ++iterator)
    {
        float temperature = iterator->first;
        fprintf(fractionBoundFile,  "%0.1fKi %0.1fKc ",temperature, temperature);
        fprintf(acceptanceRatioFile,"%0.1fKi %0.1fKc ",temperature, temperature);
    }

    fprintf(fractionBoundFile," \n");
    fprintf(acceptanceRatioFile,"\n");
}

void Simulation::closeSamplingFiles()
{
    closeSamplingFile("fractionBound", &fractionBoundFile);
    closeSamplingFile("boundconformations", &boundConformationsFile);
    closeSamplingFile("acceptance_ratios", &acceptanceRatioFile);
    closeSamplingFile("exchange_freq", &exchangeFrequencyFile);
}

void Simulation::printHelp()
{
    cout << "Usage: cgppd -f <filename> [-c] [-h] [-v] [-q] [-t x] [-s x] [-g x] [-m x ] [-e x] [-r x] [-o x] [-b x] [-n x] [-x x] [-d x]"<< endl;
    cout << "\t-h|--help: show this dialog" << endl;
    cout << "\t-f|--file <file>: Input config file" << endl;
    cout << "\t-v|--view:        use the open GL preview of this configuration, performs no simulation" << endl;
    cout << "\t-q|--nosim:       Do everything except the simulation (for use with -v)" << endl;
    cout << "The following values override those in the config file" << endl;
    cout << "\t-t|--threads x:  The number of CPU/pthreads to use" << endl;
    cout << "\t-s|--streams x:  The number of CUDA Streams" << endl;
    cout << "\t-g|--gpus x:     The number of GPUS to use" << endl;
    cout << "\t-m|--mcsteps x:  How many MC Steps to perform per replica " << endl;
    cout << "\t-e|--resteps x:  How many replica exhanges to perform" << endl;
    cout << "\t-r|--replicas x: How many replicas" << endl;
    cout << "\t-o|--output x:   The output prefix for files created by the simulation" << endl;
    cout << "\t-b|--boundary x:  The bounding box edge length" << endl;
    cout << "\t-x|--tmax x:    The temperature of the highest replica" << endl;
    cout << "\t-n|--tmin x:    The temperature of the lowest replica" << endl;
    cout << "\t-d|--blockdim x: Number of threads per CUDA block" << endl;

    exit(0);
}

void Simulation::getArgs(int argc, char **argv)
{

    if (argc <= 1)
    {
        cout << "No arguments provided." << endl;
        printHelp();
    }

    sprintf(parameters.logfile,"output/%d_logfile",parameters.pid);

    const struct option long_options[] =
    {
        {"help", no_argument, 0, 'h'},
        {"file", required_argument, 0, 'f'},
        {"view", no_argument, 0, 'v'},
        {"nosim", no_argument, 0, 'q'},

        {"threads", required_argument, 0, 't'},
        {"streams", required_argument, 0, 's'},
        {"gpus", required_argument, 0, 'g'},
        {"mcsteps", required_argument, 0, 'm'},
        {"resteps", required_argument, 0, 'e'},
        {"replicas", required_argument, 0, 'r'},

        {"output", required_argument, 0, 'o'},
        {"boundary", required_argument, 0, 'b'},
        {"tmax", required_argument, 0, 'x'},
        {"tmin", required_argument, 0, 'n'},
        {"blockdim", required_argument, 0, 'd'},

        {0, 0, 0, 0},
    };

    int opt_index = 0;

    // Important! This is a getopt variable which we need to reset if we want to do a second pass over the same options. Which we do.
    optind = 1;

    while (1)
    {
        int opt = getopt_long(argc, argv, "hf:vqt:s:g:m:e:r:o:b:x:n:d:", long_options, &opt_index);

        if (opt == -1)
            break;

        switch(opt)
        {
            case 'h':
                printHelp();
                break;
            case 'f':
                strcpy(parameters.file, optarg);
                parameters.inputFile = true;
                break;
            case 'v':
#if GLVIS
                parameters.viewConditions = true;
#else
                cout << "This build does not support OpenGL." << endl;
#endif
                break;
            case 'q':
                parameters.skipsimulation = true;
                break;
            case 't':
                parameters.threads = atoi(optarg);
                break;
            case 's':
                parameters.streams = atoi(optarg);
                break;
            case 'g':
                parameters.gpus = atoi(optarg);
                break;
            case 'm':
                parameters.MCsteps = atoi(optarg);
                break;
            case 'e':
                parameters.REsteps = atoi(optarg);
                break;
            case 'r':
                parameters.replicas = atoi(optarg);
                break;
            case 'o':
                strcpy(parameters.logfile, optarg);
                break;
            case 'b':
                parameters.bound = atof(optarg);
                break;
            case 'x':
                parameters.temperatureMax = atof(optarg);
                break;
            case 'n':
                parameters.temperatureMin = atof(optarg);
                break;
            case 'd':
#if USING_CUDA
                parameters.cuda_blockSize = atoi(optarg);
                parameters.auto_blockdim = false;
                cout << "Block size changed to: " << parameters.cuda_blockSize << endl;
#else
                cout << "This build does not support CUDA." << endl;
#endif
                break;
            default:
                cout << "Unknown option." << endl;
                printHelp();
                break;
        }
    }

// TODO: put this back when we have checkpointing again
//         else if (strcmp(argv[i],"-z")==0 || strcmp(argv[i],"--resume")==0)
//         {
//             cout << "!!! Checkpointing is not enabled in this build, cannot resume" << endl;
//             i++;
//         }
}

void Simulation::loadArgsFromFile()
{
    /* use an input file
     * files start after the "files" word in a file.
     * format is: px py pz rx ry rz r filename
     * meaning: position 3 floats x,y,z, rotation axis floats rx,ry,rz, and amount r float, filename is the rest of the line
     *
     * file can define parameters too.
     * one per line args:
     * threads 2  (2 threads)
     * streams 4  (4 streams) etc....
     * files
     * p 1 4 7 r 1 0 0 0 pdb data/conf1/1a.pdb
     * crowders
     * p 1 8 9 r 1 9 2 0.5 pdb data/conf1/x.pdb
     *
     */

    if (!parameters.inputFile)
    {// deprecated function
        cout << "No configuration file provided." << endl;
        printHelp();
    }

    ifstream input(parameters.file);

    if (!input.good())
    {
        cout << "Failed to open file: " << parameters.file << endl;
        exit(0);
    }

    char line[512] = {0};

#define PARAMETER_SECTION 0
#define MOLECULE_SECTION 1
#define CROWDER_SECTION 2

    int section = PARAMETER_SECTION;

    while (!input.eof())
    {
        input.getline(line, 512);

        if (line == NULL || line[0] == '#' || strlen(line)==0)
        {
        }
        else if (strcmp(line, "files") == 0)
        {
            section = MOLECULE_SECTION;
        }
        else if (strcmp(line, "crowders") == 0)
        {
            section = CROWDER_SECTION;
        }
        else if (section == MOLECULE_SECTION || section == CROWDER_SECTION)
        {
            moldata m;
            int result = 0;

            if (line[0]=='t')
            {
                result = sscanf(line, "t(%f,%f,%f) r(%f,%f,%f,%f) %s", &m.px, &m.py, &m.pz, &m.rx, &m.ry, &m.rz, &m.ra, m.pdbfilename);
                m.translate = true;
            }

            if (line[0]=='p')
            {
                result = sscanf(line, "p(%f,%f,%f) r(%f,%f,%f,%f) %s", &m.px, &m.py, &m.pz, &m.rx, &m.ry, &m.rz, &m.ra, m.pdbfilename);
                m.translate = false;
            }

            if (result<8)
            {
                cout << "Unable to parse molecule: " << line << endl;
            }

            else
            {
                m.crowder = (section == CROWDER_SECTION);
                // TODO: remove this later?
                if (section == MOLECULE_SECTION)
                {
                    parameters.nonCrowders++;
                }

                parameters.mdata.push_back(m);
            }
        }
        else if (section == PARAMETER_SECTION)
        {
            char * key = strtok(line," ");
            char * value = strtok(NULL," ");

            if (strcmp(key, "gpus") == 0)
            {
                parameters.gpus = atoi(value);
            }
            else if (strcmp(key, "streams") == 0)
            {
                parameters.streams = atoi(value);
            }
            else if (strcmp(key, "threads") == 0)
            {
                parameters.threads = atoi(value);
            }
            else if (strcmp(key, "mcsteps") == 0)
            {
                parameters.MCsteps = atoi(value);
            }
            else if (strcmp(key, "resteps") == 0)
            {
                parameters.REsteps = atoi(value);
            }
            else if (strcmp(key, "temperaturemax") == 0)
            {
                parameters.temperatureMax = atof(value);
            }
            else if (strcmp(key, "temperaturemin") == 0)
            {
                parameters.temperatureMin = atof(value);
            }
            else if (strcmp(key, "boundary") == 0)
            {
                parameters.bound = atof(value);
            }
            else if (strcmp(key, "replicas") == 0)
            {
                parameters.replicas = atoi(value);
            }
            else if (strcmp(key, "samplefrequency") == 0)
            {
                parameters.sampleFrequency = atoi(value);
            }
            else if (strcmp(key, "sampleafter") == 0)
            {
                parameters.sampleStartsAfter = atoi(value);
            }
            else if (strcmp(key, "prefix") == 0)
            {
                strcpy(parameters.prefix, "");
                strcpy(parameters.prefix, value);
#if REPULSIVE_CROWDING
                strcat(parameters.prefix, "_repulsive");
#else
                strcat(parameters.prefix, "_full");
#endif
            }
            else {
                cout << "Unknown parameter: "<< key << endl;
            }
        }
    }
    input.close();
}

void Simulation::check_and_modify_parameters()
{
    if (parameters.bound <= 0)
    {
        cout << "! Bounding value too small; setting equal to " << BOUNDING_VALUE << endl;
        parameters.bound = BOUNDING_VALUE;
    }

    if (parameters.temperatureMax < parameters.temperatureMin)
    {
        cout << "! Maximum temperature < minimum temperature; swapping " << parameters.temperatureMax << " <-> " << parameters.temperatureMin << endl;
        float tmp = parameters.temperatureMax;
        parameters.temperatureMax = parameters.temperatureMin;
        parameters.temperatureMin = tmp;
    }

    if (parameters.threads > parameters.replicas)
    {
        parameters.threads = parameters.replicas;
        cout << "! Threads > replicas; setting threads equal to " << parameters.threads << endl;
    }

    parameters.max_replicas_per_thread = int(ceil(float(parameters.replicas) / float(parameters.threads)));
    int spaces = parameters.replicas - parameters.max_replicas_per_thread * parameters.threads;
    int unused_threads = spaces / parameters.max_replicas_per_thread; // integer division
    if (unused_threads)
    {
        cout << "! After assignment of " << parameters.replicas << " replicas to " << parameters.threads << " threads with " << parameters.max_replicas_per_thread << " replicas per thread, " << unused_threads << " threads are left unused. You must either increase the number of replicas or decrease the number of threads. Exiting.";
        exit(0);
    }

#if USING_CUDA
    int availableGpus;
    cudaGetDeviceCount(&availableGpus);
    if (parameters.gpus > availableGpus)
    {
        parameters.gpus = availableGpus;
        cout << "! Too many GPUs; setting equal to " << parameters.gpus << endl;
    }
#endif

    if (parameters.threads > parameters.streams)
    {
        parameters.streams = parameters.threads;
        if (parameters.streams > 16 * parameters.gpus)
        {
            parameters.streams = 16 * parameters.gpus;

            if (parameters.streams > parameters.replicas)
            {
                parameters.streams = parameters.replicas;
            }
            cout << "! Too many streams; setting equal to " << parameters.streams << endl;
        }
    }

}

void Simulation::writeFileIndex()
{
    char fileindex[256];
    sprintf(fileindex,"output/%s_%d_fileindex", parameters.prefix, parameters.pid);
    FILE * fileindexf = fopen (fileindex,"w");
    fprintf(fileindexf,"index molecule_file_path crowder(Y/N)\n");

    for (int i = 0; i < parameters.mdata.size(); i++) {
        fprintf(fileindexf, "%2d %s %s", i, parameters.mdata[i].pdbfilename, (parameters.mdata[i].crowder ? "Y\n" : "N\n"));
    }

    fclose(fileindexf);
}

void Simulation::printArgs()
{
    cout << "Argument data from file:" << endl;
    cout << "-------------------------------------" << endl;
    cout << "threads " << parameters.threads << endl;
    cout << "streams " << parameters.streams << endl;
    cout << "GPUs " << parameters.gpus << endl;
    cout << "mc steps " << parameters.MCsteps << endl;
    cout << "re steps " << parameters.REsteps << endl;
    cout << "replicas " << parameters.replicas << endl;
    cout << "sampling frequency (mc steps) " << parameters.sampleFrequency << endl;
    cout << "sampling starts after (mc steps) " << parameters.sampleStartsAfter << endl;
    cout << "bounding box size " << parameters.bound << endl;
    cout << "non-crowder molecules " << parameters.nonCrowders << endl;
    cout << "maximum temperature " << parameters.temperatureMax << endl;
    cout << "minimum temperature " << parameters.temperatureMin << endl;

    cout << "Loaded: "<< endl;
    cout << "-------------------------------------------------------------"<< endl;


    for (int i = 0; i < parameters.mdata.size(); i++) {
        printf("%2d %s%s centered @ (%0.3f,%0.3f,%0.3f)\n", i, parameters.mdata[i].pdbfilename, (parameters.mdata[i].crowder ? " crowder" : ""), parameters.mdata[i].px, parameters.mdata[i].py, parameters.mdata[i].pz);
    }

    cout << "-------------------------------------------------------------"<< endl;
}
