#include "Simulation.h"

Simulation::Simulation()
{
    replicasInitialised = false;
    REMCRng = gsl_rng_alloc(gsl_rng_mt19937);

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
    if (replicasInitialised)
    {
        for (size_t i=0; i<parameters.replicas; i++)
        {
            // TODO: WTF, why are these not in the replica destructor?!
            replica[i].freeRNGs();
//             delete [] replica[i].contiguousResidues; // this is now
        }
    }
}

void Simulation::init(int argc, char **argv, int pid)
{
    int sysreturn;
    sysreturn = system("mkdir -p output output_pdb checkpoints");

    parameters.pid = pid;

    getArgs(argc, argv);
    loadArgsFromFile();
    check_and_modify_parameters();
    writeFileIndex();
    getArgs(argc, argv); // second pass to override any variables if doing performance tests

    // TODO: remove magic number; make initial array size a constant
    initialReplica.init_first_replica(parameters.mdata, aminoAcidData, parameters.bound, 30);
}

void Simulation::calibrate()
{
    cout << "Loaded: " << initialReplica.residueCount << " residues in " << initialReplica.moleculeCount << " molecules:\n";

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
    initialReplica.potential = p;
    cout << initialReplica.countpairs() << " residue interaction pairs." << endl;
    printf("CPU initial energy value: %12.8f (%12.8f) kcal/mol\n", p, pnc);


#if USING_CUDA
    // set box size
    if (parameters.auto_blockdim)
        (initialReplica.residueCount < 1024) ? parameters.cuda_blockSize = 32 : parameters.cuda_blockSize = 64;

    initialReplica.setBlockSize(parameters.cuda_blockSize);

    // set box dimensions
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

    cout << "Last CUDA error: " << cudaGetErrorString(cudaGetLastError()) << endl;
    cout.flush();
}

pthread_mutex_t waitingThreadMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t waitingCounterMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t reMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t writeFileMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t waitingThreadCond;
pthread_cond_t waitingReplicaExchangeCond;

void Simulation::run()
{
    LOG(ALWAYS, "Beginning simulation\n");

    // need at most REPLICA_COUNT threads
    cout << "Output files will be prefixed by " << parameters.prefix << "_" << parameters.pid << endl;

    if (!parameters.resume)
        // TODO: remove this; it's not actually used (and is a dupliacte of steps)
        parameters.currentStep = 0;
    // TODO: some of these should be attributes; some of them are pointless duplication of parameters.
    // TODO: also, all this parameter fixing should go in a different function.

    int waitingThreads = 0;
    int conformationsBound = 0;

    if (initialReplica.moleculeCount == 0)  // make sure something is loaded
    {
        cout << "-!- No molecules loaded. -!-" << endl;
        cout << "-!- Aborting run. -!- " << endl;
        return;
    }

    // copy the data into each replica
    // copy initial replica to other replicas and init the rngs for each
    // we can't use pthreads and CUDA at the moment, but we are going to use streams

    //basically we change the simulation to spawn N threads N = number of CPU cores

    //initialise the replicas data
    int _300kReplica = 0;

    double geometricTemperature = pow(double(parameters.temperatureMax/parameters.temperatureMin),double(1.0/double(parameters.replicas-1)));
    double geometricTranslate = pow(double(MAX_TRANSLATION/MIN_TRANSLATION),double(1.0/double(parameters.replicas-1)));
    double geometricRotation = pow(double(MAX_ROTATION/MIN_ROTATION),double(1.0/double(parameters.replicas-1)));

    for (size_t i=0; i<parameters.replicas; i++)
    {
        // changed to geometric sequence
        float temperature = parameters.temperatureMin * pow(geometricTemperature, int(i));
        float rotate_step = MIN_ROTATION * pow(geometricRotation, int(i));
        float translate_step = MIN_TRANSLATION * pow(geometricTranslate, int(i));

        // note which replica is the 300K replica for sampling
        // if there is no replica at 300K then choose the closest one.
        if ( abs(replica[i].temperature-300.0f) < abs(replica[_300kReplica].temperature-300.0f))
            _300kReplica = i;

        replica[i].init_child_replica(initialReplica, i, temperature, rotate_step, translate_step, parameters.threads);
        printf ("Replica %d %.3f %.3f %.3f\n",int(i),replica[i].temperature,replica[i].translateStep,replica[i].rotateStep);

    }

    replicasInitialised = true;

    initSamplingFiles();


    // round the ~300K replica to 300K to ensure there is at least one T/T_i == 1
    //replica[_300kReplica].temperature = 300.0f;

    // create a lookup map for replica temperatures to allow for sampling temperatures and in place exchanges
    // samples are by temperature not replica index, hence temperatures get shuffled about in replica[],
    // so use replicaTemperaturesmap[] to track where a specific temperature is

    typedef map<float, int> KRmap;
    KRmap replicaTemperatureMap;
    typedef map<int,float> RKmap;
    RKmap temperatureMap;


    // build a map of lookups for temperature lookups while sampling
    for (size_t i=0; i<parameters.replicas; i++)
    {
        // lookup of where replica with temperature x is
        replicaTemperatureMap[replica[i].temperature] = i;
        // lookup of which replica has temperature x
        temperatureMap[i] = replica[i].temperature;
    }

    fprintf(fractionBoundFile,"Iteration:  ");
    fprintf(acceptanceRatioFile,"Iteration:  ");

    for(KRmap::const_iterator iterator = replicaTemperatureMap.begin(); iterator != replicaTemperatureMap.end(); ++iterator)
    {
        float temperature = iterator->first;
        fprintf(fractionBoundFile,  "%0.1fKi %0.1fKc ",temperature,temperature);
        fprintf(acceptanceRatioFile,"%0.1fKi %0.1fKc ",temperature,temperature);
    }

    fprintf(fractionBoundFile," \n");
    fprintf(acceptanceRatioFile,"\n");


    /*#if OUTPUT_LEVEL > 1
    for (int i=0;i<parameters.replicas;i++)
    {
        sprintf (tmps,"%3d: %1.1e",i,replica[i].potential);
        strcat(startEs,tmps);
    }
    cout << "start energies (kcal/mol): " << startEs << endl;
    #endif*/

    pthread_t *thread = new pthread_t[parameters.threads];;
    //For portability, explicitly create threads in a joinable state
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_cond_init (&waitingThreadCond, NULL);
    pthread_cond_init (&waitingReplicaExchangeCond, NULL);
    pthread_mutex_init(&waitingCounterMutex, NULL);
    pthread_mutex_init(&waitingThreadMutex, NULL);
    pthread_mutex_init(&reMutex, NULL);

    SimulationData *data = new SimulationData[parameters.threads];

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(RELoopTimer) );
#endif

    printf ("--- Starting simulation ---.\n");


    pthread_mutex_lock(&waitingCounterMutex);

    //TODO: assign streams as a function of the number of GPUs and threads

    //spawn N threads to do the Monte-Carlo mutations
    for (int i=0; i<parameters.threads; i++)
    {
        data[i].replica = replica;
        data[i].replicaCount = parameters.replicas;
        data[i].index = i;
        data[i].threads = parameters.threads;
        data[i].streams = parameters.streams;///parameters.threads; //determines the number of streams available
        data[i].MCsteps = parameters.MCsteps;
        data[i].REsteps = parameters.REsteps;
        data[i].sampleFrequency = parameters.sampleFrequency;
        data[i].sampleStartsAfter = parameters.sampleStartsAfter;
        data[i].bound = parameters.bound;
        data[i].streams = parameters.streams;
        data[i].waitingThreadCount = &waitingThreads;
        data[i].conformationsBound = &conformationsBound;
        data[i].fractionBound = fractionBoundFile;
        data[i].boundConformations = boundConformationsFile;

        // assign gpus in rotation per thread, t0 = gpu0, t1 = gpu1 etc
        // % #gpus so they share if threads > gpus
        // will perform best if threads:gpus = 1:1
        data[i].GPUID = i % parameters.gpus;
        cout << "Assign thread " << i << " to GPU " << data[i].GPUID << endl;

        //start the MC loops
        pthread_create(&thread[i], &attr, MCthreadableFunction, (void*)&data[i]);
        //MCthreadableFunction(&data[i]);
    }

    int exchanges = 0;  // number of exchanges performed
    int tests = 0;      // number of replica exchange tests
    int totalExchanges = 0; // total number of exchanges
    int totalTests = 0; // total number of exchanges

    int offset = 0;  // RE offset
    int steps = 0;   // # re steps performed
    while (++steps < parameters.REsteps)  // until enough steps taken
    {
#if INCLUDE_TIMERS
        CUT_SAFE_CALL( cutStartTimer(MCLoopTimer) );
#endif

        pthread_cond_wait(&waitingReplicaExchangeCond,&waitingCounterMutex);

#if INCLUDE_TIMERS
        CUT_SAFE_CALL( cutStopTimer(MCLoopTimer) );
#endif

        //if the sampling has started
        if (parameters.MCsteps/parameters.REsteps*(steps-1) >= parameters.sampleStartsAfter)
        {
            // do fraction bound and acceptance ratio

            fprintf(fractionBoundFile,"%9d: ",parameters.MCsteps/parameters.REsteps*steps);
            fprintf(acceptanceRatioFile,"%9d: ",parameters.MCsteps/parameters.REsteps*steps);

            for(RKmap::const_iterator iterator = temperatureMap.begin(); iterator != temperatureMap.end(); ++iterator)
            {
                int index = replicaTemperatureMap[iterator->second];


                // fraction bound
                float fractionBound = float(replica[index].boundSamples)/max(1.0f,float(replica[index].samplesSinceLastExchange));
                replica[index].totalBoundSamples += replica[index].boundSamples;
                replica[index].totalSamples += replica[index].samplesSinceLastExchange;
                replica[index].boundSamples = 0;
                replica[index].samplesSinceLastExchange = 0;
                float accumulativeFractionBound = float(replica[index].totalBoundSamples)/max(1.0f,float(replica[index].totalSamples));
                fprintf(fractionBoundFile,"| %6.4f %6.4f ",fractionBound,accumulativeFractionBound);

                // acceptance rate
                float acceptanceRatio =  float(replica[index].acceptA + replica[index].accept)/float(replica[index].acceptA + replica[index].accept + replica[index].reject);
                replica[index].totalAccept += replica[index].acceptA+replica[index].accept ;
                replica[index].totalAcceptReject += replica[index].acceptA+replica[index].accept+replica[index].reject ;
                fprintf(acceptanceRatioFile,"| %6.4f %6.4f ",acceptanceRatio , float(replica[index].totalAccept)/float(replica[index].totalAcceptReject));
                replica[index].acceptA=0 ;
                replica[index].accept=0 ;
                replica[index].reject=0 ;
            }

            fprintf(fractionBoundFile,"\n");
            fflush(fractionBoundFile);
            fprintf(acceptanceRatioFile,"\n");
            fflush(acceptanceRatioFile);
        }

        int i = offset;  // replica exchange offset
        while (i < parameters.replicas-1) // loop through the replicas in pairs
        {
            int j = i+1; // replica to swap with current replica

            // i and j represent temperatures, hence we need to map temperature to position
            // z = temperatureMap[x] -> temperature of position i
            // replicaTemperatureMap[z] -> where the replica with temperature z actually is

            int t_i = replicaTemperatureMap[temperatureMap[i]];
            int t_j = replicaTemperatureMap[temperatureMap[j]];


            double delta = (1.0/replica[t_i].temperature - 1.0/replica[t_j].temperature)*(replica[t_i].potential - replica[t_j].potential)*(4184.0f/Rgas);
            if (gsl_rng_uniform(REMCRng) < min(1.0,exp(delta)))
            {
                replica[t_i].exchangeReplicas(replica[t_j]);

                replicaTemperatureMap[replica[t_i].temperature] = t_i;
                replicaTemperatureMap[replica[t_j].temperature] = t_j;

                exchanges++; // sampling
            }
            i += 2; // compare the next two neighbours
            tests++;// tests performed samplng
        }
        offset = 1-offset; // switch the swap neighbours in the replica exchange

        parameters.currentStep = mcstepsPerRE*steps;

        waitingThreads = 0;
        pthread_cond_broadcast(&waitingThreadCond);

#if GLVIS
        GLreplica = &replica[_300kReplica];
        GlutDisplay();
#endif

        int tempI = replicaTemperatureMap[300.0f];
        float frac = float(replica[tempI].totalBoundSamples)/max(1.0f,float(replica[tempI].totalSamples));
        LOG(ALWAYS, "Replica Exchange step %d of %d complete (Fraction bound @ 300K: %f)\n", steps, parameters.REsteps, frac);

        totalExchanges += exchanges;
        totalTests += tests;

        fprintf(exchangeFrequencyFile,"%10d %6.4f %6.4f\n",steps, float(totalExchanges)/float(totalTests), float(exchanges)/float(tests));
        tests = 0;
        exchanges = 0;

    }

    //printf ("Replica Exchange step %d skipped as it has no effect\n",steps);

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
    fprintf(fractionBoundFile,"%9d: ",parameters.MCsteps);
    fprintf(acceptanceRatioFile,"%9d: ",parameters.MCsteps);

    //final fraction bound
    for(RKmap::const_iterator iterator = temperatureMap.begin(); iterator != temperatureMap.end(); ++iterator)
    {
        int index = replicaTemperatureMap[iterator->second];

        // fraction bound
        float fractionBound = float(replica[index].boundSamples)/float(replica[index].samplesSinceLastExchange);
        replica[index].totalBoundSamples += replica[index].boundSamples;
        replica[index].totalSamples += replica[index].samplesSinceLastExchange;
        replica[index].boundSamples = 0;
        replica[index].samplesSinceLastExchange = 0;
        float accumulativeFractionBound = float(replica[index].totalBoundSamples)/float(replica[index].totalSamples);
        fprintf(fractionBoundFile,"| %6.4f %6.4f  ",fractionBound,accumulativeFractionBound);

        // acceptance rate
        float acceptanceRatio =  float(replica[index].acceptA)/float(replica[index].acceptA + replica[index].accept + replica[index].reject);
        fprintf(acceptanceRatioFile,"| %6.4f %6.4f ",acceptanceRatio , float(replica[index].totalAccept)/float(replica[index].totalAcceptReject));


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
    cutDeleteTimer(RELoopTimer);
    cutDeleteTimer(MCLoopTimer);
#endif

    pthread_mutex_lock(&writeFileMutex);
    closeSamplingFiles();
    pthread_mutex_unlock(&writeFileMutex);


    // Clean up/
    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&waitingCounterMutex);
    pthread_mutex_destroy(&waitingThreadMutex);
    pthread_cond_destroy(&waitingThreadCond);
    pthread_cond_destroy(&waitingReplicaExchangeCond);
    gsl_rng_free(REMCRng);
    delete [] data;
    delete [] thread;

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

void *MCthreadableFunction(void *arg)
{

    SimulationData *data = (SimulationData *) arg;
    long threadIndex = data->index;
    Replica *replica = data->replica;

    printf ("--- Monte-Carlo thread %d running. ---.\n",int(threadIndex+1));

#if GLVIS
    GLreplica = &replica[0];
#endif
    // determine the number of replicas the thread will run
    int tReplicas = int(ceil(float(data->replicaCount)/float(data->threads)));
    //to make sure we have all replicas execute, check that tReplicas*(i+1) < REPLICA_COUNT
    int replicasInThisThread = tReplicas;
    if  (tReplicas*(threadIndex+1)>data->replicaCount)
        replicasInThisThread = tReplicas - (tReplicas*data->threads)%data->replicaCount;

    if (replicasInThisThread-1+threadIndex*tReplicas > data->replicaCount)
    {
        cout << "TOO MANY THREADS. Either increase the number of replicas or decrease the thread count. " << endl;
        exit(0);
    }

#ifdef VERY_VERBOSE
    pthread_mutex_lock( &replicaMutex );
    for (int tx=0; tx<replicasInThisThread; tx++)
    {
        cout << "+ Thread " << threadIndex << " running replica " << tx+threadIndex*tReplicas << endl;
    }
    pthread_mutex_unlock( &replicaMutex );
#endif

#if USING_CUDA

    //initialise cuda for use in this thread
    cuInit(0);
    cutilCheckMsg("Failed to initialise CUDA runtime.");

    cudaSetDevice(data->GPUID);
    cutilCheckMsg("Failed to pick device for the CUDA runtime.");

    // copy the LJpotentials to gpu memory in this thread context so it can access it
    float *deviceLJpTmp;
    cudaMalloc((void**)&deviceLJpTmp,sizeof(float)*AA_COUNT*AA_COUNT);
    cutilCheckMsg("Failed to allocate contact potential memory on the GPU");

    CUDA_setBoxDimension(replica->boundingValue);
    cutilCheckMsg("Failed to copy box dimensions to GPU");

    copyLJPotentialDataToDevice(deviceLJpTmp, &replica[threadIndex*tReplicas].aminoAcids);
    cutilCheckMsg("Failed to code contact potentials to device.");


    // create streams for each subsequent MCSearch
    // ensure the replica can find the lookup table
    // initialise data on the device and copy the initial batch

#if CUDA_STREAMS

    cudaStream_t streams[16];   // reserve 16 stream slots but create only as many as needed
    // the stream/replica ration must be a whole number otherwise there will be lots of waste, ie dormant streams etc
    int streamsAvailablePerThread = data->streams/data->threads;
    int replicasPerStream = int(ceil(float(replicasInThisThread)/float(streamsAvailablePerThread)));
    int streamsPerThread  = replicasInThisThread/replicasPerStream;

    for (int is=0; is<streamsPerThread; is++)
    {
        cudaStreamCreate(&streams[is]);
    }
#endif

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    bindLJTexture(deviceLJpTmp);
#endif

    for (int tx=0; tx<replicasInThisThread; tx++)
    {
        replica[tx+threadIndex*tReplicas].device_LJPotentials = deviceLJpTmp;
        replica[tx+threadIndex*tReplicas].ReplicaDataToDevice();


        //printf("Initial potential value for replica %d (compute thread %d): %20.10f\n",int(tx+threadIndex*tReplicas) ,int(threadIndex), replica[tx+threadIndex*tReplicas].EonDevice());
#if CUDA_STREAMS
        // TODO: does this access the same replica?! I assume so.
        data->replica[tx+threadIndex*tReplicas].cudaStream = streams[tx%streamsPerThread];  // use rotation to assign replicas to streams
        replica[tx+threadIndex*tReplicas].ReserveSumSpace();                                // reserve a space for the potential summation to be stored
//         replica[tx+threadIndex*tReplicas].savedMolecule.reserveResidueSpace(replica[tx+threadIndex*tReplicas].maxMoleculeSize); // TODO: moved this to one of the init functions
#endif

    }



#endif

    int mcstepsPerRE = data->MCsteps/data->REsteps; // number of MC steps to do at a time

    // run the MC loop
#if !CUDA_MC && !CUDA_STREAMS   // if no MC on device
    int mcstep = 0;
    int samples = 0;
    while (mcstep < data->MCsteps)
    {
        //cout << "Starting mc loop at " << mcstep << " steps" << endl;
        // do all MC steps before RE
        for (int tx=0; tx<replicasInThisThread; tx++)
        {
            // to sample at the correct rate split this into another set of loops
            for (int s=0; s<mcstepsPerRE/data->sampleFrequency; s++)
            {
                replica[tx+threadIndex*tReplicas].MCSearch(data->sampleFrequency);

                if (mcstep + s*data->sampleFrequency >= data->sampleStartsAfter)
                {
                    // cout << "Sampling at step:" << mcstep+(s+1)*data->sampleFrequency << endl;
                    // if E < -1.1844 kcal/mol then its bound
                    replica[tx+threadIndex*tReplicas].sample(data->boundConformations,mcstep+s*data->sampleFrequency,BOUND_ENERGY_VALUE,&writeFileMutex);
                }
            }
        }
        mcstep += mcstepsPerRE;
        //cout << "thread " << threadIndex << " waits for mutex " << endl;

        // do replica exchange
        pthread_mutex_lock(&waitingCounterMutex);                 // lock the counter
        //cout << "thread " << threadIndex << " gets mutex " << endl;

        data->waitingThreadCount[0] = data->waitingThreadCount[0]+1;
        if (data->waitingThreadCount[0]==data->threads)                 // if all threads in waiting state
        {
            pthread_cond_signal(&waitingReplicaExchangeCond);
            //cout << "thread " << threadIndex << " signals parent" << endl;
        }
        if (mcstep<data->MCsteps)                       // wait if another MC loop must be done
        {
            //cout << "thread " << threadIndex << " waits after " << mcstep << endl;

            pthread_cond_wait(&waitingThreadCond,&waitingCounterMutex); // wait for replica exchange to finish.
            // NB! unlock the mutex locked upon the condition above
            // being met because this will unblock all threads such
            // that they will continue concurrently, if not unlocked
            // other threads will run sequentially
            //cout << "thread " << threadIndex << " releases mutex after wait" << endl;

            //  pthread_mutex_unlock(&waitingCounterMutex);
        }
        pthread_mutex_unlock(&waitingCounterMutex);


    } // continue MC
#endif


#if CUDA_STREAMS    // use streams
    // TODO: move this to Replica?
    int mcstep = 0;
    //int sampleIn = data->sampleStartsAfter+1;
    while (mcstep<data->MCsteps)
    {
        // to sample at the correct rate split this into another set of loops
        for (int mcx=0; mcx<mcstepsPerRE; mcx++) // at each mc step
        {
            for (int index=0; index<replicasInThisThread; index+=replicasPerStream)
            {
                for (int rps=0; rps<replicasPerStream; rps++)
                {
                    // batch replicas such that no stream is shared per batch
                    replica[threadIndex*tReplicas+index+rps].MCSearchMutate();
                    replica[threadIndex*tReplicas+index+rps].MCSearchEvaluate();

                }
                for (int rps=0; rps<replicasPerStream; rps++)
                {
                    replica[threadIndex*tReplicas+index+rps].MCSearchAcceptReject();
                }

                for (int rps=0; rps<replicasPerStream; rps++)
                {
                    //sampleAsync
                    if (mcstep%data->sampleFrequency == 0 && mcstep >= data->sampleStartsAfter) // when enough steps are taken && sampleFrequency steps have passed
                    {
                        replica[threadIndex*tReplicas+index+rps].sample(data->boundConformations,mcstep+mcx*data->sampleFrequency,BOUND_ENERGY_VALUE,&writeFileMutex);

                        //if (abs(replica[threadIndex*tReplicas+index+rps].temperature-300.0f)<1.0f)
                        //  cout << "Sampled: t=" << replica[threadIndex*tReplicas+index+rps].temperature << " mcstep=" << mcstep << endl;
                    }

                }
            }
            mcstep++;
        }

        // do replica exchange
        pthread_mutex_lock(&waitingCounterMutex);                 // lock the counter
        data->waitingThreadCount[0] = data->waitingThreadCount[0]+1;
        if (data->waitingThreadCount[0]==data->threads)
        {
            pthread_cond_signal(&waitingReplicaExchangeCond);
        }
        if (mcstep<data->MCsteps)                                           // wait if another MC loop must be done
        {
            pthread_cond_wait(&waitingThreadCond,&waitingCounterMutex);     // wait for replica exchange to finish.
        }
        pthread_mutex_unlock(&waitingCounterMutex);

    } // continue MC
#endif



#if USING_CUDA
    // sync all streams and free gpu memory
    for (int tx=0; tx<replicasInThisThread; tx++)
    {
#if CUDA_STREAMS
        replica[tx+threadIndex*tReplicas].FreeSumSpace();
//         delete [] replica[tx+threadIndex*tReplicas].savedMolecule.Residues;
#endif
        replica[tx+threadIndex*tReplicas].FreeDevice();
    }

#if CUDA_STREAMS
    for (int i=0; i<streamsPerThread; i++)
    {
        cudaStreamDestroy(streams[i]);
    }
#endif

    // free the memory used by the LJ table
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    unbindLJTexture();
#endif
    cudaFree(deviceLJpTmp);
#endif //using cuda
    /*pthread_mutex_lock(&waitingCounterMutex);                     // lock the counter
    data->waitingThreadCount[0] = data->waitingThreadCount[0]+1;        // wait one last time for all threads to complete before joining them.
    if (data->waitingThreadCount[0]==data->threads)
    {
        pthread_cond_signal(&waitingReplicaExchangeCond);
    }
    pthread_mutex_unlock(&waitingCounterMutex);             //unlock the counter
    */
    printf (" >>> Monte-Carlo thread %d exited.\n",int(threadIndex));
    return 0;
}
// END OF *MCthreadableFunction

// TODO: can we eliminate the copypasta? Loop?
void Simulation::initSamplingFiles()
{
    char filename[256];

    sprintf(filename,"output/%s_%d_fractionBound",parameters.prefix,parameters.pid);
    fractionBoundFile = fopen (filename,"at"); // attempt append a file of the same name
    if (!fractionBoundFile) fractionBoundFile = fopen(filename, "wt"); // create if that failed, and open and append failed (eg: file permissions r--)
    if (!fractionBoundFile)
    {
        printf("Cannot open/create file: %s to record fraction bound.\n",filename);
        return;
    }
    fprintf(fractionBoundFile,"# Fraction Bound\n#Iteration InstantaniousAve CumulativeAve\n");

    sprintf(filename,"output/%s_%d_boundconformations",parameters.prefix,parameters.pid);
    boundConformationsFile = fopen (filename,"a+");
    if (!boundConformationsFile) boundConformationsFile = fopen(filename, "wt");
    if (!boundConformationsFile)
    {
        printf("Cannot open/create file: %s to record bound conformations.\n",filename);
        return;
    }
    fprintf(boundConformationsFile,"# iteration; complex free energy (replica free energy); temperature;\n# molecule: rotation(w,x,y,z) position(x,y,z) 1 line per other molecule in the bound state\n");

    sprintf(filename,"output/%s_%d_acceptance_ratios",parameters.prefix,parameters.pid);
    acceptanceRatioFile = fopen (filename,"a+");
    if (!acceptanceRatioFile) acceptanceRatioFile = fopen(filename, "wt");
    if (!acceptanceRatioFile)
    {
        printf("Cannot open/create file: %s to record acceptance ratio.\n",filename);
        return;
    }
    fprintf(acceptanceRatioFile,"# Iteration AcceptanceRatio\n");

    sprintf(filename,"output/%s_%d_exchange_freq",parameters.prefix,parameters.pid);
    exchangeFrequencyFile = fopen (filename,"a+");
    if (!exchangeFrequencyFile) exchangeFrequencyFile = fopen(filename, "wt");
    if (!exchangeFrequencyFile)
    {
        printf("Cannot open/create file: %s to record exchanges.\n",filename);
        return;
    }
    fprintf(exchangeFrequencyFile,"# Iteration continous_ratio, inst_ratio\n");

    fflush(acceptanceRatioFile);
    fflush(fractionBoundFile);
    fflush(boundConformationsFile);
    fflush(exchangeFrequencyFile);

    return;
}

void Simulation::closeSamplingFiles()
{
    fclose(fractionBoundFile);
    fclose(boundConformationsFile);
    fclose(acceptanceRatioFile);
    fclose(exchangeFrequencyFile);

    LOG(ALWAYS, ">>> Sampling files closed.\n");
    LOG(ALWAYS, "    - output/%s_%d_fractionBound\n", parameters.prefix, parameters.pid);
    LOG(ALWAYS, "    - output/%s_%d_boundconformations\n", parameters.prefix, parameters.pid);
    LOG(ALWAYS, "    - output/%s_%d_acceptance_ratios\n", parameters.prefix, parameters.pid);
    LOG(ALWAYS, "    - output/%s_%d_exchange_freq\n", parameters.prefix, parameters.pid);
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
        cout << "! Bounding value too small, setting equal to " << BOUNDING_VALUE << endl;
        parameters.bound = BOUNDING_VALUE;
    }

    if (parameters.temperatureMax < parameters.temperatureMin)
    {
        cout << "! Maximum temperature < minimum temperature, swapping " << parameters.temperatureMax << " <-> " << parameters.temperatureMin << endl;
        float tmp = parameters.temperatureMax;
        parameters.temperatureMax = parameters.temperatureMin;
        parameters.temperatureMin = tmp;
    }

    if (parameters.threads > parameters.replicas)
    {
        parameters.threads = parameters.replicas;
        cout << "! Too many threads, setting equal to " << parameters.threads << endl;
    }

#if USING_CUDA
    int availableGpus;
    cudaGetDeviceCount(&availableGpus);
    if (parameters.gpus > availableGpus)
    {
        parameters.gpus = availableGpus;
        cout << "! Too many GPUs, setting equal to " << parameters.gpus << endl;
    }
#endif

    if (parameters.threads > parameters.streams)
    {
        parameters.streams = parameters.threads;
        if (parameters.streams > 16*parameters.gpus)
        {
            parameters.streams = 16*parameters.gpus;
            if (parameters.streams > parameters.replicas)
                parameters.streams = parameters.replicas;
            cout << "! Too many streams, setting equal to " << parameters.streams << endl;
        }
    }

    mcstepsPerRE = parameters.MCsteps/parameters.REsteps; // number of MC steps to do at a time

    // if samping does not divide into steps in this loop
    if (mcstepsPerRE % parameters.sampleFrequency > 0)// && mcstepsPerRE > data->sampleFrequency)
    {
        if (mcstepsPerRE > parameters.sampleFrequency)
        {
            // fit sampling frequency such that its a divisor or equal to the mc steps we do before RE
            while (mcstepsPerRE % parameters.sampleFrequency != 0)
                parameters.sampleFrequency++;
            cout << "-!- CHANGED: sample frequency changed as it does divide into the replica exchange frequency. -!-" << endl;
            cout << "-!- CHANGED: sample frequency <- "<< parameters.sampleFrequency <<". -!-" << endl;

        }
        if (mcstepsPerRE < parameters.sampleFrequency)
        {
            parameters.sampleFrequency = mcstepsPerRE;
            cout << "-!- CHANGED: sample frequency too long for MC loop. -!-" << endl;
            cout << "-!- CHANGED: sample frequency <- "<< mcstepsPerRE <<". -!-" << endl;
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