#include "Simulation.h"

Simulation::Simulation() : waitingThreads(0), exchanges(0), tests(0),  totalExchanges(0), totalTests(0), exchangeFraction(0), accumulativeExchangeFraction(0), offset(0), steps(0), thread_created(false), near_300k(0.0f)
{
    srand(time(NULL));
    REMCRng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set (REMCRng, random());

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

void Simulation::printHelp()
{
    cout << "Usage: cgppd -f <filename> [-h] [-p] [-q] [-v] [-t x] [-s x] [-g x] [-m x ] [-a x] [-e x] [-r x] [-o x] [-b x] [-n x] [-x x] [-d x]"<< endl;
    cout << "\t-h|--help: show this dialog" << endl;
    cout << "\t-f|--file <file>: Input config file" << endl;
    cout << "\t-p|--preview:        use the open GL preview of this configuration, performs no simulation" << endl;
    cout << "\t-q|--nosim:       Do everything except the simulation (for use with -p)" << endl;
    cout << "The following values override those in the config file" << endl;
    cout << "\t-t|--threads x:  The number of CPU/pthreads to use" << endl;
    cout << "\t-s|--streams x:  The number of CUDA streams" << endl;
    cout << "\t-g|--gpus x:     The number of GPUS to use" << endl;
    cout << "\t-m|--mcsteps x:  How many MC steps to perform per replica " << endl;
    cout << "\t-a|--sampleafter x:  How many MC steps to perform before beginning sampling" << endl;
    cout << "\t-e|--resteps x:  How many replica exchanges to perform" << endl;
    cout << "\t-r|--replicas x: How many replicas" << endl;
    cout << "\t-o|--output x:   The output prefix for files created by the simulation" << endl;
    cout << "\t-b|--boundary x:  The bounding box edge length" << endl;
    cout << "\t-x|--tmax x:    The temperature of the highest replica" << endl;
    cout << "\t-n|--tmin x:    The temperature of the lowest replica" << endl;
    cout << "\t-d|--blockdim x: Number of threads per CUDA block" << endl;
    cout << "\t-v|--v: Verbosity level (glog flag)" << endl;

    exit(0);
}

void Simulation::getArgs(int argc, char **argv, bool first_pass)
{
    const struct option long_options[] =
    {
        {"file", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"v", required_argument, 0, 'v'},
        
        {"preview", no_argument, 0, 'p'},
        {"nosim", no_argument, 0, 'q'},

        {"threads", required_argument, 0, 't'},
        {"streams", required_argument, 0, 's'},
        {"gpus", required_argument, 0, 'g'},
        {"mcsteps", required_argument, 0, 'm'},
        {"sampleafter", required_argument, 0, 'a'},
        {"resteps", required_argument, 0, 'e'},
        {"replicas", required_argument, 0, 'r'},

        {"boundary", required_argument, 0, 'b'},
        {"tmax", required_argument, 0, 'x'},
        {"tmin", required_argument, 0, 'n'},
        {"blockdim", required_argument, 0, 'd'},

        {0, 0, 0, 0},
    };

    const char short_options[] = "hf:pqt:s:g:m:a:e:r:b:x:n:d:v:w:";

    // Reset getopt's awful global variables
    int opt_index = 0;
    optind = 1;

    if (first_pass) {
        while (1)
        {
            int opt = getopt_long(argc, argv, short_options, long_options, &opt_index);

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
                    // This is a bit of a hack, but I really don't want to refactor everything to use gflags just for verbose logs.
                    FLAGS_v = atoi(optarg); 
                    break;
                default:
                    // Ignore all other options in this pass
                    break;
            }
        }

        if (!parameters.inputFile)
        {
            LOG(ERROR) << "No configuration file provided.";
            printHelp();
        }

        VLOG(1) << "Will read configuration from file: " << parameters.file;
    } else {

        VLOG(0) << "Parsing commandline parameters...";

        while (1)
        {
            int opt = getopt_long(argc, argv, short_options, long_options, &opt_index);

            if (opt == -1)
                break;

            switch(opt)
            {
                case 'h':
                    printHelp();
                    break;
                case 'f':
                    // we can ignore the file parameter on the second pass
                    break;
                case 'v':
                    // ignore
                    break;
                case 'w':
                    // ignore
                    break;
                case 'p':
    #if GLVIS
                    parameters.viewConditions = true;
                    VLOG(1) << "\tWill show OpenGL preview.";
    #else
                    LOG(WARNING) << "\tThis build does not support OpenGL.";
    #endif
                    break;
                case 'q':
                    parameters.skipsimulation = true;
                    VLOG(1) << "\tWill skip simulation.";
                    break;
                case 't':
                    parameters.threads = atoi(optarg);
                    VLOG(1) << "\tParameter threads = " << parameters.threads;
                    break;
                case 's':
                    parameters.streams = atoi(optarg);
                    VLOG(1) << "\tParameter streams = " << parameters.streams;
                    break;
                case 'g':
                    parameters.gpus = atoi(optarg);
                    VLOG(1) << "\tParameter gpus = " << parameters.gpus;
                    break;
                case 'm':
                    parameters.MCsteps = atoi(optarg);
                    VLOG(1) << "\tParameter MCsteps = " << parameters.MCsteps;
                    break;
                case 'a':
                    parameters.sampleStartsAfter = atoi(optarg);
                    VLOG(1) << "\tParameter sampleStartsAfter = " << parameters.sampleStartsAfter;
                    break;
                case 'e':
                    parameters.REsteps = atoi(optarg);
                    VLOG(1) << "\tParameter REsteps = " << parameters.REsteps;
                    break;
                case 'r':
                    parameters.replicas = atoi(optarg);
                    VLOG(1) << "\tParameter replicas = " << parameters.replicas;
                    break;
                case 'b':
                    parameters.bound = atof(optarg);
                    VLOG(1) << "\tParameter bound = " << parameters.bound;
                    break;
                case 'x':
                    parameters.temperatureMax = atof(optarg);
                    VLOG(1) << "\tParameter temperatureMax = " << parameters.temperatureMax;
                    break;
                case 'n':
                    parameters.temperatureMin = atof(optarg);
                    VLOG(1) << "\tParameter temperatureMin = " << parameters.temperatureMin;
                    break;
                case 'd':
    #if USING_CUDA
                    parameters.cuda_blockSize = atoi(optarg);
                    parameters.auto_blockdim = false;
                    VLOG(1) << "\tParameter blockdim = " << parameters.cuda_blockSize;
    #else
                    LOG(WARNING) << "\tThis build does not support CUDA.";
    #endif
                    break;
                default:
                    LOG(WARNING) << "\tUnknown parameter: " << opt;
                    printHelp();
                    break;
            }
        }
    }
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

    ifstream input(parameters.file);

    if (!input.good())
    {
        LOG(ERROR) << "Failed to open file: " << parameters.file;
        exit(0);
    }

    VLOG(0) << "Parsing config file " << parameters.file << "...";

    char line[512] = {0};

#define PARAMETER_SECTION 0
#define MOLECULE_SECTION 1
#define CROWDER_SECTION 2
#define SEGMENT_SECTION 3

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
        else if (strcmp(line, "segments") == 0)
        {
            section = SEGMENT_SECTION;
        }
        else if (section == MOLECULE_SECTION || section == CROWDER_SECTION)
        {
            moldata m;
            int result = 0;
            char pos_flag;

            result = sscanf(line, "%c(%f,%f,%f) r(%f,%f,%f,%f) %s %s", &pos_flag, &m.px, &m.py, &m.pz, &m.rx, &m.ry, &m.rz, &m.ra, m.pdbfilename, m.name);
            m.translate = (pos_flag == 't');

            if (result < 8)
            {
                LOG(WARNING) << "\tUnable to parse molecule: " << line;
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
                // TODO: if the molecule is named, add it to the name map
                // TODO: check if this comparison is right
                if (m.name && m.name[0] != '\0')
                {
                    parameters.mdata_map[string(m.name)] = parameters.mdata.size() - 1;
                }

                VLOG(1) << "\tAdded molecule from file " << m.pdbfilename << ".";
            }
        }
        else if (section == SEGMENT_SECTION)
        {
            char * molname = strtok(line, " ");            
            char * word = strtok(NULL, " ");
                        
            if (strcmp(word, "all") == 0)
            {
                parameters.mdata[parameters.mdata_map[string(molname)]].all_flexible = true;
            }
            else
            {
                segdata s;

                while (word != NULL)
                {
                    s.residue_indices.push_back(atoi(word) - 1);
                    word = strtok(NULL, " ");
                }
                
                parameters.mdata[parameters.mdata_map[string(molname)]].segments.push_back(s);
            }
        }
        else if (section == PARAMETER_SECTION)
        {
            char * key = strtok(line, " ");
            char * value = strtok(NULL, " ");

            VLOG(1) << "\tParameter " << key << " = " << value;

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
                sprintf(parameters.prefix, "%s_%d", value, parameters.pid);
            }
            else if (strcmp(key, "title") == 0)
            {
                int length = 0;
                length += sprintf(parameters.title + length, "%s", value);
                value = strtok(NULL," ");
                while (value != NULL)
                {
                    length += sprintf(parameters.title + length, " %s", value);
                    value = strtok(NULL," ");
                }
            }
            else {
                LOG(WARNING) << "\tUnknown parameter: " << key;
            }
        }
    }

    input.close();
}

void Simulation::check_and_modify_parameters()
{
    VLOG(0) << "Checking parameters for sanity...";

    if (parameters.bound <= 0)
    {
        LOG(WARNING) << "\tWARNINGING: Bounding value too small; setting equal to " << BOUNDING_VALUE << ".";
        parameters.bound = BOUNDING_VALUE;
    }

    if (parameters.temperatureMax < parameters.temperatureMin)
    {
        LOG(WARNING) << "\tWARNINGING: Maximum temperature < minimum temperature; swapping " << parameters.temperatureMax << " and " << parameters.temperatureMax << ".";
        float tmp = parameters.temperatureMax;
        parameters.temperatureMax = parameters.temperatureMin;
        parameters.temperatureMin = tmp;
    }

    if (parameters.threads > parameters.replicas)
    {
        parameters.threads = parameters.replicas;
        LOG(WARNING) << "\tWARNINGING: Threads > replicas; setting threads equal to " << parameters.threads << ".";
    }

    parameters.max_replicas_per_thread = int(ceil(float(parameters.replicas) / float(parameters.threads)));
    int spaces = parameters.replicas - parameters.max_replicas_per_thread * parameters.threads;
    int unused_threads = spaces / parameters.max_replicas_per_thread; // integer division
    if (unused_threads)
    {
        LOG(ERROR) << "\tERROR: After assignment of " << parameters.replicas << " replicas to " << parameters.threads << " threads with " << parameters.max_replicas_per_thread << " replicas per thread, " << unused_threads << " threads are left unused. You must either increase the number of replicas or decrease the number of threads. Exiting.";
        exit(0);
    }

#if USING_CUDA
    int availableGpus;
    cudaGetDeviceCount(&availableGpus);
    if (parameters.gpus > availableGpus)
    {
        parameters.gpus = availableGpus;
        LOG(WARNING) << "\tWARNINGING: Too many GPUs; setting equal to " << parameters.gpus << ".";
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
            LOG(WARNING) << "\tWARNINGING: Too many streams; setting equal to " << parameters.streams << ".";
        }
    }

    int length = strlen(parameters.prefix);
#if REPULSIVE_CROWDING
    length += sprintf(parameters.prefix + length, "_repcrowding");
#endif

#if LJ_REPULSIVE
    length += sprintf(parameters.prefix + length, "_repLJ");
#endif

#if LJ_OFF
    length += sprintf(parameters.prefix + length, "_LJoff");
#endif
}

void Simulation::init(int argc, char **argv, int pid)
{
    // Get the parameters

    parameters.pid = pid;

    getArgs(argc, argv, true);
    loadArgsFromFile();
    getArgs(argc, argv, false);

    // sanity check and calculation of some secondary parameters
    check_and_modify_parameters();

    printArgs();

    VLOG(0) << "Loading amino acid data " << AMINOACIDDATASOURCE << " and pair lookup table " << LJPDSOURCE;
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);

#ifdef displayLJpotentials
    printPotentials(&aminoAcidData);
#endif

    VLOG(0) << "Initialising first replica...";

    // Create the initial replica, and compare CPU and GPU potential

    // TODO: remove magic number; make initial array size a constant
    initialReplica.init_first_replica(parameters, aminoAcidData, 30);
    

//     VLOG(1) << "Performing potential test...\n");
#if INCLUDE_TIMERS
    initialReplica.initTimers();
#endif
    
    initialReplica.potential = initialReplica.E().total() + initialReplica.internal_molecule_E(true).total();

    if (initialReplica.nonCrowderCount < initialReplica.moleculeCount)
    {
#if REPULSIVE_CROWDING
        VLOG(0) << "\tCrowding is modelled using: u(r) = (6/r)^(12).";
#else
        VLOG(0) << "\tCrowding is modelled using the full potential calculations.";
#endif
    }

    // File stuff
    VLOG(0) << "Output files will be written to output " << parameters.prefix;

    // We need to create these directories in order to open files for writing
    // We need the files now so we can store them in the simulation data
    // TODO: maybe move all this stuff to run
    char mkdir_command[256];
    memset(mkdir_command, 0, 256);
    sprintf(mkdir_command, "mkdir -p output/%s/pdb",  parameters.prefix);
    int make_dirs = system(mkdir_command);

    initSamplingFiles();

    // now set up all the replicas

    VLOG(0) << "Initialising child replicas...";

    geometricTemperature = pow(double(parameters.temperatureMax/parameters.temperatureMin),double(1.0/double(parameters.replicas-1)));
    geometricTranslate = pow(double(MAX_TRANSLATION/MIN_TRANSLATION),double(1.0/double(parameters.replicas-1)));
    geometricRotation = pow(double(MAX_ROTATION/MIN_ROTATION),double(1.0/double(parameters.replicas-1)));

    for (size_t i = 0; i < parameters.replicas; i++)
    {
        replica[i].init_child_replica(initialReplica, int(i), geometricTemperature, geometricRotation, geometricTranslate, parameters);
        VLOG(1) << "\tReplica " << int(i) << "; temperature: " << replica[i].temperature << "; translate step: " << replica[i].translateStep << "; rotate step: " << replica[i].rotateStep;
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
    
    // Find the replica which is closest to 300K
    
    near_300k = temperature[0];
    
    for (map<float, int>::iterator p = position_of_temperature.begin(); p != position_of_temperature.end(); p++) {
        if (fabs(p->first - 300.0f) < fabs(near_300k - 300.0f)) {
            near_300k = p->first;
        }
    }
    
    // TODO: print out the correct temperature instead of 300K later (save the value on the class).
    _300kReplica = &replica[position_of_temperature[near_300k]];

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
        memset(data[i].prefix, 0, 256);
        strcpy(data[i].prefix, parameters.prefix);
        memset(data[i].title, 0, 256);
        strcpy(data[i].title, parameters.title);
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
        VLOG(1) << "\tAssigning thread " << i << " to GPU " << data[i].GPUID;

        data[i].max_replicas_per_thread = parameters.max_replicas_per_thread;

        // we have already checked that the replicas can be distributed like this
        data[i].replicas_in_this_thread = parameters.max_replicas_per_thread;
        if (i == parameters.threads - 2) // last thread
        {
            data[i].replicas_in_this_thread -= parameters.replicas - parameters.max_replicas_per_thread * parameters.threads;
        }

#if CUDA_STREAMS
        // the stream/replica ration must be a whole number otherwise there will be lots of waste, ie dormant streams etc
        data[i].replicas_per_stream = int(ceil(float(parameters.max_replicas_per_thread)/float(parameters.streams / parameters.threads)));
        data[i].streams_per_thread  = data[i].replicas_in_this_thread/data[i].replicas_per_stream;
#endif
    }

    cout.flush();
}

void Simulation::run()
{
    // we can't use pthreads and CUDA at the moment, but we are going to use streams

    VLOG(0) << "Beginning simulation...";

    // TODO: add stuff for resuming here

    if (initialReplica.moleculeCount == 0)  // make sure something is loaded
    {
        LOG(ERROR) << "No molecules loaded. Aborting run.";
        return;
    }

    writeFileIndex();
    writeSamplingFileHeaders();

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(RELoopTimer) );
#endif

    VLOG(0) << "--- Launching threads ---";

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

        if (parameters.viewConditions)
        {
            GlutDisplay();
        }
#endif

        VLOG(0) << "Replica Exchange step " << steps << " of " << parameters.REsteps << " complete (Fraction bound @ " << near_300k << "K: " << _300kReplica->accumulativeFractionBound;

        exchange_frequency();
    }

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStartTimer(MCLoopTimer) );
#endif

    pthread_cond_broadcast(&waitingThreadCond);
    pthread_mutex_unlock(&waitingCounterMutex);  // release the mutex so MC threads can continue.

    VLOG(0) << "--- Replica Exchanges Complete.---";
    VLOG(0) << "--- Waiting for threads to exit. ---";

    // join the threads that have finished
    for (int i=0; i<parameters.threads; i++)
        pthread_join(thread[i],NULL);

#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(MCLoopTimer) );
#endif

    VLOG(0) << "--- All threads complete.---";
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

    VLOG(0) << "--- Simulation finished.---";


#if INCLUDE_TIMERS
    CUT_SAFE_CALL( cutStopTimer(RELoopTimer) );
    VLOG(0) << "Simulation Timers";
    VLOG(0) << "MC Loop: Tot  " << cutGetTimerValue(MCLoopTimer) << " ms  Ave " << cutGetTimerValue(MCLoopTimer)/float(parameters.REsteps) << "ms (" << parameters.MCsteps << " steps, " << parameters.replicas << " replicas, " << parameters.threads << " threads, " << parameters.streams << " streams)";
    VLOG(0) << "Simulation: " << cutGetTimerValue(RELoopTimer) << " ms  (" << parameters.REsteps << " exchanges)";
#endif

    pthread_mutex_lock(&writeFileMutex);
    closeSamplingFiles();
    pthread_mutex_unlock(&writeFileMutex);

    VLOG(0) << "Simulation done.";

#if INCLUDE_TIMERS
    for (size_t i=0; i<parameters.replicas; i++)
    {
        VLOG(0) << "Replica " << i << " timers";
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

#if USING_CUDA
#if CUDA_STREAMS
void setup_CUDA(int device_id, float box_dimension, float * device_LJ_potentials, AminoAcids * amino_acid_data, cudaStream_t * streams, int streams_per_thread)
#else
void setup_CUDA(int device_id, float box_dimension, float * device_LJ_potentials, AminoAcids * amino_acid_data)
#endif // CUDA_STREAMS
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

#if CUDA_STREAMS
    // TODO: check whether this needs to happen before any of the things above.
    for (int i = 0; i < streams_per_thread; i++)
    {
        cudaStreamCreate(&streams[i]);
    }
#endif // CUDA_STREAMS
}

#if CUDA_STREAMS
void teardown_CUDA(float * device_LJ_potentials, cudaStream_t * streams, int streams_per_thread)
#else
void teardown_CUDA(float * device_LJ_potentials)
#endif // CUDA_STREAMS
{
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
    unbindLJTexture();
    cutilCheckMsg("Error freeing texture");
#endif

    cudaFree(device_LJ_potentials);

#if CUDA_STREAMS
    // TODO: check whether this needs to happen before any of the things above.
    for (int i = 0; i < streams_per_thread; i++)
    {
        cudaStreamDestroy(streams[i]);
    }
#endif // CUDA_STREAMS
}

#endif // USING_CUDA

void *MCthreadableFunction(void *arg)
{
    SimulationData *data = (SimulationData *) arg;

    int replica_offset = data->index * data->max_replicas_per_thread;

    VLOG(0) << "--- Monte-Carlo thread " << int(data->index + 1) << " running. ---";

// TODO: why do we do this?
#if GLVIS
    GLreplica = &data->replica[0];
#endif

#ifdef VERBOSE_THREAD_LOGGING
    pthread_mutex_lock( data->logMutex );
    for (int tx = 0; tx < data->replicas_in_this_thread; tx++)
    {
        VLOG(0) << "+ Thread " << data->index << " running replica " << tx + data->index * data->max_replicas_per_thread;
    }
    pthread_mutex_unlock( data->logMutex );
#endif

#if USING_CUDA
    // create streams for each subsequent MCSearch
    // ensure the replica can find the lookup table
    // initialise data on the device and copy the initial batch

    float * deviceLJpTmp;
#if CUDA_STREAMS
    cudaStream_t streams[16];   // reserve 16 stream slots but create only as many as needed
    setup_CUDA(data->GPUID, data->bound, deviceLJpTmp, &data->replica[replica_offset].aminoAcids, streams, data->streams_per_thread);
#else
    setup_CUDA(data->GPUID, data->bound, deviceLJpTmp, &data->replica[replica_offset].aminoAcids);
#endif // CUDA_STREAMS

    for (int tx = 0; tx < data->replicas_in_this_thread; tx++)
    {
#if CUDA_STREAMS
        data->replica[tx + replica_offset].setup_CUDA(deviceLJpTmp, streams, tx % data->streams_per_thread);
#else
        data->replica[tx + replica_offset].setup_CUDA(deviceLJpTmp);
#endif
    }
#endif // USING_CUDA

    // run the MC loop
#if !CUDA_MC && !CUDA_STREAMS   // if no MC on device
    int mcstep = 0;
    while (mcstep < data->MCsteps)
    {
        // do all MC steps before RE
        for (int tx = 0; tx < data->replicas_in_this_thread; tx++)
        {
            // to sample at the correct rate split this into another set of loops
            for (int s = 0; s < data->MC_steps_per_RE / data->sampleFrequency; s++)
            {
                data->replica[tx + replica_offset].MCSearch(data->sampleFrequency, mcstep);

                if (mcstep + s * data->sampleFrequency >= data->sampleStartsAfter)
                {
                    // if E < -1.1844 kcal/mol then its bound
                    data->replica[tx + replica_offset].sample(data, mcstep + s*data->sampleFrequency, BOUND_ENERGY_VALUE);
                }
            }
        }
        mcstep += data->MC_steps_per_RE;

        // do replica exchange
        pthread_mutex_lock(data->waitingCounterMutex);                 // lock the counter

        data->waitingThreadCount[0] = data->waitingThreadCount[0] + 1;
        if (data->waitingThreadCount[0] == data->threads)                 // if all threads in waiting state
        {
            pthread_cond_signal(data->waitingReplicaExchangeCond);
        }
        if (mcstep < data->MCsteps)                       // wait if another MC loop must be done
        {
            pthread_cond_wait(data->waitingThreadCond, data->waitingCounterMutex); // wait for replica exchange to finish.
            // NB! unlock the mutex locked upon the condition above
            // being met because this will unblock all threads such
            // that they will continue concurrently, if not unlocked
            // other threads will run sequentially

            //  pthread_mutex_unlock(data->waitingCounterMutex);
        }
        pthread_mutex_unlock(data->waitingCounterMutex);


    } // continue MC
#endif


#if CUDA_STREAMS    // use streams
    int mcstep = 0;
    while (mcstep < data->MCsteps)
    {
        // to sample at the correct rate split this into another set of loops
        for (int mcx = 0; mcx < data->MC_steps_per_RE; mcx++) // at each mc step
        {
            for (int index = 0; index < data->replicas_in_this_thread; index += data->replicas_per_stream)
            {
                for (int rps = 0; rps < data->replicas_per_stream; rps++)
                {
                    // batch replicas such that no stream is shared per batch
                    data->replica[replica_offset + index + rps].MCSearchMutate(mcstep);
                    data->replica[replica_offset + index + rps].MCSearchEvaluate(mcstep);
                }

                for (int rps = 0; rps < data->replicas_per_stream; rps++)
                {
                    data->replica[replica_offset + index + rps].MCSearchAcceptReject(mcstep);
                }

                for (int rps = 0; rps < data->replicas_per_stream; rps++)
                {
                    //sampleAsync
                    if (mcstep % data->sampleFrequency == data->sampleFrequency - 1 && mcstep >= data->sampleStartsAfter)
                    {
                        // We check if step % frequency == frequency - 1 so that sampling is done after the *last* step of each block of steps, not the first.
                        // This matches the behaviour of the synchronous CUDA version.
                        data->replica[replica_offset + index + rps].sample(data, mcstep - data->sampleFrequency + 1, BOUND_ENERGY_VALUE);
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
    for (int tx = 0; tx < data->replicas_in_this_thread; tx++)
    {
        data->replica[tx + replica_offset].teardown_CUDA();
    }

#if CUDA_STREAMS
    teardown_CUDA(deviceLJpTmp, streams, data->streams_per_thread);
#else
    teardown_CUDA(deviceLJpTmp);
#endif // CUDA_STREAMS

#endif // USING_CUDA
    VLOG(0) << "--- Monte-Carlo thread " << int(data->index + 1) << " exited. ---";
    return 0;
}
// END OF *MCthreadableFunction

void Simulation::initSamplingFile(const char * name, FILE ** file_addr)
{
    char filename[256];
    FILE * file;

    sprintf(filename, "output/%s/%s", parameters.prefix, name);
    file = fopen (filename,"a+"); // attempt append a file of the same name
    if (!file)
    {
        file = fopen(filename, "wt"); // create if that failed, and open and append failed (eg: file permissions r--)
    }
    if (!file)
    {
        LOG(ERROR) << "Cannot open/create file: " << filename;
        return;
    }

    *file_addr = file;
}

void Simulation::closeSamplingFile(const char * name, FILE ** file_addr)
{
    fclose(*file_addr);
    VLOG(0) << "Closing sampling file: output/" << parameters.prefix << "/" <<  name;
}

void Simulation::initSamplingFiles()
{
    initSamplingFile("fractionBound", &fractionBoundFile);
    initSamplingFile("boundconformations", &boundConformationsFile);
    initSamplingFile("acceptance_ratios", &acceptanceRatioFile);
    initSamplingFile("exchange_freq", &exchangeFrequencyFile);
}

void Simulation::writeSamplingFileHeaders()
{
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

void Simulation::writeFileIndex()
{
    char fileindex[256];
    sprintf(fileindex,"output/%s/fileindex", parameters.prefix);
    FILE * fileindexf = fopen (fileindex,"w");
    fprintf(fileindexf,"index molecule_file_path crowder(Y/N)\n");

    for (int i = 0; i < parameters.mdata.size(); i++) {
        fprintf(fileindexf, "%2d %s %s", i, parameters.mdata[i].pdbfilename, (parameters.mdata[i].crowder ? "Y\n" : "N\n"));
    }

    fclose(fileindexf);
}

void Simulation::printArgs()
{
    VLOG(1) << "Final parameters:";
    VLOG(1) << "\tThreads: " << parameters.threads;
    VLOG(1) << "\tStreams: " << parameters.streams;
    VLOG(1) << "\tGPUs: " << parameters.gpus;
    VLOG(1) << "\tReplicas: " << parameters.replicas;

    VLOG(1) << "\tMC steps: " << parameters.MCsteps;
    VLOG(1) << "\tRE steps: " << parameters.REsteps;
    VLOG(1) << "\tSampling frequency (MC steps): " << parameters.sampleFrequency;
    VLOG(1) << "\tSampling starts after (MC steps): " << parameters.sampleStartsAfter;

    VLOG(1) << "\tBounding box size: " << parameters.bound;
    VLOG(1) << "\tNon-crowder molecules: " << parameters.nonCrowders;
    VLOG(1) << "\tMaximum temperature: " << parameters.temperatureMax;
    VLOG(1) << "\tMinimum temperature: " << parameters.temperatureMin;
}
