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
    RELoopTimer = NULL;
    MCLoopTimer = NULL;
    sdkCreateTimer(&RELoopTimer);
    sdkCreateTimer(&MCLoopTimer);
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
    sdkDeleteTimer(&RELoopTimer);
    sdkDeleteTimer(&MCLoopTimer);
#endif
}

void Simulation::printHelp()
{
    printSettings();
    cout << "Usage: cgppd -f <filename> [-h] [-p] [-q] [-v x] [-t x] [-s x] [-g x] [-z x] [-m x ] [-a x] [-e x] [-r x] [-o x] [-b x] [-n x] [-x x] [-d x]"<< endl;
    cout << "\t-h|--help: show this dialog" << endl;
    cout << "\t-f|--file <file>:    Input config file" << endl;
    cout << "\t-p|--preview:        Use the open GL preview of this configuration, performs no simulation" << endl;
    cout << "\t-q|--nosim:          Do everything except the simulation (for use with -p)" << endl;
    cout << "\t-v|--v:              Verbosity level" << endl;
    cout << "\nThe following values override those in the config file" << endl;
    cout << "\t-t|--threads x:      The number of CPU/pthreads to use" << endl;
    cout << "\t-s|--streams x:      The number of CUDA streams" << endl;
    cout << "\t-g|--gpus x:         The number of GPUs to use" << endl;
    cout << "\t-z|--gpuoffset x:    Use GPUs starting with this one" << endl;
    cout << "\t-m|--mcsteps x:      How many MC steps to perform per replica" << endl;
    cout << "\t-a|--sampleafter x:  How many MC steps to perform before beginning sampling" << endl;
    cout << "\t-e|--resteps x:      How many replica exchanges to perform" << endl;
    cout << "\t-r|--replicas x:     How many replicas" << endl;
    cout << "\t-o|--output x:       The output prefix for files created by the simulation" << endl;
    cout << "\t-b|--boundary x:     The bounding box edge length" << endl;
    cout << "\t-x|--tmax x:         The temperature of the highest replica" << endl;
    cout << "\t-n|--tmin x:         The temperature of the lowest replica" << endl;
    cout << "\t-d|--blockdim x:     Number of threads per CUDA block" << endl;
}

void Simulation::getArgs(int argc, char **argv, bool first_pass)
{
    const struct option long_options[] =
    {
        {"file", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"verbosity", required_argument, 0, 'v'},
        
        {"preview", no_argument, 0, 'p'},
        {"nosim", no_argument, 0, 'q'},

        {"threads", required_argument, 0, 't'},
        {"streams", required_argument, 0, 's'},
        {"gpus", required_argument, 0, 'g'},
        {"gpuoffset", required_argument, 0, 'z'},
        
        {"mcsteps", required_argument, 0, 'm'},
        {"sampleafter", required_argument, 0, 'a'},
        {"resteps", required_argument, 0, 'e'},
        {"replicas", required_argument, 0, 'r'},
        
        {"output", required_argument, 0, 'o'},

        {"boundary", required_argument, 0, 'b'},
        {"tmax", required_argument, 0, 'x'},
        {"tmin", required_argument, 0, 'n'},
        {"blockdim", required_argument, 0, 'd'},

        {0, 0, 0, 0},
    };

    const char short_options[] = "hf:pqt:s:g:z:m:a:e:r:o:b:x:n:d:v:w:";

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
                    throw "help";
                    break;
                case 'f':
                    strcpy(parameters.file, optarg);
                    parameters.inputFile = true;
                    break;
                case 'v':
                    switch(atoi(optarg))
                    {
                        case 0:
                            LOGOG_SET_LEVEL(LOGOG_LEVEL_NONE);
                            break;
                        case 1:
                            // keep LOGOG_LEVEL_ERROR
                            break;
                        case 2:
                            LOGOG_SET_LEVEL(LOGOG_LEVEL_WARN);
                            break;
                        case 3:
                            LOGOG_SET_LEVEL(LOGOG_LEVEL_INFO);
                            break;
                        case 4:
                            LOGOG_SET_LEVEL(LOGOG_LEVEL_ALL);
                            break;
                    }
                    
                    break;
                default:
                    // Ignore all other options in this pass
                    break;
            }
            
        }

        if (!parameters.inputFile)
        {
            printHelp();
            throw "No configuration file provided.";
        }
        
        LOGOG_DEBUG("Will read configuration from file: %s", parameters.file);
    } else {
        LOGOG_INFO("Parsing commandline parameters...");

        while (1)
        {
            int opt = getopt_long(argc, argv, short_options, long_options, &opt_index);

            if (opt == -1)
                break;

            switch(opt)
            {
                case 'h':
                    // ignore
                case 'f':
                    // ignore
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
                    LOGOG_DEBUG("\tWill show OpenGL preview.");
    #else
                    LOGOG_WARN("\tThis build does not support OpenGL.");
    #endif
                    break;
                case 'q':
                    parameters.skipsimulation = true;
                    LOGOG_DEBUG("\tWill skip simulation.");
                    break;
                case 't':
                    parameters.threads = atoi(optarg);
                    LOGOG_DEBUG("\tParameter threads = %d", parameters.threads);
                    break;
                case 's':
                    parameters.streams = atoi(optarg);
                    LOGOG_DEBUG("\tParameter streams = %d", parameters.streams);
                    break;
                case 'g':
                    parameters.gpus = atoi(optarg);
                    LOGOG_DEBUG("\tParameter gpus = %d", parameters.gpus);
                    break;
                case 'z':
                    parameters.gpuoffset = atoi(optarg);
                    LOGOG_DEBUG("\tParameter gpuoffset = %d", parameters.gpuoffset);
                    break;
                case 'm':
                    parameters.MCsteps = atoi(optarg);
                    LOGOG_DEBUG("\tParameter MCsteps = %d", parameters.MCsteps);
                    break;
                case 'a':
                    parameters.sampleStartsAfter = atoi(optarg);
                    LOGOG_DEBUG("\tParameter sampleStartsAfter = %d", parameters.sampleStartsAfter);
                    break;
                case 'e':
                    parameters.REsteps = atoi(optarg);
                    LOGOG_DEBUG("\tParameter REsteps = %d", parameters.REsteps);
                    break;
                case 'r':
                    parameters.replicas = atoi(optarg);
                    LOGOG_DEBUG("\tParameter replicas = %d", parameters.replicas);
                    break;
                case 'o':
                    sprintf(parameters.prefix, "%s_%d", optarg, parameters.pid);
                    LOGOG_DEBUG("\tParameter prefix = %s", parameters.prefix);
                    break;
                case 'b':
                    parameters.bound = atof(optarg);
                    LOGOG_DEBUG("\tParameter bound = %.1f", parameters.bound);
                    break;
                case 'x':
                    parameters.temperatureMax = atof(optarg);
                    LOGOG_DEBUG("\tParameter temperatureMax = %.1f", parameters.temperatureMax);
                    break;
                case 'n':
                    parameters.temperatureMin = atof(optarg);
                    LOGOG_DEBUG("\tParameter temperatureMin = %.1f", parameters.temperatureMin);
                    break;
                case 'd':
    #if USING_CUDA
                    parameters.cuda_blockSize = atoi(optarg);
                    parameters.auto_blockdim = false;
                    LOGOG_DEBUG("\tParameter cuda_blockSize = %d", parameters.cuda_blockSize);
    #else
                    LOGOG_WARN("\tThis build does not support CUDA.");
    #endif
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
        throw "Failed to open config file.";
    }
    
    LOGOG_INFO("Parsing config file %s...", parameters.file);

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
                LOGOG_WARN("\tUnable to parse molecule: %s", line);
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
                
                LOGOG_DEBUG("\tAdded molecule from file %s.", m.pdbfilename);
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
            
            LOGOG_DEBUG("\tParameter %s = %s", key, value);

            if (strcmp(key, "gpus") == 0)
            {
                parameters.gpus = atoi(value);
            }
            else if (strcmp(key, "gpuoffset") == 0)
            {
                parameters.gpuoffset = atoi(value);
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
                LOGOG_WARN("\tUnknown parameter: %s", key);
            }
        }
    }

    input.close();
}

void Simulation::check_and_modify_parameters()
{
    LOGOG_INFO("Checking parameters for sanity...");

    if (parameters.bound <= 0)
    {
        LOGOG_WARN("\tBounding value too small; setting equal to %.1f.", BOUNDING_VALUE);
        parameters.bound = BOUNDING_VALUE;
    }

    if (parameters.temperatureMax < parameters.temperatureMin)
    {
        LOGOG_WARN("\tMaximum temperature < minimum temperature; swapping %.1f and %.1f.", parameters.temperatureMin, parameters.temperatureMax);
        float tmp = parameters.temperatureMax;
        parameters.temperatureMax = parameters.temperatureMin;
        parameters.temperatureMin = tmp;
    }

    if (parameters.threads > parameters.replicas)
    {
        parameters.threads = parameters.replicas;
        LOGOG_WARN("\tThreads > replicas; setting threads equal to %d.", parameters.threads);
    }

    parameters.max_replicas_per_thread = int(ceil(float(parameters.replicas) / float(parameters.threads)));
    int spaces = parameters.max_replicas_per_thread * parameters.threads - parameters.replicas;
    int unused_threads = spaces / parameters.max_replicas_per_thread; // integer division
    if (unused_threads)
    {
        LOGOG_ERROR("\tAfter assignment of %d replicas to %d threads with %d replicas per thread, %d threads are left unused. You must either increase the number of replicas or decrease the number of threads.", parameters.replicas, parameters.threads, parameters.max_replicas_per_thread, unused_threads);
        throw "Wasted threads.";
    }

#if USING_CUDA
    int availableGpus;
    cudaGetDeviceCount(&availableGpus);
    if (parameters.gpuoffset >= availableGpus) {
        LOGOG_ERROR("\tCannot allocate %d GPUs with an offset of %d because only %d are available.", parameters.gpus, parameters.gpuoffset, availableGpus);
        throw "GPU allocation impossible.";
    } else if (parameters.gpus + parameters.gpuoffset > availableGpus) {
        LOGOG_WARN("\tCannot allocate %d GPUs with an offset of %d because only %d are available. Reducing number of GPUs used to %d.", parameters.gpus, parameters.gpuoffset, availableGpus, (availableGpus - parameters.gpuoffset));
        parameters.gpus = availableGpus - parameters.gpuoffset;
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
            LOGOG_WARN("\tToo many streams; setting equal to %d.", parameters.streams);
        }
    }
}

void Simulation::init(int argc, char **argv, int pid)
{
    // Get the parameters

    parameters.pid = pid;

    getArgs(argc, argv, true);
    loadArgsFromFile();
    getArgs(argc, argv, false);

    // File stuff
    
    // TODO: move this back where it was (into run) so that no directories are made until the simulation is actually running

    // We need to create these directories in order to open files for writing
    // We need the files now so we can store them in the simulation data
    // Move this earlier, so we can put the logs here
    char mkdir_command[256];
    memset(mkdir_command, 0, 256);
    sprintf(mkdir_command, "mkdir -p output/%s/pdb",  parameters.prefix);
    int make_dirs = system(mkdir_command);

    LOGOG_INFO("Created output directory: output/%s", parameters.prefix);

    initSamplingFiles();
    
    ostringstream settings;
    printSettings(settings);
    LOGOG_INFO(settings.str().c_str());

    // sanity check and calculation of some secondary parameters
    check_and_modify_parameters();

    printArgs();
    
    LOGOG_DEBUG("Loading amino acid data %s and pair lookup table %s", AMINOACIDDATASOURCE, LJPDSOURCE);
    
    aminoAcidData.init(AMINOACIDDATASOURCE, LJPDSOURCE);

#ifdef displayLJpotentials
    aminoAcidData.printPotentials();
#endif

    LOGOG_DEBUG("Initialising first replica...");

    // Create the initial replica, and compare CPU and GPU potential

    // TODO: remove magic number; make initial array size a constant
    initialReplica.init_first_replica(parameters, aminoAcidData, 30);
    
    if (initialReplica.moleculeCount == 0)  // make sure something is loaded
    {
        throw "No molecules loaded.";
    }
    
#if INCLUDE_TIMERS
    initialReplica.initTimers();
#endif
    
    initialReplica.potential = initialReplica.E().total() + initialReplica.internal_molecule_E(true).total();

    if (initialReplica.nonCrowderCount < initialReplica.moleculeCount)
    {
#if REPULSIVE_CROWDING
        LOGOG_INFO("\tCrowding is modelled using: u(r) = (6/r)^(12).");
#else
        LOGOG_INFO("\tCrowding is modelled using the full potential calculations.");
#endif
    }

    // now set up all the replicas

    LOGOG_DEBUG("Initialising child replicas...");

    geometricTemperature = pow(double(parameters.temperatureMax/parameters.temperatureMin),double(1.0/double(parameters.replicas-1)));
    geometricTranslate = pow(double(MAX_TRANSLATION/MIN_TRANSLATION),double(1.0/double(parameters.replicas-1)));
    geometricRotation = pow(double(MAX_ROTATION/MIN_ROTATION),double(1.0/double(parameters.replicas-1)));

    for (size_t i = 0; i < parameters.replicas; i++)
    {
        replica[i].init_child_replica(initialReplica, int(i), geometricTemperature, geometricRotation, geometricTranslate, parameters);
        LOGOG_DEBUG("\tReplica %d; temperature: %.1f; translate step: %.1f; rotate step: %.1f", i, replica[i].temperature, replica[i].translateStep, replica[i].rotateStep);
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
        LOGOG_DEBUG("Assigning parameters for thread %d", i);
        data[i].replica = replica; // common pointer to all the replicas!
        data[i].replicaCount = parameters.replicas;
        LOGOG_DEBUG("\tReplicas: %d", data[i].replicaCount);
        data[i].index = i;
        data[i].threads = parameters.threads;
        LOGOG_DEBUG("\tthreads: %d", data[i].threads);
        data[i].streams = parameters.streams;
        LOGOG_DEBUG("\tstreams: %d", data[i].streams);
        data[i].MCsteps = parameters.MCsteps;
        LOGOG_DEBUG("\tMCsteps: %d", data[i].MCsteps);
        data[i].REsteps = parameters.REsteps;
        LOGOG_DEBUG("\tREsteps: %d", data[i].REsteps);
        data[i].MC_steps_per_RE = parameters.MCsteps/parameters.REsteps;
        LOGOG_DEBUG("\tMC_steps_per_RE: %d", data[i].MC_steps_per_RE);
        data[i].sampleFrequency = parameters.sampleFrequency;
        LOGOG_DEBUG("\tsampleFrequency: %d", data[i].sampleFrequency);
        data[i].sampleStartsAfter = parameters.sampleStartsAfter;
        LOGOG_DEBUG("\tsampleStartsAfter: %d", data[i].sampleStartsAfter);
        data[i].bound = parameters.bound;
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
        data[i].endBarrier = &endBarrier;

        // assign gpus in rotation per thread, t0 = gpu0, t1 = gpu1 etc
        // NEW: added optional offset
        // % #gpus so they share if threads > gpus
        // will perform best if threads:gpus = 1:1
        data[i].GPUID = i % parameters.gpus + parameters.gpuoffset;
        LOGOG_DEBUG("\tGPUID: %d", data[i].GPUID);
        data[i].max_replicas_per_thread = parameters.max_replicas_per_thread;
        LOGOG_DEBUG("\tmax_replicas_per_thread: %d", data[i].max_replicas_per_thread);

        // we have already checked that the replicas can be distributed like this
        data[i].replicas_in_this_thread = parameters.max_replicas_per_thread;
        if (i == parameters.threads - 2) // last thread
        {
            data[i].replicas_in_this_thread -= parameters.replicas - parameters.max_replicas_per_thread * parameters.threads;
        }
        LOGOG_DEBUG("\treplicas_in_this_thread: %d", data[i].replicas_in_this_thread);

#if CUDA_STREAMS
        // the stream/replica ration must be a whole number otherwise there will be lots of waste, ie dormant streams etc
        data[i].replicas_per_stream = int(ceil(float(parameters.max_replicas_per_thread)/float(parameters.streams / parameters.threads)));
        LOGOG_DEBUG("\treplicas_per_stream: %d", data[i].replicas_per_stream);
        data[i].streams_per_thread  = data[i].replicas_in_this_thread/data[i].replicas_per_stream;
        LOGOG_DEBUG("\tstreams_per_thread: %d", data[i].streams_per_thread);
#endif
    }

    cout.flush();
}

void Simulation::resume_from(char * directory_name)
{
    // do stuff
}


void Simulation::run()
{
    LOGOG_INFO("Beginning simulation...");

    // TODO: add stuff for resuming here

    writeFileIndex();
    writeSamplingFileHeaders();

#if INCLUDE_TIMERS
    sdkStartTimer(&RELoopTimer);
#endif

    LOGOG_INFO("--- Launching threads ---");

    pthread_mutex_lock(&waitingCounterMutex);
    
    pthread_barrier_init(&endBarrier, NULL, parameters.threads);

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
        sdkStartTimer(&MCLoopTimer);
#endif

        pthread_cond_wait(&waitingReplicaExchangeCond, &waitingCounterMutex);

#if INCLUDE_TIMERS
        sdkStopTimer(&MCLoopTimer);
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
        
        LOGOG_INFO("Replica Exchange step %d of %d complete (Fraction bound @ %.1fK: %.1f", steps, parameters.REsteps, near_300k, _300kReplica->accumulativeFractionBound);

        exchange_frequency();
    }

#if INCLUDE_TIMERS
    sdkStartTimer(&MCLoopTimer);
#endif

    pthread_cond_broadcast(&waitingThreadCond);
    pthread_mutex_unlock(&waitingCounterMutex);  // release the mutex so MC threads can continue.
    
    LOGOG_INFO("--- Replica Exchanges Complete.---");
    LOGOG_INFO("--- Waiting for threads to exit. ---");

    // join the threads that have finished
    for (int i=0; i<parameters.threads; i++)
        pthread_join(thread[i],NULL);

#if INCLUDE_TIMERS
    sdkStopTimer(&MCLoopTimer);
#endif

    LOGOG_INFO("--- All threads complete.---");
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

#if INCLUDE_TIMERS
    sdkStopTimer(&RELoopTimer);
    LOGOG_INFO("Simulation Timers");
    LOGOG_INFO("MC Loop: Tot %.1fms  Ave %.1fms (%d MC steps, %d replicas, %d threads, %d streams", sdkGetTimerValue(&MCLoopTimer), sdkGetTimerValue(&MCLoopTimer)/float(parameters.REsteps), parameters.MCsteps, parameters.replicas, parameters.threads, parameters.streams);
    LOGOG_INFO("Simulation: %1.fms  (%d exchanges)", sdkGetTimerValue(&RELoopTimer), parameters.REsteps);
#endif

    pthread_mutex_lock(&writeFileMutex);
    closeSamplingFiles();
    pthread_mutex_unlock(&writeFileMutex);

#if INCLUDE_TIMERS
    for (size_t i=0; i<parameters.replicas; i++)
    {
        LOGOG_INFO("Replica %d timers", i);
        replica[i].printTimers();
    }
#endif

    LOGOG_INFO("--- Simulation finished.---");
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
    getLastCudaError("Failed to initialise CUDA runtime.");

    cudaSetDevice(device_id);
    getLastCudaError("Failed to pick device for the CUDA runtime.");

    // TODO: free this afterwards using CUDA_freeBoxDimension?
    CUDA_setBoxDimension(box_dimension);
    getLastCudaError("Failed to copy box dimensions to GPU.");

    // copy the LJpotentials to gpu memory in this thread context so it can access it
    cudaMalloc((void**)&device_LJ_potentials, LJArraySize);
    getLastCudaError("Failed to allocate contact potential memory on the GPU.");

    copyLJPotentialDataToDevice(device_LJ_potentials, amino_acid_data);
    getLastCudaError("Failed to copy contact potentials to device.");

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
    getLastCudaError("Error freeing texture");
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

// TODO: why do we do this?
#if GLVIS
    GLreplica = &data->replica[0];
#endif

#ifdef VERBOSE_THREAD_LOGGING
    pthread_mutex_lock( data->logMutex );
    for (int tx = 0; tx < data->replicas_in_this_thread; tx++)
    {
        LOGOG_INFO("+ Thread %d running replica %d", data->index, tx + data->index * data->max_replicas_per_thread);
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

    int res = pthread_barrier_wait(data->endBarrier);

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
        throw "Could not open / create sampling file.";
    }

    *file_addr = file;
}

void Simulation::closeSamplingFile(const char * name, FILE ** file_addr)
{
    fclose(*file_addr);
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
    memset(fileindex, 0, 256);
    sprintf(fileindex,"output/%s/fileindex", parameters.prefix);
    FILE * fileindexf = fopen (fileindex,"w");
    fprintf(fileindexf,"index molecule_file_path crowder(Y/N)\n");

    for (int i = 0; i < parameters.mdata.size(); i++) {
        fprintf(fileindexf, "%2d %s %s", i, parameters.mdata[i].pdbfilename, (parameters.mdata[i].crowder ? "Y\n" : "N\n"));
    }

    fclose(fileindexf);
}

void Simulation::printSettings(ostream & os)
{    
    os << "CGPPD version: " << HGVERSION << endl;
    os << "Compiled with:" << endl;

#if GLVIS
    os << "\tOpenGL support" << endl;
#endif

#if USING_CUDA
    os << "\tCUDA support" << endl;
#if CUDA_STREAMS
    os << "\t\tAsynchronous GPU calls (CUDA capability 1.1+ required)" << endl;
#endif // CUDA_STREAMS
    os << "\t\tTile size: " << TILE_DIM << endl;
    string mem_type;
#if LJ_LOOKUP_METHOD == SHARED_MEM
    mem_type = "Shared";
#elif LJ_LOOKUP_METHOD == CONST_MEM
    mem_type = "Constant";
#elif LJ_LOOKUP_METHOD == GLOBAL_MEM
    mem_type = "Global";
#elif LJ_LOOKUP_METHOD == TEXTURE_MEM
    mem_type = "Texture";
#endif // LJ_LOOKUP_METHOD
    os << "\t\tLJ lookup memory type: " << mem_type << endl;
#endif // USING_CUDA

#if COMPENSATE_KERNEL_SUM
    os << "\tKahan summation in kernels" << endl;
#endif
#if FLEXIBLE_LINKS
    os << "\tFlexible linkers" << endl;
#endif
#if LJ_REPULSIVE
    os << "\tLennard-Jones potentials always repulsive" << endl;
#endif
#if LJ_OFF
    os << "\tLennard-Jones potentials off" << endl;
#endif
#if ASSUME_POLYMER_FOLDING_TEST
    os << "\tAll residues assumed to be alanine (polymer test)" << endl;
#endif
}

void Simulation::printArgs()
{
    LOGOG_DEBUG("Final parameters:");
    LOGOG_DEBUG("\tThreads: %d", parameters.threads);
    LOGOG_DEBUG("\tStreams: %d", parameters.streams);
    LOGOG_DEBUG("\tGPUs: %d", parameters.gpus);
    LOGOG_DEBUG("\tGPU offset: %d", parameters.gpuoffset);
    LOGOG_DEBUG("\tReplicas: %d", parameters.replicas);
    LOGOG_DEBUG("\tMC steps: %d", parameters.MCsteps);
    LOGOG_DEBUG("\tRE steps: %d", parameters.REsteps);
    LOGOG_DEBUG("\tSampling frequency (MC steps): %d", parameters.sampleFrequency);
    LOGOG_DEBUG("\tSampling starts after (MC steps): %d", parameters.sampleStartsAfter);
    LOGOG_DEBUG("\tBounding box size: %.1f", parameters.bound);
    LOGOG_DEBUG("\tNon-crowder molecules: %d", parameters.nonCrowders);
    LOGOG_DEBUG("\tMaximum temperature: %.1f", parameters.temperatureMax);
    LOGOG_DEBUG("\tMinimum temperature: %.1f", parameters.temperatureMin);
}
