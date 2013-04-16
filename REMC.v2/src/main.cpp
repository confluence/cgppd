#include "main.h"

void printHelp(bool badArg)
{
    if (badArg)
        cout << " * Bad arguments" << endl;
    cout << "Usage: REMCDockingXXX -f <filename> [-c] [-h] [-v] [-0] [-t x] [-s x] [-g x] [-m x ] [-e x] [-r x] [-o x] [-b x] [-t0 x] [-t1 x] [-bx x]"<< endl;
    cout << "\t-h|--help: show this dialog" << endl;
    cout << "\t-c: Check, perform 10 reference potential sums, gpu vs cpu" << endl;
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
    cout << "\t-t0|--tmax x: 	The temperature of the highest replica" << endl;
    cout << "\t-t1|--tmin x:    The temperature of the lowest replica" << endl;
    cout << "\t-t1|--tmin x:    The temperature of the lowest replica" << endl;
    cout << "\t-bx|--blockdim x: Number of threads per CUDA block" << endl;

    exit(0);
}

// TODO: COPYPASTA! Surely, there is a standard library for commandline parameters?
bool getArgs(argdata * d, int argc, char **argv)
{
    d->resume = false;
    d->nonCrowders = 0;
    d->viewConditions = false;
    d->skipsimulation = false;

#if USING_CUDA
    d->cuda_blockSize = TILE_DIM;
    d->auto_blockdim = true;
#endif

    if (!d->inputFile)
    {
        d->gpus = 1;
        d->threads = THREAD_COUNT;
        d->streams = STREAM_COUNT;
        d->MCsteps = MC_STEPS;
        d->REsteps = REMC_STEPS;
        d->sampleFrequency = SAMPLE_FREQ;
        d->sampleStartsAfter = STEPS_BEFORE_SAMPLE;
        d->inputFile = false;
        d->replicas = REPLICA_COUNT;
        d->bound = BOUNDING_VALUE;
        d->temperatureMin = LOWTEMP;
        d->temperatureMax = HIGHTEMP;
    }

    memset(&d->logfile,0,256);
    sprintf(d->logfile,"output/%d_logfile",d->pid);

    if (argc <= 1)
    {
        cout << "No arguments." << endl;
        printHelp(false);
        return true;  // use default values
    }

    int i = 1;
    while (i<argc)
    {
        // threads to be used
        if (strcmp(argv[i],"-v")==0 || strcmp(argv[i],"--view")==0)
        {
#ifndef EnableOPENGL
            cout << "!!! This build does not support OpenGL and the -v option" << endl;
#else
            d->viewConditions = true;
#endif
            i++;
        }

        else if (strcmp(argv[i],"-q")==0 || strcmp(argv[i],"--nosim")==0)
        {
            d->skipsimulation = true;
            i++;
        }

#if USING_CUDA
        else if (strcmp(argv[i],"-bx")==0 || strcmp(argv[i],"--blockdim")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            d->cuda_blockSize = atoi(argv[i+1]);
            cout << "Block size changed to: " << cuda_blockSize << endl;
            i+=2;
            d->auto_blockdim = false;
        }
#endif
        // threads to be used
        else if (strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--threads")==0)
        {
            if (i+1>=argc)
                printHelp(true);

            d->threads = atoi(argv[i+1]);
            i+=2;
        }

        // streams to be used
        else if (strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--streams")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            d->streams = atoi(argv[i+1]);
            i+=2;
        }

        // gpus to be used
        else if (strcmp(argv[i],"-g")==0 || strcmp(argv[i],"--gpus")==0)
        {
            if (i+1>=argc)
                printHelp(true);

            d->gpus = atoi(argv[i+1]);
            i+=2;
        }

        //input file for molecules
        else if (strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--file")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            strcpy(d->file,argv[i+1]);
            d->inputFile = true;
            i+=2;
        }

        else if (strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--mcsteps")==0)
        {
            if (i+1>=argc)
                printHelp(true);

            d->MCsteps = atoi(argv[i+1]);
            i+=2;
        }

        else if (strcmp(argv[i],"-e")==0 || strcmp(argv[i],"--resteps")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            d->REsteps = atoi(argv[i+1]);
            i+=2;
        }

        else if (strcmp(argv[i],"-r")==0 || strcmp(argv[i],"--replicas")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            d->replicas = atoi(argv[i+1]);
            i+=2;
        }

        else if (strcmp(argv[i],"-o")==0 || strcmp(argv[i],"--output")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            strcpy(d->logfile,argv[i+1]);
            i+=2;
        }
        else if (strcmp(argv[i],"-b")==0 || strcmp(argv[i],"--boundary")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            d->bound = atof(argv[i+1]);
            i+=2;
        }
        else if (strcmp(argv[i],"-t0")==0 || strcmp(argv[i],"--tmin")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            d->temperatureMin = atof(argv[i+1]);
            i+=2;
        }
        else if (strcmp(argv[i],"-t1")==0 || strcmp(argv[i],"--tmax")==0)
        {
            if (i+1>=argc)
                printHelp(true);
            d->temperatureMax = atof(argv[i+1]);
            i+=2;
        }
        else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0)
        {
            printHelp(false);
        }
        else if (strcmp(argv[i],"-z")==0 || strcmp(argv[i],"--resume")==0)
        {
            cout << "!!! Checkpointing is not enabled in this build, cannot resume" << endl;
            i++;
        }
        else
        {
            i++;
        }
    }
    return false;
}

void loadArgsFromFile(argdata * parameters)
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
    if (parameters->inputFile)
    {
        ifstream input(parameters->file);
        if (!input.good())
        {
            cout << "Failed to open file: " << parameters->file << "\n";
            exit(0);
        }

        char line[512] = {0};
        char saveline[512];

        bool parameterSection = true; // parameters come first
        bool moleculeSection = false;
        bool crowdersSection = false;
        int lc = 0;
        // for each line in the file
        while (!input.eof())
        {
            input.getline(line,512);
            lc++;
            strcpy(saveline,line);
            char *token = strtok(line," ");

            if (token==NULL || line[0] == '#' || strlen(token)==0)
            {
                // comment or blank line
            }
            else if (strcmp(token,"files")==0)
            {
                parameterSection = false;
                moleculeSection = true;
            }
            else if (strcmp(token,"crowders")==0)
            {
                moleculeSection = false;
                crowdersSection = true;
            }
            else if (parameterSection)  // parameters
            {
                //gpus,streams,threads,mcsteps,resteps,replicas
                if (strcmp(token,"gpus")==0)
                {
                    token = strtok(NULL," ");
                    parameters->gpus = atoi(token);
                }
                if (strcmp(token,"streams")==0)
                {
                    token = strtok(NULL," ");
                    parameters->streams = atoi(token);
                }
                if (strcmp(token,"threads")==0)
                {
                    token = strtok(NULL," ");
                    parameters->threads = atoi(token);
                }
                if (strcmp(token,"mcsteps")==0)
                {
                    token = strtok(NULL," ");
                    parameters->MCsteps = atoi(token);
                }
                if (strcmp(token,"resteps")==0)
                {
                    token = strtok(NULL," ");
                    parameters->REsteps = atoi(token);
                }
                if (strcmp(token,"temperaturemax")==0)
                {
                    token = strtok(NULL," ");
                    parameters->temperatureMax = atof(token);
                }
                if (strcmp(token,"temperaturemin")==0)
                {
                    token = strtok(NULL," ");
                    parameters->temperatureMin = atof(token);
                }
                if (strcmp(token,"boundary")==0)
                {
                    token = strtok(NULL," ");
                    parameters->bound = atof(token);
                }
                if (strcmp(token,"replicas")==0)
                {
                    token = strtok(NULL," ");
                    parameters->replicas = atoi(token);
                }
                if (strcmp(token,"samplefrequency")==0)
                {
                    token = strtok(NULL," ");
                    parameters->sampleFrequency = atoi(token);
                }
                if (strcmp(token,"sampleafter")==0)
                {
                    token = strtok(NULL," ");
                    parameters->sampleStartsAfter = atoi(token);
                }
                if (strcmp(token,"prefix")==0)
                {
                    token = strtok(NULL," ");
                    strcpy(parameters->prependageString,"");
                    strcpy(parameters->prependageString,token);
#if REPULSIVE_CROWDING
                    strcat(parameters->prependageString,"_repulsive");
#else
                    strcat(parameters->prependageString,"_full");
#endif
                }
            }
            else if (moleculeSection||crowdersSection) // files
            {
                // TODO: read this into a vector of moldata; initialise replica afterwards
//                 float px,py,pz,rx,ry,rz,ra;
//                 char *pdbfilename = new char [256];
//                 bool translate = true;
                moldata m;
                int result = 0;
                if (saveline[0]=='t')
                {
                    result = sscanf(saveline,"t(%f,%f,%f) r(%f,%f,%f,%f) %s", &m.px,&m.py,&m.pz,&m.rx,&m.ry,&m.rz,&m.ra,m.pdbfilename);
                    m.translate = true;
                }
                if (saveline[0]=='p')
                {
                    result = sscanf(saveline,"p(%f,%f,%f) r(%f,%f,%f,%f) %s", &m.px,&m.py,&m.pz,&m.rx,&m.ry,&m.rz,&m.ra,m.pdbfilename);
                    m.translate = false;
                }
                if (result<8)
                {
                    cout << "Failed to parse line " << lc << ": " << saveline << endl;
                }
                else
                {
//                     int moleculeId = initialReplica->loadMolecule(pdbfilename);
//                     if (translate)
//                     {
//                         initialReplica->molecules[moleculeId].translate(Vector3f(px,py,pz));
//                     }
//                     else
//                     {
//                         initialReplica->molecules[moleculeId].setPosition(Vector3f(px,py,pz));
//
//                     }
//
//                     if (ra >= 0.000)
//                     {
//                         Vector3double v = Vector3double(rx,ry,rz);
//                         v.normalizeInPlace();
//                         initialReplica->molecules[moleculeId].rotate(v,ra);
//                     }
                    m.crowder = crowdersSection;
                    // TODO: remove this later
                    if (!crowdersSection)
                    {
                        parameters->nonCrowders++;
                    }
//                     if(crowdersSection)
//                     {
//                         initialReplica->molecules[moleculeId].setMoleculeRoleIdentifier(CROWDER_IDENTIFIER);
//                     }
//                     else
//                     {
//                         parameters->nonCrowders++;
//                         initialReplica->nonCrowderResidues += initialReplica->molecules[moleculeId].residueCount;
//                     }
                    parameters->mdata.push_back(m);
                }

//                 delete [] pdbfilename;
            }
        }
        input.close();
    }
    if (parameters->bound <= 0)
    {
        cout << "! Bounding value too small, setting equal to " << BOUNDING_VALUE << endl;
        parameters->bound = BOUNDING_VALUE;
    }

    if (parameters->temperatureMax < parameters->temperatureMin)
    {
        cout << "! Maximum temperature < minimum temperature, swapping " << parameters->temperatureMax << " <-> " << parameters->temperatureMin << endl;
        float tmp = parameters->temperatureMax;
        parameters->temperatureMax = parameters->temperatureMin;
        parameters->temperatureMin = tmp;
    }

    if (parameters->threads > parameters->replicas)
    {
        cout << "! Too many threads, setting equal to " << parameters->replicas << endl;
        parameters->threads = parameters->replicas;
    }
    if (parameters->threads > parameters->streams)
    {
        // TODO: is this right? Should it be set to the default threads?
        parameters->streams = parameters->threads;
        if (parameters->streams > 16*parameters->gpus)
        {
            parameters->streams = 16*parameters->gpus;
            if (parameters->streams > parameters->replicas)
                parameters->streams = parameters->replicas;
            cout << "! Too many streams, setting equal to " << parameters->gpus << endl;
        }
    }

//     initialReplica->nonCrowderCount = parameters->nonCrowders;
    char * fileindex = new char[256];
    sprintf(fileindex,"output/%s_%d_fileindex",parameters->prependageString,parameters->pid);
    FILE * fileindexf = fopen (fileindex,"w");
    fprintf(fileindexf,"index molecule_file_path crowder(Y/N)\n");
    delete [] fileindex;

// #if OUTPUT_LEVEL > 0
// TODO: move this outside, use initialReplica, remove nonCrowders from parameters; use new logging

//     cout << "Argument data from file:" << endl;
//     cout << "-------------------------------------" << endl;
//     cout << "threads " << parameters->threads << endl;
//     cout << "streams " << parameters->streams << endl;
//     cout << "GPUs " << parameters->gpus << endl;
//     cout << "mc steps " << parameters->MCsteps << endl;
//     cout << "re steps " << parameters->REsteps << endl;
//     cout << "replicas " << parameters->replicas << endl;
//     cout << "sampling frequency (mc steps) " << parameters->sampleFrequency << endl;
//     cout << "sampling starts after (mc steps) " << parameters->sampleStartsAfter << endl;
//     cout << "bounding box size " << parameters->bound << endl;
//     cout << "non-crowder molecules " << parameters->nonCrowders << endl;
//     cout << "maximum temperature " << parameters->temperatureMax << endl;
//     cout << "minimum temperature " << parameters->temperatureMin << endl;
//
//     cout << "Loaded: "<< endl;
//     cout << "-------------------------------------------------------------"<< endl;
//     for (int z=0; z<parameters->nonCrowders; z++)
//     {
//         printf("%2d %s centered @ (%0.3f,%0.3f,%0.3f)\n",z,initialReplica->molecules[z].filename,initialReplica->molecules[z].center.x,initialReplica->molecules[z].center.y,initialReplica->molecules[z].center.z);
//         fprintf(fileindexf,"%2d %s N\n",z,initialReplica->molecules[z].filename,initialReplica->molecules[z].center.x,initialReplica->molecules[z].center.y,initialReplica->molecules[z].center.z);
//     }
//     for (int z=parameters->nonCrowders; z<initialReplica->moleculeCount; z++)
//     {
//         printf("%2d %s crowder centered @ (%0.3f,%0.3f,%0.3f)\n",z,initialReplica->molecules[z].filename,initialReplica->molecules[z].center.x,initialReplica->molecules[z].center.y,initialReplica->molecules[z].center.z);
//         fprintf(fileindexf,"%2d %s Y\n",z,initialReplica->molecules[z].filename);
//     }
//     cout << "-------------------------------------------------------------"<< endl;
// #endif
    fclose(fileindexf);
}

int main(int argc, char **argv)
{
    cout.precision(8);
    int sysreturn;
    // TODO: move this inside Simulation init?
    sysreturn = system("mkdir -p checkpoints");
    sysreturn = system("mkdir -p output");
    cout << "Version: " << HGVERSION << endl;
    cout << "Compiled with:" << endl;
#ifdef GLVIS
    cout << "  OpenGL support" << endl;
#endif
#if USING_CUDA
    cout << "  CUDA support" << endl;
#endif
#if CUDA_STREAMS
    cout << "  Asynchronous GPU calls (CUDA capability 1.1+ required)" << endl;
#endif
#if COMPENSATE_KERNEL_SUM
    cout << "  Kahan summation in kernels" << endl;
#endif

    // get global options for the simulation
    argdata parameters;
    parameters.pid = int(getpid());
    parameters.inputFile = false;
    memset(parameters.prependageString,0,256);

    bool use_defaults = getArgs(&parameters, argc, argv);

    //if there is an input file load its contents here

    if (parameters.inputFile)
    {
        loadArgsFromFile(&parameters);
        getArgs(&parameters, argc, argv); // second pass to override any variables if doing performance tests
    }
    else
    {
        // deprecated function
        cout << "EXIT: Program requires a configuration file." << endl;
        cout << "use: REMCDockingGPU -f filename" << endl;
        exit(0);
    }

#if USING_CUDA
    if (strcmp(argv[argc-1],"-c")==0)
    {
        simulation.run_check();
    }
    else
    {
#endif
        simulation.init(parameters);

#if GLVIS
        if (parameters.viewConditions)
        {
            glutInit(&argc, argv);
            camera.setPosition(Vector3f(-15,15,15),Vector3f(1,-1,-1),Vector3f(0,1,0));
            char windowName[64] = {"REMC Protein Docker"};
            GlutInit(WIDTH,HEIGHT,windowName);
            gl_replicaCount = parameters.replicas;
            gl_boundingValue = int(parameters.bound);
            GLreplica = &simulation.initialReplica;
        }
#endif

        if (!parameters.skipsimulation)
        {
            simulation.run();
        }

#if GLVIS
        if (parameters.viewConditions)
        {
            // TODO: this is overwritten inside run. Why do we set it back here?
            GLreplica = &simulation.initialReplica;
            cout << "Entering free gl viewing mode." << endl;
            glutMainLoop();
        }
#endif

#if USING_CUDA
    }
#endif
    cout << "Finished." << endl;
    return 0;
}
