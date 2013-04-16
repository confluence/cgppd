#ifndef MAIN_H_
#define MAIN_H_
#include "Simulation.h"

using namespace std;

// Simulation simulation;

void printHelp(bool badArg); // TODO: wrapper for bad args + help; rename to usage
bool getArgs(argdata * d, int argc, char **argv);
void loadArgsFromFile(argdata * parameters);

#endif // MAIN_H_