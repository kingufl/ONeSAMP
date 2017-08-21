/*
 * Header file for refactored version of ONeSAMP that handles command line
 * arguments. See the implementation for more details.
 */

// INCLUDE GUARD

#ifndef REFACTOR_ARGUMENTS_H
#define REFACTOR_ARGUMENTS_H

// HEADERS

#include <stdlib.h>
#include <stdio.h>
#include "../macro/refactor_macro.h"

// PROTOTYPES

void reportParseError(char *message, long row, long column);
void reportError(char *message);
void resetArguments();
void parseArguments(int argc, char **argv);
void flushArguments();
char *parseProgramName();
int parseFillInAbsentData();
int parseFormFlag();
int parseIterations();
int parseNLoci();
void setNLoci();
int parseNLociAllocation();
int parseInputSamples();
void setInputSamples();
int parseInputSamplesAllocation();
int parseBottleneck(int samp);
int parseBottleneckMin();
int parseBottleneckMax();
int parseBottleneckLength(int samp);
int parseBottleneckLengthMin();
int parseBottleneckLengthMax();
int parseFinalSize();
int parseIndividuals();
double parseMinAlleleFrequency();
double parseOmitLocusThreshold();
double parseTheta(int samp);
double parseThetaMin();
double parseThetaMax();
double parseMRate(int samp);
double parseMRateMin();
double parseMRateMax();
int parseRFlag();
int parseSyntaxCheck();
int parseExample();
int parseExamplePop();
int parseRawSample();
int parseSingleGeneration();
double parsePositiveDouble(int argNum, char **argv);
double *parsePositiveDoublePair(int argNum, char **argv);
int parsePositiveInt(int argNum, char **argv);
int *parsePositiveIntPair(int argNum, char **argv);
double getProportionMissingData();
void setProportionMissingData(double val);
int *getMotifLengths();
void setMotifLengths(int *new_motif_lengths);
void freeMotifLengths();
int initializeMicrosat1(int locus);
int initializeMicrosat2(int locus);
int readMicrosatelliteMotifLengths(int argc, char **argv);
// INCLUDE GUARD

#endif
