#ifndef REFACTOR_PARSER_H
#define REFACTOR_PARSER_H

// Parser functions
void parse(int syntax_check_flag, const char *dataFileName, int *syntax_results);
void parseFromFile(int syntax_check_flag, FILE *dataFileName, int *syntax_results);
int skipPastWhitespace();
int skipPastDigits();
int skipPastComma();
void enqueueParserToken(int token);
void parseGenotypeToken(int *returnGene1, int *returnGene2);

// Parser queue functions
void initializeParserQueue();
int *accessPriorElement(int ago);

// String match functions
int matchPop();

// Match functions
int matchSlashR(int *c);
int matchSlashN(int *c);
int matchComma(int *c);
int matchEOF(int *c);
int matchP(int *c);
int matchO(int *c);
int matchDigit(int *c);
int matchEndOfLine(int *c);
int matchWhitespace(int *c);
int matchHardDelimiter(int *c);
int matchSoftDelimiter(int *c);

// Logical function
int zeroOne(int i);

// Input filtering functions
void filterLowCoverageLoci();
void filterLowCoverageIndividuals();
void fillInMissingData();
void obtainMissingProportion();
void filterMonomorphicLoci();

// Determine variants
int variantsOfFinalLocus(int sample, int locusI); 

#endif
