#ifndef REFACTOR_PARSER_C
#define REFACTOR_PARSER_C

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "../macro/refactor_macro.h"

#define EOF_ERROR_ESCAPE() reportParseError((char *) eof_error_string, *lptr, *cptr)
#define COMMA_ERROR_ESCAPE() reportParseError((char *) comma_error_string, *lptr, *cptr)

#define QUEUE_LENGTH 8
#define MAX_NO_LINES 1000000
#define MAX_NO_CHARACTERS_PER_LINE 1000000
#define ACCEPTCHARACTER() nextChar = fgetc(dataFile); enqueueParserToken(nextChar);

char eof_error_string[] = "Unexpected EOF";
char comma_error_string[] = "Unexpected comma";
char genePop_error_string[] = "Genotype data does not conform to GenePop 4.0 standard: malformed number.";
char blankLine_error_string[] = "GenePop 4.0 format does not allow blank lines.";
char unexpected_input_string[] = "Unexpected data in input file: incorrect number of individuals or loci specified to program, or inconsistent number of alleles specified for each individual.";
char invalid_SNP_string[] = "Invalid SNP: must be a number that is 0, 1, 2, 3, or 4.";

long lineNum; // Which line we are currently at in the file.
long charNum; // Which character we are currently at in the file past last line.
long *lptr; // Pointer to lineNum variable.
long *cptr; // Pointer to charNum variable.
int nextChar; // Next character in file.
int queue[QUEUE_LENGTH]; // Data stored in queue
int queueIndex; // Index of most recent element.
int assertEOF; // Detect if the next token should be EOF: (1 if so, 0 if not)
int parseNumGenes; // Detect the number of genes per line.
int individual; // Detect the number of individuals in the file.
FILE *dataFile; // Variable to store index in file.

const char modeRead[] = "r\0";

/*! \brief Parses data from an input filename
 */
void parse(int syntax_check_flag, const char *dataFileName, int *syntax_results){
  // Read in data file and parse
  FILE *dataFile = fopen(dataFileName, modeRead);
  parseFromFile(syntax_check_flag, dataFile, syntax_results);
  fclose(dataFile);
}

/*! \brief Main parsing engine
 */
void parseFromFile(int syntax_check_flag, FILE *curData, int *syntax_results){
  dataFile = curData;
  initializeParserQueue();
  nextChar = fgetc(dataFile);

  syntax_results[0] = 0; // Number of loci per sample
  syntax_results[1] = 0; // Number of individuals

  // Parse first line (skip)
  while(!matchEndOfLine(&nextChar) && !matchEOF(&nextChar)) nextChar = fgetc(dataFile);

  assert(*lptr == 1);

  (*lptr)++; (*cptr) = 0;

  // Check for unexpected EOF.
  if(matchEOF(&nextChar)) EOF_ERROR_ESCAPE();

  assert(*lptr == 2);

  // Parse second line, count commas as proxy for number of genes, but if there are no commas, count identifiers until POP.

  parseNumGenes = 0;
  while(!matchPop()){
    if(matchEOF(&nextChar)) EOF_ERROR_ESCAPE();   // Check for unexpected EOF.
    if(matchComma(&nextChar) || matchEndOfLine(&nextChar)) parseNumGenes++;
    ACCEPTCHARACTER();
  }
  parseNumGenes--;

  // Ensure we have a reasonable number of genes.
  assert(parseNumGenes > 0);
  //assert(parseNumGenes < 100000);

  // Process genotype information
  int locusCounter;
  int returnGene1;
  int returnGene2;
  int s;
  
  while(!matchEOF(&nextChar)){
    while(!matchEndOfLine(&nextChar) && !matchEOF(&nextChar)){
      skipPastComma();
      for(locusCounter = 0; locusCounter < parseNumGenes || matchEOF(&nextChar); locusCounter++){
        if(nextChar == ' ' || nextChar == '\t') {ACCEPTCHARACTER();}
        if(matchEndOfLine(&nextChar)) break;
        while(matchSoftDelimiter(&nextChar) && !matchEOF(&nextChar)){ ACCEPTCHARACTER();}
        ACCEPTCHARACTER();
        skipPastDigits();
        parseGenotypeToken(&returnGene1, &returnGene2);
        //printf("%d %d\n", returnGene1, locusCounter);
        if((returnGene1 > 4 || returnGene2 > 4) && parseFormFlag() == 0)
          reportParseError(invalid_SNP_string, *lptr, *cptr - 1);
        if(matchEOF(&nextChar)) break;
//        if(individual >= parseInputSamples() || locusCounter >= parseNLoci())
//          reportParseError(unexpected_input_string, *lptr, *cptr);
        if(!syntax_check_flag) {
          storeInitialGenotype(individual, locusCounter, &returnGene1, &returnGene2);
        }
      }
      individual++;
      while(matchWhitespace(&nextChar)) {ACCEPTCHARACTER();}
    }
    syntax_results[1] = individual;
    individual = 0;
    ACCEPTCHARACTER();
  }
  syntax_results[0] = parseNumGenes;
  return;
}

/*! \brief Skips the parser beyond numeric digits.
 */
int skipPastDigits(){
  while(matchDigit(&nextChar) && !matchEOF(&nextChar)){ ACCEPTCHARACTER();}
  return 0;
}

/*! \brief Skips the parser beyond intervening whitespace.
 */
int skipPastWhitespace(){
  while(matchWhitespace(&nextChar) && !matchEOF(&nextChar)){ ACCEPTCHARACTER();}
  return 0;
}

/*! \brief Skips the parser beyond the next comma.
 */
int skipPastComma(){
  while(!matchComma(&nextChar) && !matchEOF(&nextChar)){ ACCEPTCHARACTER();}
  return 0;
}

/////////////////////////////////////
// PARSER STATE QUEUE AND LISTENER //
/////////////////////////////////////

// Define a queue to store the last seven tokens from the parser.
// The current queue length was chosen to correctly identify the token
// XPOPX, where X is a delimiter: \r \n comma or EOF, and nothing is case
// sensitive.

// For example, if the queueIndex is 3 and it reads in ABCDEFG, the queue
// function accessPriorElement will return G, F, E, D, C

// QUEUE_LENGTH should be long enough to facilitate the parser handling these
// operations:

// 1. Detect POP delimiter (5 characters)
// 2. Read in a 4 digit genotype and detect a soft delimiter on both ends
//    (6 characters)
// 3. Read in a 6 digit genotype and detect a soft delimiter on both ends
//    (8 characters)

/*! \brief Set up a queue to detect the state of the parser.
 */
void initializeParserQueue(){
  queueIndex = 0;
  assertEOF = 0;
  parseNumGenes = 0;

  // Zero out queue to EOF tokens
  int i;
  for(i = 0; i < QUEUE_LENGTH; i++) queue[i] = EOF;

  // Keep track of where we are in the file: line and column numbers.
  lineNum = 1;
  charNum = 0;
  lptr = &lineNum;
  cptr = &charNum;
}

/*! \brief Places a character into the parser queue
 */
void enqueueParserToken(int token){

//  printf("%c\n", token); fflush(stdout);

  nextChar = token;
  // Place the next element in the queue
  queueIndex = (queueIndex + 1) % QUEUE_LENGTH;
  queue[queueIndex] = token;

  // Prior patterns in the code might make us expect an end of file to comply
  // with GenePop 4.0. Enforce these restrictions here.

  if(assertEOF && token != EOF)
    reportParseError(blankLine_error_string, *lptr, *cptr);

  // Add 1 to character position.
  charNum++;

  // Restrict length of line.
  //assert(charNum < MAX_NO_CHARACTERS_PER_LINE);

  int prior1 = *accessPriorElement(1);
  int prior2 = *accessPriorElement(2);

  // Compare the prior two elements with end of line characters.
  int endOfLinePrior1 = matchEndOfLine(&prior1);
  int endOfLinePrior2 = matchEndOfLine(&prior2);

  // If we don't have an end of line, we're done.
  if(!matchEndOfLine(&token)) return;

  charNum = 0;
  (*lptr)++;

  // If this point is reached, must do end of line tasks.
  // But, don't duplicate them if the last two characters are \r and \n.
  if(endOfLinePrior1 && prior1 != token){
    charNum = 0;
    (*lptr)++;
  }

  // Don't want very large numbers of lines.
  assert(*lptr < MAX_NO_LINES);

  // A set of three consecutive matches to \r or \n in a row is not allowed
  // in GenePop 4.0 unless we are at the end of the file.
  // The parser will enforce this condition.
  if(endOfLinePrior1 && endOfLinePrior2) {assertEOF = 1; lineNum;}

  // A consecutive pair of \n characters or \r characters is not allowed in GenePop 4.0 unless we are at the end of the file.
  // The parser will enforce this condition.
  if(prior1 == token) {assertEOF = 1; lineNum--;}
}

/*! \brief Determines the element in the queue how many steps ago (0 is most recent).
 */
int *accessPriorElement(int ago){
  assert(ago >= 0);
  assert(ago < QUEUE_LENGTH);
  return queue + ((queueIndex + (QUEUE_LENGTH - ago)) % QUEUE_LENGTH);
}

/*! \brief Reads numeric genotypes from the parser queue
 */
void parseGenotypeToken(int *returnGene1, int *returnGene2){
  int *prior7 = accessPriorElement(7);
  int *prior6 = accessPriorElement(6);
  int *prior5 = accessPriorElement(5);
  int *prior4 = accessPriorElement(4);
  int *prior3 = accessPriorElement(3);
  int *prior2 = accessPriorElement(2);
  int *prior1 = accessPriorElement(1);
  int *prior0 = accessPriorElement(0);

  // The fifth to last character determines whether the number has 4 or 6 digits.
  if(!matchDigit(prior5) && !matchSoftDelimiter(prior5)){
    reportParseError(genePop_error_string, *lptr, *cptr - 5);
  }

  // Check for legal 4 digit number
  if(matchSoftDelimiter(prior0) && matchDigit(prior1) && matchDigit(prior2)
            && matchDigit(prior3) && matchDigit(prior4)
            && matchSoftDelimiter(prior5)){

    *returnGene1 = 10 * (*prior4 - ASCII_ZERO) + (*prior3 - ASCII_ZERO);
    *returnGene2 = 10 * (*prior2 - ASCII_ZERO) + (*prior1 - ASCII_ZERO);
    return; // Successfully stored the genotype pair.
  }

  // Check for legal 6-digit number
  if(matchSoftDelimiter(prior0) && matchDigit(prior1) && matchDigit(prior2)
          && matchDigit(prior3) && matchDigit(prior4) && matchDigit(prior5)
          && matchDigit(prior6) && matchSoftDelimiter(prior7)){
    *returnGene1 = 100 * (*prior6 - ASCII_ZERO) + 10 * (*prior5 - ASCII_ZERO)
                      + (*prior4 - ASCII_ZERO);
    *returnGene2 = 100 * (*prior3 - ASCII_ZERO) + 10 * (*prior2 - ASCII_ZERO)
                      + (*prior1 - ASCII_ZERO);
    return; // Successfully stored the genotype pair.
  }

  // Check for EOF
  if(matchEOF(prior0) || matchEOF(prior1) || matchEOF(prior2) || matchEOF(prior3) || matchEOF(prior4) || matchEOF(prior5) || matchEOF(prior6) || matchEOF(prior7)) return;

  // Otherwise, data is malformed and report error
  reportParseError(genePop_error_string, *lptr, *cptr - 1);
}

/*! \brief Determines if the letters POP were read in to denote the beginning of the genotypes
 */
int matchPop(){
  return zeroOne(matchEndOfLine(accessPriorElement(0)) &&
                 matchP(accessPriorElement(1)) &&
                 matchO(accessPriorElement(2)) &&
                 matchP(accessPriorElement(3)) &&
                 matchEndOfLine(accessPriorElement(4)));
}

/*! \brief Matches carriage return
 */
int matchSlashR(int *c){
  return zeroOne(((char) *c) == '\r');
}

/*! \brief Matches newline character
 */
int matchSlashN(int *c){
  return zeroOne(((char) *c) == '\n');
}

/*! \brief Matches a comma
 */
int matchComma(int *c){
  return zeroOne(((char) *c) == ',');
}

/*! \brief Matches an end of file
 */
int matchEOF(int *c){
  return zeroOne(*c == EOF);
}

/*! \brief Matches the letter P
 */
int matchP(int *c){
  return zeroOne(((char) *c == 'P') || ((char) *c) == 'p');
}

/*! \brief Matches the letter O
 */
int matchO(int *c){
  return zeroOne(((char) *c) == 'O' || ((char) *c) == 'o');
}

/*! \brief Matches a digit
 */
int matchDigit(int *d){
  return zeroOne((char)*d == '0' || (char)*d == '1' || (char)*d == '2' ||
                 (char)*d == '3' || (char)*d == '4' || (char)*d == '5' ||
                 (char)*d == '6' || (char)*d == '7' || (char)*d == '8' ||
                 (char)*d == '9');
}

/*! \brief Matches an end of line character
 */
int matchEndOfLine(int *c){
  return zeroOne(matchSlashN(c) || matchSlashR(c));
}

/*! \brief Matches whitespace
 */
int matchWhitespace(int *c){
  return zeroOne(((char) *c) == ' '
              || ((char) *c) == '\t'
              || matchEndOfLine(c));
}

/*! \brief Matches a break between lines or other explicit delimiter such as a comma
 */
int matchHardDelimiter(int *c){
  return zeroOne(matchEndOfLine(c) || matchComma(c) || matchEOF(c));
}

/*! \brief Matches a break between tokens in the parser
 */
int matchSoftDelimiter(int *c){
  return zeroOne(matchHardDelimiter(c) || matchWhitespace(c));
}

/*! \brief Converts a boolean value to a zero or a one.
 */
int zeroOne(int i){
  return i == 0 ? 0 : 1;
}

/*! \brief Returns Hamming distance between two pairs of genes.
 */
int pairDistance(int geneA1, int geneA2, int geneB1, int geneB2){
  if(geneA1 > geneA2){
    return pairDistance(geneA2, geneA1, geneB1, geneB2);
  }
  if(geneB1 > geneB2){
    return pairDistance(geneA1, geneA2, geneB2, geneB1);
  }
  if(geneA2 == 0) return 2;
  if(geneB2 == 0) return 2;
  if(geneA1 == 0) return (geneA2 == geneB1 || geneA2 == geneB2) ? 1 : 2;
  if(geneB1 == 0) return (geneB2 == geneA1 || geneB2 == geneA2) ? 1 : 2;
  return (geneA1 == geneB1 ? 0 : 1) + (geneA2 == geneB2 ? 0 : 1);
}

/*! \brief Computes the Hamming distance between two individuals.
 */
int hammingDistance(int indivI, int indivJ){
  int geneA1;
  int geneA2;
  int geneB1;
  int geneB2;
  int result = 0;
  int num_loci = parseNLoci();
  int k;
  for(k = 0; k < num_loci; k++){
    loadInitialGenotype(indivI, k, &geneA1, &geneA2);
    loadInitialGenotype(indivJ, k, &geneB1, &geneB2);
    result += pairDistance(geneA1, geneA2, geneB1, geneB2);
  }
  return result;
}

#define DISTANCE(i, j) distances[(i) * numIndivs + (j)]

/*! \brief Fills in data with its nearest neighbor.
 */
void fillInMissingData(){
  int geneA1;
  int geneA2;
  int geneB1;
  int geneB2;
  int i, j, k, l;
  int num_loci = parseNLoci();
  int numIndivs = parseInputSamples();
  int *distances = malloc(numIndivs * numIndivs * sizeof(int));

  // Collect Hamming distances between all pairs of individuals
  for(i = 0; i < numIndivs; i++){
    for(j = 0; j < numIndivs; j++){
      DISTANCE(i, j) = hammingDistance(i, j);
    }
  }

  int smallestDistance;

  // Replace data with nearest neighbor according to Hamming distance
  for(i = 0; i < numIndivs; i++){
    for(j = 0; j < num_loci; j++){
      loadInitialGenotype(i, j, &geneA1, &geneA2);
      if(geneA1 == 0 || geneA2 == 0){

        // Fill in individual i, locus j
        smallestDistance = 2 * num_loci;
        geneB1 = 0;
        geneB2 = 0;
        // Loop through candidate loci to find smallest distance
        for(k = 0; k < numIndivs; k++){
          if(i == k) continue;
          loadInitialGenotype(k, j, &geneB1, &geneB2);
          if(geneB1 == 0 || geneB2 == 0) continue;
          if(DISTANCE(i, k) < smallestDistance){
            smallestDistance = DISTANCE(i, k);
          }
        }
        // Select most frequently occurring nearest neighbor
        int *frequencies = malloc(numIndivs * sizeof(int));
        for(k = 0; k < numIndivs; k++){
          frequencies[k] = 0;
        }
        for(k = 0; k < numIndivs; k++){
          if(DISTANCE(i, k) == smallestDistance){
            for(l = 0; l <= k; l++){
              if(pairDistance(geneA1, geneA2, geneB1, geneB2) == 0){
                frequencies[l]++;
                break;
              }
            }
          }
        }
        // Locate maximum frequency
        int maxIndex = 0;
        for(k = 0; k < numIndivs; k++){
          if(frequencies[maxIndex] < frequencies[k]){
            maxIndex = k;
          }
        }

        // Fill in missing data
        loadInitialGenotype(maxIndex, j, &geneB1, &geneB2);
        storeInitialGenotype(maxIndex, j, &geneB1, &geneB2);

        free(frequencies);
        // Select most common locus
      }
    }
  }
  free(distances);
}

/*! \brief Removes a locus from consideration.
 */
void deleteLocus(int locusI, int num_loci){
  int geneA1;
  int geneA2;
  int numIndivs = parseInputSamples();
  int i;
  for(i = 0; i < numIndivs; i++){
    loadInitialGenotype(i, num_loci - 1, &geneA1, &geneA2);
    storeInitialGenotype(i, locusI, &geneA1, &geneA2);
  }
  // Fix repeat motif lengths
  if(parseFormFlag())
    getMotifLengths()[locusI] = getMotifLengths()[num_loci - 1];
}

/*! \brief Removes an individual from the simulation.
 */
void deleteIndividual(int indivI, int numIndivs){
  int geneA1;
  int geneA2;
  int num_loci = parseNLoci();
  int i;
  for(i = 0; i < num_loci; i++){
    loadInitialGenotype(numIndivs - 1, i, &geneA1, &geneA2);
    storeInitialGenotype(indivI, i, &geneA1, &geneA2);
  }
}

/*! \brief Computes the coverage of a locus.
 */
int locusCoverage(int locusI){
  int geneA1;
  int geneA2;
  int totalIndividualsWithoutMissingData = 0;
  int i;
  for(i = 0; i < parseInputSamples(); i++){
    loadInitialGenotype(i, locusI, &geneA1, &geneA2);
    if(geneA1 != 0 && geneA2 != 0) totalIndividualsWithoutMissingData++;
  }
  return totalIndividualsWithoutMissingData;
}

/*! \brief Computes the number of variants of a locus (returns zero, one, or two+)
 * If the number of alleles is at least 2, this method returns 2
 */
int variantsOfLocus(int locusI){
  int geneA1;
  int geneA2;
  int nonZeroVariant = 0;
  int i;
  for(i = 0; i < parseInputSamples(); i++){
    loadInitialGenotype(i, locusI, &geneA1, &geneA2);
    if(nonZeroVariant == 0) nonZeroVariant = geneA1;
    if(nonZeroVariant == 0) nonZeroVariant = geneA2;
    if(nonZeroVariant == 0) continue;
    if(nonZeroVariant != geneA1 || nonZeroVariant != geneA2) return 2;
  }
  if(nonZeroVariant == 0) return 0;
  return 1;
}

/*! \brief Computes the coverage of an individual.
 */
int indivCoverage(int indivI){
  int geneA1;
  int geneA2;
  int totalLociWithoutMissingData = 0;
  int i;
  for(i = 0; i < parseNLoci(); i++){
    loadInitialGenotype(indivI, i, &geneA1, &geneA2);
    if(geneA1 != 0 && geneA2 != 0) totalLociWithoutMissingData++;
  }
  return totalLociWithoutMissingData;
}

/*! \brief Removes loci with coverage below a threshold from consideration.
 */
void filterLowCoverageLoci(){
  int numIndivs = parseInputSamples();
  int num_loci = parseNLoci();
  int i = 0;
  while(i < num_loci){
    if(locusCoverage(i) < parseOmitLocusThreshold() * numIndivs){
      deleteLocus(i, num_loci);
      num_loci--;
    } else i++;
  }
  setNLoci(num_loci);
}

/*! \brief Removes individuals with coverage below a threshold from consideration.
 */
void filterLowCoverageIndividuals(){
  int numIndivs = parseInputSamples();
  int num_loci = parseNLoci();
  int i = 0;
  while(i < numIndivs){
    if(indivCoverage(i) < parseOmitLocusThreshold() * num_loci){
      deleteIndividual(i, numIndivs);
      numIndivs--;
    } else i++;
  }
  setNLoci(num_loci);
}

/*! \brief Removes loci that are monomorphic.
 */
void filterMonomorphicLoci(){
  int numIndivs = parseInputSamples();
  int num_loci = parseNLoci();
  int i = 0;
  while(i < num_loci){
    if(variantsOfLocus(i) < 2){
      deleteLocus(i, num_loci);
      num_loci--;
    } else i++;
  }
  setNLoci(num_loci);
}

/*! \brief Sets the proportion of missing data.
 */
void obtainMissingProportion(){
  int numIndivs = parseInputSamples();
  int num_loci = parseNLoci();
  int totalCoverage = 0;
  int j = 0;
  while(j < numIndivs){
    totalCoverage += indivCoverage(j);
    j++;
  }
  setProportionMissingData(((double) totalCoverage) / (num_loci * numIndivs));
}

/*! \brief Computes the number of variants of a locus (returns zero, one, or two+)
 * If the number of alleles is at least 2, this method returns 2
 */
int variantsOfFinalLocus(int sample, int locusI){
  int geneA1;
  int geneA2;
  int nonZeroVariant = 0;
  int i;
  for(i = 0; i < parseInputSamples(); i++){
    loadFinalGenotype(sample, i, locusI, &geneA1, &geneA2);
    if(nonZeroVariant == 0) nonZeroVariant = geneA1;
    if(nonZeroVariant == 0) nonZeroVariant = geneA2;
    if(nonZeroVariant == 0) continue;
    if(nonZeroVariant != geneA1 || nonZeroVariant != geneA2) return 2;
  }
  if(nonZeroVariant == 0) return 0;
  return 1;
}

#endif
