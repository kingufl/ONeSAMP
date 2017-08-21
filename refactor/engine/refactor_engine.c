// INCLUDE GUARD

#ifndef REFACTOR_ENGINE_C
#define REFACTOR_ENGINE_C

// HEADERS

#include "../macro/refactor_macro.h"
#include <stdio.h>
#include <stdlib.h>

// FUNCTIONS

/*! \brief Runs main engine for OneSamp.
 *
 */
int onesamp_engine(int argc, char **argv){

  // Read in command line arguments
  int syntax_results[2];
  resetArguments();
  parseArguments(argc, argv);

  // If checking syntax, run the parser and exit (argc == SYNTAX_ARGS)
  if(parseSyntaxCheck()){
    // Describes parameters of input data: numloci, individuals
    parseFromFile(TRUE, stdin, syntax_results);
    printf("-l%d -i%d\n", syntax_results[0], syntax_results[1]);
    return 0;
  }

  // Read in motifs
  if(parseFormFlag())
    if(readMicrosatelliteMotifLengths(argc, argv) != parseNLoci())
      reportError("Mangled list of motif lengths passed after -m flag: must be lists consisting of 2, 3, 4, and 6 equal to number of loci.");

  // Allocate all memory
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;
  int num_loci = parseNLoci();
  int initial_indivs_count = parseInputSamples();
  int final_indivs_count = parseInputSamples();

  // Allocate space to store results of statistics compuation
  allocateOneSampMemory(parseInputSamplesAllocation(), parseBottleneckMax(), parseInputSamples(), parseIterations(), parseNLociAllocation(), numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  double *mnals = doubleData;
  double *m = doubleData + 1 * parseIterations();
  double *iis = doubleData + 2 * parseIterations();
  double *lnbeta = doubleData + 3 * parseIterations();
  double *hetx = doubleData + 4 * parseIterations();
  double *mnehet = doubleData + 5 * parseIterations();
  double *mhomo = doubleData + 6 * parseIterations();
  double *varhomo = doubleData + 7 * parseIterations();
  double *skhomo = doubleData + 8 * parseIterations();
  double *kurhomo = doubleData + 9 * parseIterations();
  double *ne = doubleData + 10 * parseIterations();

  // Simulate the generations
  int i;
  int j;
  int k;
  int current;
  int next;
  int temp;

  // Allocate arrays to store intermediate generations.
  struct gtype_type *females[2], *males[2];
  for(i = 0; i < 2; i++){
    females[i] = (struct gtype_type *)malloc(parseBottleneckMax() * STRUCT_GTYPE_SIZE);
    males[i] = (struct gtype_type *)malloc(parseBottleneckMax() * STRUCT_GTYPE_SIZE);
    for(j = 0; j < parseBottleneckMax(); j++){
      females[i][j].pgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
      males[i][j].pgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
      females[i][j].mgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
      males[i][j].mgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
    }
  }

  if(!parseExamplePop()){
    parseFromFile(FALSE, stdin, syntax_results);
    // Read in the data from a file.
    // Parse input file for initial condition
    if(syntax_results[0] != parseNLoci())
      reportError("Supplied number of loci on command line is incorrect.");
    if(syntax_results[1] != parseInputSamples())
      reportError("Supplied number of individuals on command line is incorrect.");

    // If simulating one generation as specified by flag, immediately simulate
    // and return results
  }

  if(parseSingleGeneration()){
    // Copy SNPs over to final generation unchanged if we're computing
    // stats without simulating populations
  final_indivs_data = (struct gtype_type **)malloc(parseIterations() * STRUCT_GTYPE_STAR_SIZE);
    for(i = 0; i < 1; i++) {
      final_indivs_data[i] = (struct gtype_type *)malloc(2 * parseBottleneck(0) * STRUCT_GTYPE_SIZE);
      for(j = 0; j < 2 * parseBottleneck(0); j++){
        final_indivs_data[i][j].pgtype = (ALLELE_TYPE *)malloc(parseNLoci() * sizeof(ALLELE_TYPE));
        final_indivs_data[i][j].mgtype = (ALLELE_TYPE *)malloc(parseNLoci() * sizeof(ALLELE_TYPE));
      }
    }
    // We are doing just one iteration, so pick that value
    assort(2 * parseBottleneck(0), final_indivs_data[0], initial_indivs_data, initial_indivs_data, parseInputSamples(), 0, num_loci);
    writeoutput(final_indivs_data, 2 * parseBottleneck(0));
    for(i = 0; i < 1; i++) {
      for(j = 0; j < 2 * parseBottleneck(0); j++){
        free(final_indivs_data[i][j].pgtype);
        free(final_indivs_data[i][j].mgtype);
      }
      free(final_indivs_data[i]);
    }
    free(final_indivs_data);
    return 0;
  }

  if(!parseExamplePop()){
    filterMonomorphicLoci();
    filterLowCoverageLoci();
    filterLowCoverageIndividuals();
    if(parseFillInAbsentData()) fillInMissingData();
    //printf("individuals = %d, loci = %d", parseInputSamples(), parseNLoci());
    //fflush(stdout);
    //exit(1);
    //obtainMissingProportion();
  }

/*! \brief Allocates arrays to store sampling generations.
 */
  final_indivs_data = (struct gtype_type **)malloc(parseIterations() * STRUCT_GTYPE_STAR_SIZE);
  for(i = 0; i < parseIterations(); i++) {
    final_indivs_data[i]=(struct gtype_type *)malloc(parseInputSamples() * STRUCT_GTYPE_SIZE);
    for(j = 0; j < parseInputSamples(); j++){
      final_indivs_data[i][j].pgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci() * sizeof(ALLELE_TYPE));
      final_indivs_data[i][j].mgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci() * sizeof(ALLELE_TYPE));
    }
  }

  // Simulate intermediate bottleneck generations
  // Make several data sets (parseIterations() of them) of final generations.
  if(!parseRawSample()) for(i = 0; i < parseIterations(); i++){
    // If the flag is set, make up a sample data set to complete an example.
    // Generate a simulated case of initial conditions
    #include "../engine/refactor_coalescent_engine.txt"

    current = 0;
    next = 1;

    // Simulate bottleneck generations from random mating
    for(j = 0; j < parseBottleneckLength(i); j++){
      assort(parseBottleneck(i), females[next], females[current], males[current], parseBottleneck(i), i, extraProportionOfBufferLoci * parseNLoci());
      assort(parseBottleneck(i), males[next], females[current], males[current], parseBottleneck(i), i, extraProportionOfBufferLoci * parseNLoci());
      temp = current;
      current = next;
      next = temp;
    }

    // Then, generate a set of genotypes for final generation
    assort(final_indivs_count, final_indivs_data[i], females[current], males[current], parseBottleneck(i), i, extraProportionOfBufferLoci * parseNLoci());

    // Remove monomorphic loci if possible by replacing them with the extra loci
    int extraLociIndex = parseNLoci();
    for(j = 0; j < parseNLoci(); j++){
      int tempGene1;
      int tempGene2;
      if(variantsOfFinalLocus(i, j) == 2) continue;
      while(extraLociIndex < extraProportionOfBufferLoci * parseNLoci() && variantsOfFinalLocus(i, extraLociIndex) < 2){
        extraLociIndex++;
      }
      if(!(extraLociIndex < extraProportionOfBufferLoci * parseNLoci())) break;
      for(k = 0; k < parseInputSamples(); k++){
        loadFinalGenotype(i, k, extraLociIndex, &tempGene1, &tempGene2);
        storeFinalGenotype(i, k, j, &tempGene1, &tempGene2);
      }
      extraLociIndex++;
    }

    // Add in missing data in coalescent population to mimic input population.
    //double missingDataProbability = getProportionMissingData();

    if(!parseExamplePop()){
      for(k = 0; k < parseInputSamples(); k++){
        int maskFromInputSample = disrand(0, parseInputSamples() - 1);
        for(j = 0; j < parseNLoci(); j++){
          int mother;
          int father;
          loadInitialGenotype(maskFromInputSample, j, &father, &mother);
          if(father == 0) final_indivs_data[i][k].pgtype[j] = 0;
          if(mother == 0) final_indivs_data[i][k].mgtype[j] = 0;
        }
      }
    }

    //writeoutput(final_indivs_data, final_indivs_count);

  } else { // We are using a raw sample
    // Copy SNPs over to final generation unchanged if we're computing
    // stats without simulating populations

    for(i = 0; i < parseInputSamples(); i++){
      for(j = 0; j < parseNLoci(); j++){
        final_indivs_data[0][i].pgtype[j] = initial_indivs_data[i].pgtype[j];
        final_indivs_data[0][i].mgtype[j] = initial_indivs_data[i].mgtype[j];
      }
    }
  }

  // Deallocate intermediate genotype arrays
  for(i = 0; i < 2; i++){
    for(j = 0; j < parseBottleneckMax(); j++){
      free(females[i][j].pgtype);
      free(males[i][j].pgtype);
      free(females[i][j].mgtype);
      free(males[i][j].mgtype);
    }
    free(females[i]);
    free(males[i]);
  }

  if(parseExamplePop()){
    // Dump population
    writeoutput(final_indivs_data, parseInputSamples());
  } else {
    // Compute statistics

    //writeoutput(final_indivs_data, parseInputSamples());

    // Statistic 6: mnals
    // Summarize information about alleles
    counts(numberOfAllelesPtr, final_indivs_data, mnals, gType, gcountPtr);

    // Statistic 1: m
    // Only do next call if we're not using SNPs
    // Calculate range, m, change in frequency of alleles 
    sortM(numberOfAlleles, parseIterations(), m, gType, gcountPtr);

    // Statistic 3: lnbeta
    // Only do next call if we're not using SNPs
    // Calculate beta statistic
    beta(numberOfAlleles, parseIterations(), lnbeta, gType, gcountPtr);

    // Statistics 4 and 5: hetx, mnehet
    // Calculate the excess heterozygosity
    hetexcess(numberOfAlleles, parseIterations(), final_indivs_data, hetx, mnehet, gType, gcountPtr);

    // Statistics 7 and 8 (and 9 and 10): mhomo, varhomo, skhomo, kurhomo
    // Calculate mean, variance, skew, and kurtosis of heterozygosity
    multih(parseIterations(), final_indivs_data, mhomo, varhomo, skhomo, kurhomo, gType);

    // Statistic 2: iis
    // Calculate Burrows Weir stat from Vitalis and Couvet
    twolocusiis(numberOfAlleles, parseIterations(), final_indivs_data, iis, gType, gcountPtr);

    // Statistic 0: ne
    // Calculate Ne
    for(i = 0; i < parseIterations(); i++){
      ne[i] = parseRawSample() ? -1 : (2*parseBottleneck(i)+(double)(1.0/(2*parseBottleneck(i)))+0.5);
    }

    for(i = 0; i < parseIterations(); i++){
      printf("%f %f %f %f %f %f %f %f %f\n", ne[i], iis[i], hetx[i], mnehet[i], mnals[i], mhomo[i], varhomo[i], m[i], lnbeta[i]);
    }

  }

  // Deallocate structure 5
  for(i = 0; i < parseIterations(); i++){
    for(j = 0; j < parseInputSamples(); j++){
      free(final_indivs_data[i][j].pgtype);
      free(final_indivs_data[i][j].mgtype);
    }
    free(final_indivs_data[i]);
  }
  free(final_indivs_data);
  deallocateOneSampMemory(initial_indivs_count, parseBottleneckMax(), final_indivs_count, parseIterations(), num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);
  flushArguments();

  // Stop the random number table.
  closegfsr();
  return 0;
}

#endif
