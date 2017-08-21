#include <gtest/gtest.h>

extern "C"{
#include "../macro/refactor_macro.h"
}

#define CRASH assert(0);

TEST(stats, SNPtest1){

  int i, j;
  int argc = 12;
  char a0[] = {'o', 'n', 'e', 's', 'a', 'm', 'p', '\0'};
  char a1[] = {'-', 'l', '6', '\0'};
  char a2[] = {'-', 'i', '3', '\0'};
  char a3[] = {'-', 's', '\0'};
  char a4[] = {'-', 't', '1', '\0'};
  char a5[] = {'-', 'b', '8', '\0'};
  char a6[] = {'-', 'w', '\0'};
  char a7[] = {'-', 'r', 'C', '\0'};
  char a8[] = {'-', 'd', '0', '\0'};
  char a9[] = {'-', 'v', '1', '\0'};
  char a10[] = {'-', 'u', '0', '.', '5', '\0'};
  char a11[] = {'-', 'o', '0', '\0'};
  char *argv[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11};
  int initial_indivs_count = 0;
  int bottleneck_indivs_count = 0;
  int final_indivs_count = 3;
  int num_samples = 3;
  int num_loci = 6;
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;

  parseArguments(argc, argv);
  allocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

/*! \brief Allocates arrays to store sampling generations.
 */
  final_indivs_data = (struct gtype_type **)malloc(parseIterations() * STRUCT_GTYPE_STAR_SIZE);
  for(i = 0; i < parseIterations(); i++) {
    final_indivs_data[i]=(struct gtype_type *)malloc(parseInputSamples() * STRUCT_GTYPE_SIZE);
    for(j = 0; j < parseInputSamples(); j++){
      final_indivs_data[i][j].pgtype = (ALLELE_TYPE *)malloc(parseNLoci() * sizeof(ALLELE_TYPE));
      final_indivs_data[i][j].mgtype = (ALLELE_TYPE *)malloc(parseNLoci() * sizeof(ALLELE_TYPE));
    }
  }

  double *mnals = doubleData;
  double *m = doubleData + 1 * num_samples;
  double *iis = doubleData + 2 * num_samples;
  double *lnbeta = doubleData + 3 * num_samples;
  double *hetx = doubleData + 4 * num_samples;
  double *mnehet = doubleData + 5 * num_samples;
  double *mhomo = doubleData + 6 * num_samples;
  double *varhomo = doubleData + 7 * num_samples;
  double *skhomo = doubleData + 8 * num_samples;
  double *kurhomo = doubleData + 9 * num_samples;

  const int current_sample = 0;
  const int snpN = 0;
  const int snpA = 1;
  const int snpC = 2;
  const int snpG = 3;
  const int snpT = 4;

// SNP genotype for test case

// CC TT AC NN TT AT
// CC TT AA AG GT AT
// AA GT AN AC GT AT

  storeFinalGenotype(current_sample, 0, 0, &snpC, &snpC);
  storeFinalGenotype(current_sample, 0, 1, &snpT, &snpT);
  storeFinalGenotype(current_sample, 0, 2, &snpA, &snpC);
  storeFinalGenotype(current_sample, 0, 3, &snpN, &snpN);
  storeFinalGenotype(current_sample, 0, 4, &snpT, &snpT);
  storeFinalGenotype(current_sample, 0, 5, &snpA, &snpT);
  storeFinalGenotype(current_sample, 1, 0, &snpC, &snpC);
  storeFinalGenotype(current_sample, 1, 1, &snpT, &snpT);
  storeFinalGenotype(current_sample, 1, 2, &snpA, &snpA);
  storeFinalGenotype(current_sample, 1, 3, &snpA, &snpG);
  storeFinalGenotype(current_sample, 1, 4, &snpG, &snpT);
  storeFinalGenotype(current_sample, 1, 5, &snpA, &snpT);
  storeFinalGenotype(current_sample, 2, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 2, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 2, 2, &snpA, &snpN);
  storeFinalGenotype(current_sample, 2, 3, &snpA, &snpC);
  storeFinalGenotype(current_sample, 2, 4, &snpG, &snpT);
  storeFinalGenotype(current_sample, 2, 5, &snpA, &snpT);

  // Statistic 6: mnals
  counts(numberOfAllelesPtr, final_indivs_data, mnals, gType, gcountPtr);

  // Statistic 1: m
  // Doesn't work with SNPs
  // sortM(numberOfAlleles, parseIterations(), m, gType, gcountPtr);

  // Statistic 2: iis
  // Doesn't work on this example
  twolocusiis(numberOfAlleles, parseIterations(), final_indivs_data, iis, gType, gcountPtr);

  // Statistic 3: lnbeta
  // Doesn't work with SNPs
  // beta(numberOfAlleles, parseIterations(), lnbeta, gType, gcountPtr);

  // Statistics 4 and 5: hetx, mnehet
  hetexcess(numberOfAlleles, parseIterations(), final_indivs_data, hetx, mnehet, gType, gcountPtr);

  // Statistics 7 and 8 (and 9 and 10): mhomo, varhomo, skhomo, kurhomo
  multih(parseIterations(), final_indivs_data, mhomo, varhomo, skhomo, kurhomo, gType);

  EXPECT_DOUBLE_EQ(hetx[current_sample], -0.11464968152866262);
  EXPECT_DOUBLE_EQ(mnehet[current_sample], 0.52333333333333332);
  EXPECT_DOUBLE_EQ(mnals[current_sample], 2.5);
  EXPECT_DOUBLE_EQ(mhomo[current_sample], 8. / 3.);
  EXPECT_DOUBLE_EQ(varhomo[current_sample], 7. / 3.);
  EXPECT_DOUBLE_EQ(skhomo[current_sample], -0.20782656212951636);
  EXPECT_DOUBLE_EQ(kurhomo[current_sample], -7. / 3.);

  deallocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  // Deallocate structure 5
  for(i = 0; i < parseIterations(); i++){
    for(j = 0; j < parseInputSamples(); j++){
      free(final_indivs_data[i][j].pgtype);
      free(final_indivs_data[i][j].mgtype);
    }
    free(final_indivs_data[i]);
  }
  free(final_indivs_data);

}

TEST(stats, SNPtest2){
  int i, j;
  int argc = 12;
  char a0[] = {'o', 'n', 'e', 's', 'a', 'm', 'p', '\0'};
  char a1[] = {'-', 'l', '3', '\0'};
  char a2[] = {'-', 'i', '2', '0', '\0'};
  char a3[] = {'-', 's', '\0'};
  char a4[] = {'-', 't', '1', '\0'};
  char a5[] = {'-', 'b', '8', '\0'};
  char a6[] = {'-', 'w', '\0'};
  char a7[] = {'-', 'r', 'C', '\0'};
  char a8[] = {'-', 'd', '0', '\0'};
  char a9[] = {'-', 'v', '1', '\0'};
  char a10[] = {'-', 'u', '0', '.', '5', '\0'};
  char a11[] = {'-', 'o', '0', '\0'};
  char *argv[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11};
  int initial_indivs_count = 0;
  int bottleneck_indivs_count = 0;
  int final_indivs_count = 20;
  int num_samples = 20;
  int num_loci = 3;
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;

  // Simulate parsing command line arguments
  parseArguments(argc, argv);

  // Allocate memory
  allocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

/*! \brief Allocates arrays to store sampling generations.
 */
  final_indivs_data = (struct gtype_type **)malloc(parseIterations() * STRUCT_GTYPE_STAR_SIZE);
  for(i = 0; i < parseIterations(); i++) {
    final_indivs_data[i]=(struct gtype_type *)malloc(parseInputSamples() * STRUCT_GTYPE_SIZE);
    for(j = 0; j < parseInputSamples(); j++){
      final_indivs_data[i][j].pgtype = (ALLELE_TYPE *)malloc(parseNLoci() * sizeof(ALLELE_TYPE));
      final_indivs_data[i][j].mgtype = (ALLELE_TYPE *)malloc(parseNLoci() * sizeof(ALLELE_TYPE));
    }
  }

  double *mnals = doubleData;
  double *m = doubleData + 1 * num_samples;
  double *iis = doubleData + 2 * num_samples;
  double *lnbeta = doubleData + 3 * num_samples;
  double *hetx = doubleData + 4 * num_samples;
  double *mnehet = doubleData + 5 * num_samples;
  double *mhomo = doubleData + 6 * num_samples;
  double *varhomo = doubleData + 7 * num_samples;
  double *skhomo = doubleData + 8 * num_samples;
  double *kurhomo = doubleData + 9 * num_samples;

  const int current_sample = 0;
  const int snpN = 0;
  const int snpA = 1;
  const int snpC = 2;
  const int snpG = 3;
  const int snpT = 4;

  // SNP genotype for test case
  storeFinalGenotype(current_sample, 0, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 0, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 0, 2, &snpA, &snpG);
  storeFinalGenotype(current_sample, 1, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 1, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 1, 2, &snpA, &snpT);
  storeFinalGenotype(current_sample, 2, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 2, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 2, 2, &snpG, &snpG);
  storeFinalGenotype(current_sample, 3, 0, &snpG, &snpA);
  storeFinalGenotype(current_sample, 3, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 3, 2, &snpA, &snpT);
  storeFinalGenotype(current_sample, 4, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 4, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 4, 2, &snpG, &snpG);
  storeFinalGenotype(current_sample, 5, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 5, 1, &snpG, &snpG);
  storeFinalGenotype(current_sample, 5, 2, &snpG, &snpT);
  storeFinalGenotype(current_sample, 6, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 6, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 6, 2, &snpG, &snpG);
  storeFinalGenotype(current_sample, 7, 0, &snpG, &snpG);
  storeFinalGenotype(current_sample, 7, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 7, 2, &snpC, &snpT);
  storeFinalGenotype(current_sample, 8, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 8, 1, &snpG, &snpG);
  storeFinalGenotype(current_sample, 8, 2, &snpG, &snpT);
  storeFinalGenotype(current_sample, 9, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 9, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 9, 2, &snpG, &snpC);
  storeFinalGenotype(current_sample, 10, 0, &snpA, &snpN);
  storeFinalGenotype(current_sample, 10, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 10, 2, &snpG, &snpC);
  storeFinalGenotype(current_sample, 11, 0, &snpA, &snpN);
  storeFinalGenotype(current_sample, 11, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 11, 2, &snpC, &snpC);
  storeFinalGenotype(current_sample, 12, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 12, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 12, 2, &snpC, &snpC);
  storeFinalGenotype(current_sample, 13, 0, &snpG, &snpG);
  storeFinalGenotype(current_sample, 13, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 13, 2, &snpC, &snpC);
  storeFinalGenotype(current_sample, 14, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 14, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 14, 2, &snpG, &snpC);
  storeFinalGenotype(current_sample, 15, 0, &snpA, &snpN);
  storeFinalGenotype(current_sample, 15, 1, &snpG, &snpG);
  storeFinalGenotype(current_sample, 15, 2, &snpG, &snpC);
  storeFinalGenotype(current_sample, 16, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 16, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 16, 2, &snpG, &snpC);
  storeFinalGenotype(current_sample, 17, 0, &snpN, &snpN);
  storeFinalGenotype(current_sample, 17, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 17, 2, &snpT, &snpC);
  storeFinalGenotype(current_sample, 18, 0, &snpA, &snpA);
  storeFinalGenotype(current_sample, 18, 1, &snpG, &snpG);
  storeFinalGenotype(current_sample, 18, 2, &snpG, &snpC);
  storeFinalGenotype(current_sample, 19, 0, &snpN, &snpA);
  storeFinalGenotype(current_sample, 19, 1, &snpG, &snpT);
  storeFinalGenotype(current_sample, 19, 2, &snpG, &snpC);

  // Statistic 6: mnals
  counts(numberOfAllelesPtr, final_indivs_data, mnals, gType, gcountPtr);

  // Statistic 1: m
  // Doesn't work with SNPs
  // sortM(numberOfAlleles, parseIterations(), m, gType, gcountPtr);

  // Statistic 2: iis
  twolocusiis(numberOfAlleles, parseIterations(), final_indivs_data, iis, gType, gcountPtr);

  // Statistic 3: lnbeta
  // Doesn't work with SNPs
  // beta(numberOfAlleles, parseIterations(), lnbeta, gType, gcountPtr);

  // Statistics 4 and 5: hetx, mnehet
  hetexcess(numberOfAlleles, parseIterations(), final_indivs_data, hetx, mnehet, gType, gcountPtr);

  // Statistics 7 and 8 (and 9 and 10): mhomo, varhomo, skhomo, kurhomo
  multih(parseIterations(), final_indivs_data, mhomo, varhomo, skhomo, kurhomo, gType);

  EXPECT_DOUBLE_EQ(iis[current_sample], 0.074121580056172948);
  EXPECT_DOUBLE_EQ(hetx[current_sample], -0.088891959566905987);
  EXPECT_DOUBLE_EQ(mnehet[current_sample], 0.47959048428452905);
  EXPECT_DOUBLE_EQ(mnals[current_sample], 3);
  EXPECT_DOUBLE_EQ(mhomo[current_sample], 1.25);
  EXPECT_DOUBLE_EQ(varhomo[current_sample], 0.51315789473684215);
  EXPECT_DOUBLE_EQ(skhomo[current_sample], -0.3570448637697346);
  EXPECT_DOUBLE_EQ(kurhomo[current_sample], -1.1220167652859963);

  deallocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  // Deallocate structure 5
  for(i = 0; i < parseIterations(); i++){
    for(j = 0; j < parseInputSamples(); j++){
      free(final_indivs_data[i][j].pgtype);
      free(final_indivs_data[i][j].mgtype);
    }
    free(final_indivs_data[i]);
  }
  free(final_indivs_data);
}
