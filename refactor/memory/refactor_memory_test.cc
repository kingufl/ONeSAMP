#include <gtest/gtest.h>

extern "C"{
#include "../macro/refactor_macro.h"
}

// Generally, want to make sure the data structures are large enough. Any
// invalid memory access should show up on valgrind.
TEST(memory, memory_number_of_alleles){
  int initial_indivs_count = 40;
  int bottleneck_indivs_count = 10;
  int final_indivs_count = 25;
  int num_samples = 30;
  int num_loci = 40;
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;

  int g1 = 1;
  int *g1ptr = &g1;
  int g2 = 2;
  int *g2ptr = &g2;
  int g3 = 3;
  int *g3ptr = &g3;
  int g4 = 4;
  int *g4ptr = &g4;

  allocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  numberOfAlleles[29][0] = 4;

  deallocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);
}

TEST(memory, memory_final_genotype_type_count){
  int initial_indivs_count = 40;
  int bottleneck_indivs_count = 10;
  int final_indivs_count = 25;
  int num_samples = 30;
  int num_loci = 40;
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;

  int g1 = 1;
  int *g1ptr = &g1;
  int g2 = 2;
  int *g2ptr = &g2;
  int g3 = 3;
  int *g3ptr = &g3;
  int g4 = 4;
  int *g4ptr = &g4;

  allocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  gType[29][0][0] = 1;
  gType[29][num_loci - 1][MAX_NO_ALLELES - 1] = 2;

  deallocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);
}

TEST(memory, memory_final_genotype_frequency_count){
  int initial_indivs_count = 40;
  int bottleneck_indivs_count = 10;
  int final_indivs_count = 25;
  int num_samples = 30;
  int num_loci = 40;
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;

  int g1 = 1;
  int *g1ptr = &g1;
  int g2 = 2;
  int *g2ptr = &g2;
  int g3 = 3;
  int *g3ptr = &g3;
  int g4 = 4;
  int *g4ptr = &g4;

  allocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  gcount[0][0][0] = 1;
  gcount[final_indivs_count - 1][num_loci - 1][MAX_NO_ALLELES - 1] = 2;

  deallocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);
}

TEST(memory, memory_initial_genotype_data){
  int initial_indivs_count = 40;
  int bottleneck_indivs_count = 10;
  int final_indivs_count = 25;
  int num_samples = 30;
  int num_loci = 40;
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;

  int g1 = 1;
  int *g1ptr = &g1;
  int g2 = 2;
  int *g2ptr = &g2;
  int g3 = 3;
  int *g3ptr = &g3;
  int g4 = 4;
  int *g4ptr = &g4;

  allocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  storeInitialGenotype(30, 3, g1ptr, g2ptr);
  loadInitialGenotype(30, 3, g3ptr, g4ptr);
  ASSERT_EQ(g1, g3);
  ASSERT_EQ(g2, g4);

  storeInitialGenotype(35, 4, g1ptr, g2ptr);
  loadInitialGenotype(35, 4, g3ptr, g4ptr);
  ASSERT_EQ(g1, g3);
  ASSERT_EQ(g2, g4);

  deallocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);
}

TEST(memory, memory_double_array_access){
  int initial_indivs_count = 40;
  int bottleneck_indivs_count = 10;
  int final_indivs_count = 25;
  int num_samples = 30;
  int num_loci = 40;
  int **numberOfAlleles;
  int ***numberOfAllelesPtr = &numberOfAlleles;
  double *doubleData;
  double **doubleDataPtr = &doubleData;
  int ***gType;
  int ****gTypePtr = &gType;
  int ***gcount;
  int ****gcountPtr = &gcount;

  int g1 = 1;
  int *g1ptr = &g1;
  int g2 = 2;
  int *g2ptr = &g2;
  int g3 = 3;
  int *g3ptr = &g3;
  int g4 = 4;
  int *g4ptr = &g4;

  allocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  double *kurhomo = doubleData + 9 * num_samples;
  kurhomo[num_samples - 1] = 0.5;

  deallocateOneSampMemory(initial_indivs_count, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);
}


