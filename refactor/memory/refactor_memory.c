#include "../macro/refactor_macro.h"

// TODO: struct gtype_type will have static size when allocated with malloc
#ifndef ONESAMP_GLOBALS
#define ONESAMP_GLOBALS

struct gtype_type *initial_indivs_data, **final_indivs_data;
#endif

// Allocate structure 1
// ALLOCATE DYNAMIC MEMORY TO TRACK INFO EACH ITERATION

/*! \brief Allocates arrays to store number of alleles.
 */
void allocateStruct1(int initial_inidivs_count_allocation, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci_allocation, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr){
  int i;
  int j;
  *numberOfAllelesPtr = (int **)malloc(num_samples * sizeof(int *));
  for(j = 0; j < num_samples; j++) (*numberOfAllelesPtr)[j] = (int *)malloc(num_loci_allocation * sizeof(int));
}

/*! \brief Allocates arrays to track number of alleles at each locus.
 */
void allocateStruct2(int initial_inidivs_count_allocation, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci_allocation, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr){
  int i;
  int j;
  *gTypePtr = (int ***)malloc(num_samples * sizeof(int **));
  for(i = 0; i < num_samples; i++){
    (*gTypePtr)[i] = (int **)malloc(num_loci_allocation*sizeof(int *));
    for(j = 0; j < num_loci_allocation; j++)(*gTypePtr)[i][j] = (int *)malloc(MAX_NO_ALLELES*sizeof(int));
  }
}

/*! \brief Allocates arrays to track frequency of alleles at each locus.
 */
void allocateStruct3(int initial_inidivs_count_allocation, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci_allocation, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr){
  int i;
  int j;
  *gcountPtr = (int ***)malloc(num_samples * sizeof(int **));
  for(i = 0; i < num_samples; i++){
    (*gcountPtr)[i] = (int **)malloc(num_loci_allocation*sizeof(int *));
    for(j = 0; j < num_loci_allocation; j++)(*gcountPtr)[i][j] = (int *)malloc(MAX_NO_ALLELES*sizeof(int));
  }
}

/*! \brief Allocates arrays to store initial data.
 */
void allocateStruct4(int initial_indivs_count_allocation, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci_allocation, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr){
  int i;
  int j;
  initial_indivs_data = (struct gtype_type *)malloc(initial_indivs_count_allocation * STRUCT_GTYPE_SIZE);
  for(j = 0; j < initial_indivs_count_allocation; j++){
    initial_indivs_data[j].pgtype = (ALLELE_TYPE *)malloc(num_loci_allocation*sizeof(ALLELE_TYPE));
    initial_indivs_data[j].mgtype = (ALLELE_TYPE *)malloc(num_loci_allocation*sizeof(ALLELE_TYPE));
  }
}

/*! \brief Allocates arrays to store generation of statistics.
 */
void allocateStruct6(int initial_inidivs_count_allocation, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci_allocation, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr){
  int i;
  int j;
  *doubleDataPtr = (double *)malloc(11 * num_samples * sizeof(double));
}

/*! \brief Invokes calls to allocate all global data structures.
 */
void allocateOneSampMemory(int initial_inidivs_count_allocation, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci_allocation, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr){
  allocateStruct1(initial_inidivs_count_allocation, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci_allocation, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  allocateStruct2(initial_inidivs_count_allocation, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci_allocation, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  allocateStruct3(initial_inidivs_count_allocation, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci_allocation, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  allocateStruct4(initial_inidivs_count_allocation, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci_allocation, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  allocateStruct6(initial_inidivs_count_allocation, bottleneck_indivs_count, final_indivs_count, num_samples, num_loci_allocation, numberOfAllelesPtr, doubleDataPtr, gTypePtr, gcountPtr);

  // Remaining random arrays handled in arguments file.
}

/*! \brief Invokes calls to deallocate all global data structures.
 */
void deallocateOneSampMemory(int initial_inidivs_count_allocation, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci_allocation, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr){
  int i;
  int j;

  // Deallocate structure 6
  free(*doubleDataPtr);

  // Deallocate structure 4
  for(j = 0; j < initial_inidivs_count_allocation; j++){
    free(initial_indivs_data[j].pgtype);
    free(initial_indivs_data[j].mgtype);
  }
  free(initial_indivs_data);

  // Deallocate structure 3
  for(i = 0; i < num_samples; i++) {
    for(j = 0; j < num_loci_allocation; j++) free((*gcountPtr)[i][j]);
    free((*gcountPtr)[i]);
  }
  free(*gcountPtr);

  // Deallocate structure 3
  for(i = 0; i < num_samples; i++) {
    for(j = 0; j < num_loci_allocation; j++) free((*gTypePtr)[i][j]);
    free((*gTypePtr)[i]);
  }
  free(*gTypePtr);

  // Deallocate structure 2
  for(i = 0; i < num_samples; i++) free((*numberOfAllelesPtr)[i]);
  free(*numberOfAllelesPtr);
}

/*! \brief Loads an initial genotype from memory.
 */
void loadInitialGenotype(int individual, int index, int *genotype1, int *genotype2){
  *genotype1 = initial_indivs_data[individual].pgtype[index];
  *genotype2 = initial_indivs_data[individual].mgtype[index];
}

/*! \brief Stores an initial genotype to memory.
 */
void storeInitialGenotype(int individual, int index, int *genotype1, int *genotype2){
  initial_indivs_data[individual].pgtype[index] = *genotype1;
  initial_indivs_data[individual].mgtype[index] = *genotype2;
}

/*! \brief Loads a genotype from the final sample.
 */
void loadFinalGenotype(int sample, int individual, int index, int *genotype1, int *genotype2){
  *genotype1 = final_indivs_data[sample][individual].pgtype[index];
  *genotype2 = final_indivs_data[sample][individual].mgtype[index];
}

/*! \brief Stores a genotype to st.
 */
void storeFinalGenotype(int sample, int individual, int index, const int *genotype1, const int *genotype2){
  final_indivs_data[sample][individual].pgtype[index] = *genotype1;
  final_indivs_data[sample][individual].mgtype[index] = *genotype2;
}
