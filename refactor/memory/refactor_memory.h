#include "../macro/refactor_macro.h"

#ifndef REFACTOR_MEMORY_H
#define REFACTOR_MEMORY_H
void allocateStructure1(int initial_indivs_count, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr);
void allocateOneSampMemory(int initial_indivs_count, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr);
void deallocateOneSampMemory(int initial_indivs_count, int bottleneck_indivs_count, int final_indivs_count, int num_samples, int num_loci, int ***numberOfAllelesPtr, double **doubleDataPtr, int ****gTypePtr, int ****gcountPtr);
void storeInitialGenotype(int individual, int index, int *genotype1, int *genotype2);
void loadInitialGenotype(int individual, int index, int *genotype1, int *genotype2);
void storeFinalGenotype(int sample, int individual, int index, const int *genotype1, const int *genotype2);
void loadFinalGenotype(int sample, int individual, int index, int *genotype1, int *genotype2);
#endif
