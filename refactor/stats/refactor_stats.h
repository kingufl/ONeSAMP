#include "../macro/refactor_macro.h"

#ifndef REFACTOR_STATS_H
#define REFACTOR_STATS_H

extern int ***gcount;

double fallingQuotient(double s, double t1, double t2, int c);
double allelePr(int val1, int val2, double theta);
int minAlleleCount();

void mutate(ALLELE_TYPE *gene, int motif);
void mutateSNP(ALLELE_TYPE *gene);
void mutateMicroSat(ALLELE_TYPE *gene, int motif);

void assort(int nextgen, gtype_type *offvec, gtype_type *mothers, gtype_type *fathers, int indivs, int samp, int num_loci);
void counts(int ***numberOfAllelesPtr, gtype_type **samp_data, double mnals[], int ***gType, int ****gcountPtr);
void sortM(int **numberOfAlleles, int num_samples, double m[], int ***gType, int ****gcountPtr);
void twolocusiis(int **numberOfAlleles,int num_samples, gtype_type **samp_data, double iis[], int ***gType, int ****gcountPtr);  // Weir composite LD estimator using all alleles
void beta(int **numberOfAlleles, int num_samples, double lnbeta[], int ***gType, int ****gcountPtr);
void hetexcess(int **numberOfAlleles,int num_samples, gtype_type **samp_data, double hetx[], double mnehet[], int ***gType, int ****gcountPtr); // Need skiploc operational for missing data
void multih(int num_samples, gtype_type **samp_data, double mhomo[], double varhomo[], double skhomo[], double kurhomo[], int ***gType);

#endif
