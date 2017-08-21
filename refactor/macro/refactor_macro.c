/*! \file refactor_macro.c
 *  \brief Functions and variables for utility macros in OneSamp.
 */

#include "refactor_macro.h"

/*! \var GFSR_STYPE rand_table[]
 *  \brief Random number table structure that stores state of random numbers.
 */
GFSR_STYPE rand_table[P()];

/*! \var int jindic
 *  \brief Current index of number in GFSR table.
 */
int jindic;

/*! \var const char filename[]
 *  \brief File name of GFSR table.
 */
const char filenameGFSR[] = "INITFILE";

/*! \def writeoutput(struct gtype_type **samp_data)
 *  \brief displays genotype data from final generation
 */
void writeoutput(struct gtype_type **samp_data, int final_indivs_count)
{
  int num_loci = parseNLoci();
  int num_samples = parseIterations();

  int j,k,i;
  printf("Auto-generated genotype output.\n");
    for (j = 0; j < num_loci;++j) {
      printf("%1d\n",j+1);
  }
  printf("Pop\n");
  for(k = 0; k < num_samples; k++) {
    for(j = 0; j < final_indivs_count; j++){
      printf("%d , ", j+1);
      for(i = 0; i < num_loci;++i){
        printf("%02d%02d ", samp_data[k][j].mgtype[i], samp_data[k][j].pgtype[i]);
      }
      printf("\n");
    }
  }
}

/*! \def numberOfAllelesDump(int **numberOfAlleles)
 *  \brief displays count of kinds of alleles at each locus in final generation
 */
void numberOfAllelesDump(int **numberOfAlleles){
  int num_loci = parseNLoci();
  int num_samples = parseIterations();
  int k, i;
  printf("Auto-generated number of alleles output.\n");
  for(k = 0; k < num_samples; k++) {
    for(i = 0; i < num_loci;++i){
      printf("%d ", numberOfAlleles[k][i]);
    }
    printf("\n");
  }
}

/*! \def TypeDump(int ***gType, int **numberOfAlleles)
 *  \brief displays genotypes of final generation
 */
void gTypeDump(int ***gType, int **numberOfAlleles){
  int num_loci = parseNLoci();
  int num_samples = parseIterations();
  int k, i, j;
  printf("Auto-generated gType output.\n");
  for(k = 0; k < num_samples; k++) {
    for(i = 0; i < num_loci;++i){
      printf("(");
      for(j = 0; j < numberOfAlleles[k][i]; j++){
        printf("%2d ", gType[k][i][j]);
      }
      printf(")");
    }
    printf("\n");
  } 
}

/*! \def gcountDump(int ***gCount, int **numberOfAlleles)
 *  \brief displays counts of alleles in final generation
 */
void gcountDump(int ***gCount, int **numberOfAlleles){
  int num_loci = parseNLoci();
  int num_samples = parseIterations();
  int k, i, j;
  printf("Auto-generated gcount output.\n");
  for(k = 0; k < num_samples; k++) {
    for(i = 0; i < num_loci; i++){
      printf("(");
      for(j = 0; j < numberOfAlleles[k][i]; j++){
        printf("%2d ", gCount[k][i][j]);
      }
      printf(")");
    }
    printf("\n");
  } 
}
