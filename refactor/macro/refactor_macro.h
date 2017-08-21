/*! \file refactor_macro.h
 *  \brief Macro commands for the OneSamp program.
 */

#ifndef REFACTOR_MACRO_H
#define REFACTOR_MACRO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>

/* \def MAX_NO_ALLELES
 * \brief Defines maximum number of alleles.
 *
 * This value must float based on the number of input loci.
 * Because all of the input alleles could differ, we must be careful.
 * If we have SNPs, there are at most 5 options.
 * Otherwise, there could be as many options as there are alleles in the locus.
 */
#define MAX_NO_ALLELES (parseFormFlag() == 0 ? 5 : 999)

/* \def ALLELE_TYPE
 * \brief Type of variable in which to store allele.
 *
 */
#define ALLELE_TYPE short

#ifndef TYPEDEFS
#define TYPEDEFS
struct gtype_type {
    ALLELE_TYPE *pgtype;
    ALLELE_TYPE *mgtype;
};
typedef struct gtype_type gtype_type;
#endif

/* \def DOUBLE_TOLERANCE
 * \brief Tolerance for doubles in unit tests.
 */ 
#define DOUBLE_TOLERANCE (0.00001)

// Standard true and false macros

/* \def TRUE
 * \brief Expresses the boolean value TRUE.
 */
#define TRUE (0 == 0)

/* \def FALSE
 * \brief Expresses the boolean value FALSE.
 */
#define FALSE (0 == 1)

/*! \def C_RANDOM_FLAG
 *  \brief If true, use C random numbers. If false, use GFSR random numbers.
 */
#define C_RANDOM_FLAG (parseRFlag())

/* \def ASCII_ZERO
 * \brief Value of zero in ASCII: used for character conversion in parsing.
 */
#define ASCII_ZERO (48)

/* \def extraProportionOfBufferLoci
 * \brief Create extra loci so that we can throw out monomorphic ones when they arise.
 * Larger values are more "safe", but will use more memory and will cause
 * the program to run more slowly.
 */
#define extraProportionOfBufferLoci 2

// GLOBAL VARIABLES
// TODO Try to get rid of these extern statements.
extern struct gtype_type *initial_indivs_data, **final_indivs_data;

// Only one of the four options may be uncommented
//
// OPTION 1 (the rand_table array is of type int[])
//
#define GFSR_STYPE int
#define GFSR_STYPE_UNSIGNED unsigned int
#define GFSR_STYPE_UNSIGNED_MAX UINT_MAX
//
// OPTION 2 (the rand_table array is of type unsigned int[])
//
//#define GFSR_STYPE unsigned int
//#define GFSR_STYPE_UNSIGNED unsigned int
//#define GFSR_STYPE_UNSIGNED_MAX UINT_MAX
//
// OPTION 3 (the rand_table array is of type long[])
//
//#define GFSR_STYPE long
//#define GFSR_STYPE_UNSIGNED unsigned long
//#define GFSR_STYPE_UNSIGNED_MAX ULONG_MAX
//
// OPTION 4 (the rand_table array is of type unsigned long[])
//
//#define GFSR_STYPE unsigned long
//#define GFSR_STYPE_UNSIGNED unsigned long
//#define GFSR_STYPE_UNSIGNED_MAX ULONG_MAX

/*! \def int STRUCT_GTYPE_SIZE 
 *  \brief Size of genotype structure.
 */
#define STRUCT_GTYPE_SIZE (2 * sizeof(int *))

/*! \def int STRUCT_GTYPE_STAR_SIZE
 *  \brief Size of pointer to genotype structure.
 */
#define STRUCT_GTYPE_STAR_SIZE (sizeof(int **))

// RANDOMIZATION

// This contains a group of pseudorandom number macros for OneSamp. The
// source of these algorithms originally came from:
//
// T. G. Lewis , W. H. Payne, Generalized Feedback Shift Register
//   Pseudorandom Number Algorithm, Journal of the ACM (JACM), v.20 n.3,
//   p.456-468, July 1973  [doi>10.1145/321765.321777] 
//
// A partial implementation includes code from the R statistical package
// circa 1995, and further algorithms from the paper were added to this
// codebase.

/*! \def GFSR_SBITS()
 *  \brief The number of bits contained in the GFSR register
 */
#define GFSR_SBITS() (8 * sizeof(GFSR_STYPE))

/*! \def int GFSR_RESET_PERIOD()
 *  \brief Length of time between calls to GFSR bits.
 */
#define GFSR_RESET_PERIOD() (100 * (P()))

/*! \def int GFSR_FLUSH_ITERATIONS()
 *  \brief Number of iterations used to flush GFSR register.
 */
#define GFSR_FLUSH_ITERATIONS() (5000 * (P()))

/*! \def const char *GFSR_INIT_NAME()
 *  \brief Name of variable used to represent file name of init variable.
 */
#define GFSR_INIT_NAME() (filenameGFSR)

/*! \def P()
 *  \brief Leading exponent of irreducible polynomial for pseudorandom GFSR.
 *
 *  The irreducible polynomial used is x^p + x^q + 1.
 */
#define P() (98)

/*! \def Q()
 *  \brief Second exponent of irreducible polynomial for pseudorandom GFSR.
 *
 *  The irreducible polynomial used is x^p + x^q + 1.
 */
#define Q() (27)

/*! \def int branchWithProbability(float p)
 *  \brief Sets bit with input probability based on pseudorandom number.
 */
#define branchWithProbability(p) (p > gfsr4() ? 1 : 0)

/*! \def int uniformRandomDistribution(l, t)
 *  \brief Returns a random or pseudorandom between l and t.
 */
#define uniformRandomDistribution(l, t) disrand(l, t)

/*! \def int randBit()
 *  \brief Returns a random or pseudorandom bit.
 */
#define randBit() uniformRandomDistribution(0, 1)

/*! \def randomQuantizedInterval(float min, float max, float step)
 *  \brief Returns a value occurring at a discrete step in an interval
 */
#define randomQuantizedIntervalSelection(min, max, step) \
  (step * disrand((int)(min/step),(int)(max/step)))

/*! \def void resetgfsr()
 *  \brief Resets the pseudorandom GFSR to standard values.
 *
 * Initializes the GFSR as described in (Lewis and Payne 1973). This is a C
 * equivalent of the FORTRAN code given in that paper.
 * After each call, the intrand() calls are reset.
 *
 * If the GFSR is initialized with the macros above as follows:
 *
 * P = 98
 * Q = 27
 * reset period = 100 * P
 * flush iterations = 5000 * P
 *
 * we can reproduce the values in (Lewis and Payne 1973 Fig. 6), and these
 * are given here (to within working precision).
 *
 * After the appropriate initialization macros are called, these
 * statements return the floats from the paper:
 *
 * float f1 = gfsr4(); // 0.369633
 * float f2 = gfsr4(); // 0.406314
 * float f3 = gfsr4(); // 0.428778
 * float f4 = gfsr4(); // 0.474114
 * float f5 = gfsr4(); // 0.953158
 * printf("%f\n%f\n%f\n%f\n%f\n", f1, f2, f3, f4, f5);
 *
 * Per Lewis and Payne, the general idea of the initialization operation is
 * three loops:
 *
 * 1. Initialize GFSR memory to zero.
 *
 * 2. As we range over bits of working precision i: set the ith bit of all
 * words to 1, extract GFSR_RESET_PERIOD() random values and discard them.
 *
 * 3. To flush initial values out of the register that do not have good
 * statistical properties, extract GFSR_FLUSH_ITERATIONS() more intrand values.
 */ 
#define resetgfsr() { \
  assert((GFSR_SBITS()) >= 1); \
  assert((GFSR_RESET_PERIOD()) >= 0); \
  assert((GFSR_FLUSH_ITERATIONS()) >= 0); \
  int temp_macro_resetgfsr_counter; \
  int temp_macro_resetgfsr_maskindex; \
\
\
  for(temp_macro_resetgfsr_counter = 0; \
      temp_macro_resetgfsr_counter < P(); \
      temp_macro_resetgfsr_counter++) \
      rand_table[temp_macro_resetgfsr_counter] &= 0; \
\
\
  for(temp_macro_resetgfsr_maskindex = 0; \
      temp_macro_resetgfsr_maskindex < (GFSR_SBITS()); \
      temp_macro_resetgfsr_maskindex++) { \
    for(temp_macro_resetgfsr_counter = 0; \
        temp_macro_resetgfsr_counter < P(); \
        temp_macro_resetgfsr_counter++) \
        rand_table[temp_macro_resetgfsr_counter] \
        |= 1 << temp_macro_resetgfsr_maskindex; \
    for(temp_macro_resetgfsr_counter = 0; \
        temp_macro_resetgfsr_counter < (GFSR_RESET_PERIOD()); \
        temp_macro_resetgfsr_counter++) intrand(); \
  } \
\
\
  for(temp_macro_resetgfsr_counter = 0; \
      temp_macro_resetgfsr_counter < (GFSR_FLUSH_ITERATIONS()); \
      temp_macro_resetgfsr_counter++) intrand(); \
}

// GLOBALS
// See refactor_macro.c for definitions.

extern GFSR_STYPE rand_table[P()];
extern int jindic;
extern const char filenameGFSR[];

// I copied explicitly the functionality of the original code, modifying some
// aspects of how it is interpereted (which I believe is okay because it is
// GPL), but in some cases modified the formatting or other internal aspects
// may or may not have been modified. The following listing of the
// pseudorandom number generater has the following license:

/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

// BEGIN GPL v2 R code insertion

/*! \def int intrand()
 *  \brief Generates random int if C_RANDOM_FLAG is set, else GFSR random int.
 */
#define intrand() \
  (C_RANDOM_FLAG ? rand() : \
  (++jindic, jindic %= (P()), \
  rand_table[jindic] \
    = rand_table[jindic]^rand_table[(jindic + (Q())) % (P())]))

/*! \def int disrand(l, t)
 *  \brief Generates random int between l and t. NOTE: must buffer inputs.
 */
#define disrand(l, t) \
  ((int) ((((unsigned int) intrand())%((t) - (l) + 1)) + (l)))

/*! \def float gfsr4()
 *  \brief Returns a random float in the range [0,1).
 */
#define gfsr4() (C_RANDOM_FLAG ? ((float) rand()) / INT_MAX : (float)(( GFSR_STYPE_UNSIGNED ) intrand()) / \
                   (GFSR_STYPE_UNSIGNED_MAX + 1.0))

/*! \def void opengfsr()
 *  \brief Opens the random number table.
 *
 * More specifically, if the C_RANDOM_FLAG is set, then standard C library
 * routines will be used to generate random numbers.
 *
 * If not, then if the GFSR filename exists, the old state of the GFSR
 * will be used to generate pseudorandom numbers.
 *
 * In the last case, the standard algorithm for resetting a GFSR as used by
 * Lewis and Payne will be used to reset the pseudorandom number generator.
 * 
 * EDIT: this algorithm no longer reads jindic from the file. By carefully
 * setting up the algorithm to write to the file, we can assume jindic is
 * always zero and thereby simplify the input file.
 */
#define opengfsr() \
{ \
  if(C_RANDOM_FLAG) \
  { \
    srand(time(NULL)); \
  } \
  else \
  { \
    int temp_macro_opengfsr_counter; \
    if(strlen(GFSR_INIT_NAME()) != 0) \
    {\
      FILE *temp_macro_opengfsr_file = fopen(GFSR_INIT_NAME(),"r"); \
      if(temp_macro_opengfsr_file==NULL) \
        { printf("I need %s! Where is %s?\n", \
          GFSR_INIT_NAME(), GFSR_INIT_NAME()); exit(1); } \
      for(temp_macro_opengfsr_counter = 0; \
          temp_macro_opengfsr_counter < P(); \
          temp_macro_opengfsr_counter++) \
        fscanf(temp_macro_opengfsr_file,"%d", \
          &rand_table[temp_macro_opengfsr_counter]);\
      fclose(temp_macro_opengfsr_file); \
      jindic = 0; \
    } \
    else resetgfsr(); \
  } \
}

/*! \def void closegfsr()
 *  \brief Closes the GFSR random number table.
 *
 * EDIT: this algorithm no longer writes jindic to the file. We will write
 * the output of the file so that we can assume jindic is zero.
 */
#define closegfsr() \
{ \
  if(!C_RANDOM_FLAG) \
  { \
  int temp_macro_closegfsr_j; \
  FILE *temp_macro_closegfsr_rt = fopen(GFSR_INIT_NAME(),"w"); \
  for(temp_macro_closegfsr_j = jindic; \
      temp_macro_closegfsr_j < P() + jindic; \
      temp_macro_closegfsr_j++) \
  fprintf(temp_macro_closegfsr_rt,"%d\n", \
    rand_table[temp_macro_closegfsr_j % P()]); \
  fclose(temp_macro_closegfsr_rt); \
  } \
}

// END R code insertion

#include "../arguments/refactor_arguments.h"
#include "../engine/refactor_engine.h"
#include "../parser/refactor_parser.h"
#include "../memory/refactor_memory.h"
#include "../stats/refactor_stats.h"

void writeoutput(gtype_type **samp_data, int final_indivs_count);
void numberOfAllelesDump(int **numberOfAlleles);
void gTypeDump(int ***gType, int **numberOfAlleles);
void gcountDump(int ***gCount, int **numberOfAlleles);

#endif
