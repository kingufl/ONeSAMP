#include <gtest/gtest.h>

extern "C"{
#include "refactor_arguments.h"
}

// onesamp -l10 -i11 -s -t3 -b8 -f13 -rC -u0.5 -v1
TEST(arguments, testA){
  int argc = 12;
  char a0[] = {'o', 'n', 'e', 's', 'a', 'm', 'p', '\0'};
  char a1[] = {'-', 'l', '1', '0', '\0'};
  char a2[] = {'-', 'i', '1', '1', '\0'};
  char a3[] = {'-', 's', '\0'};
  char a4[] = {'-', 't', '3', '\0'};
  char a5[] = {'-', 'b', '8', '\0'};
  char a6[] = {'-', 'u', '0', '.', '5', '\0'};
  char a7[] = {'-', 'r', 'C', '\0'};
  char a8[] = {'-', 'd', '4', '\0'};
  char a9[] = {'-', 'v', '1', '\0'};
  char a10[] = {'-', 'o', '1', '\0'};
  char a11[] = {'-', 'g', '\0'};
  char *argv[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11};
  parseArguments(argc, argv);
  EXPECT_EQ(parseProgramName()[0], 'o');  // Verify program name
  EXPECT_EQ(parseProgramName()[1], 'n');
  EXPECT_EQ(parseProgramName()[2], 'e');
  EXPECT_EQ(parseProgramName()[3], 's');
  EXPECT_EQ(parseProgramName()[4], 'a');
  EXPECT_EQ(parseProgramName()[5], 'm');
  EXPECT_EQ(parseProgramName()[6], 'p');
  EXPECT_EQ(parseProgramName()[7], '\0');
  EXPECT_EQ(parseFormFlag(), FALSE);      // Verify interpreted SNPs
  EXPECT_EQ(parseIterations(), 1);        // Verify number of iterations
  EXPECT_EQ(parseNLoci(), 10);            // Verify number of loci
  EXPECT_EQ(parseInputSamples(), 11);     // Verify input size
  EXPECT_EQ(parseBottleneckMin(), 4);        // Verify bottleneck size
  EXPECT_EQ(parseBottleneckLengthMin(), 4);  // Verify bottleneck length
  EXPECT_EQ(parseThetaMin(), 1);             // Verify theta parameter
  EXPECT_EQ(parseMRateMin(), 0.5);           // Verify mutation rate
  EXPECT_EQ(parseSyntaxCheck(), FALSE);   // Not doing a syntax check
  flushArguments();
}

TEST(arguments, testB){
  int argc = 11;
  char a0[] = {'o', 'n', 'e', 's', 'a', 'm', 'p', '\0'};
  char a1[] = {'-', 'i', '7', '1', '\0'};
  char a2[] = {'-', 'u', '0', '.', '1', '2', '\0'};
  char a3[] = {'-', 't', '5', '6', '\0'};
  char a4[] = {'-', 'd', '7', '\0'};
  char a5[] = {'-', 'v', '0', '.', '1', '\0'};
  char a6[] = {'-', 'b', '1', '0', '\0'};
  char a7[] = {'-', 'r', 'R', 'E', 'S', 'E', 'T', '\0'};
  char a8[] = {'-', 's', '\0'};
  char a9[] = {'-', 'l', '4', '0', '\0'};
  char a10[] = {'-', 'x', '\0'};
  char *argv[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10};
  parseArguments(argc, argv);
  EXPECT_EQ(randBit(), 1);                 // If the GFSR is reset correctly,
  EXPECT_EQ(randBit(), 1);                 // the random bits will be
  EXPECT_EQ(randBit(), 0);                 // standard.
  EXPECT_EQ(randBit(), 0);
  EXPECT_EQ(randBit(), 0);
  EXPECT_EQ(randBit(), 0);
  EXPECT_EQ(randBit(), 1);
  EXPECT_EQ(randBit(), 0);
  EXPECT_EQ(randBit(), 0);
  EXPECT_EQ(parseProgramName()[0], 'o');   // Verify program name
  EXPECT_EQ(parseProgramName()[1], 'n');
  EXPECT_EQ(parseProgramName()[2], 'e');
  EXPECT_EQ(parseProgramName()[3], 's');
  EXPECT_EQ(parseProgramName()[4], 'a');
  EXPECT_EQ(parseProgramName()[5], 'm');
  EXPECT_EQ(parseProgramName()[6], 'p');
  EXPECT_EQ(parseProgramName()[7], '\0');
  EXPECT_EQ(parseFormFlag(), FALSE);      // Verify interpreted microsats
  EXPECT_EQ(parseIterations(), 56);       // Verify number of iterations
  EXPECT_EQ(parseNLoci(), 40);            // Verify number of loci
  EXPECT_EQ(parseInputSamples(), 71);     // Verify input size
  EXPECT_EQ(parseBottleneckMin(), 5);       // Verify bottleneck size
  EXPECT_EQ(parseBottleneckLengthMin(), 7);  // Verify bottleneck length
  EXPECT_EQ(parseThetaMin(), 0.1);           // Verify theta parameter
  EXPECT_EQ(parseMRateMin(), 0.12);          // Verify mutation rate
  EXPECT_EQ(parseSyntaxCheck(), TRUE);    // Verify doing a syntax check
  EXPECT_EQ(parseExample(), FALSE);       // Verify not doing an example
  flushArguments();
  closegfsr();
}

// onesamp -i71 -u0.12 -72 -b10 -v0.1 -f13 -rGFSR -m -l19
TEST(arguments, testC){
  int argc = 12;
  char a0[] = {'o', 'n', 'e', 's', 'a', 'm', 'p', '\0'};
  char a1[] = {'-', 'd', '1', '0', '\0'};
  char a2[] = {'-', 'u', '0', '.', '1', '2', '\0'};
  char a3[] = {'-', 't', '7', '2', '\0'};
  char a4[] = {'-', 'b', '1', '0', '\0'};
  char a5[] = {'-', 'v', '0', '.', '1', '\0'};
  char a6[] = {'-', 'i', '7', '1', '\0'};
  char a7[] = {'-', 'r', 'G', 'F', 'S', 'R', '\0'};
  char a8[] = {'-', 'm', '\0'};
  char a9[] = {'-', 'l', '1', '9', '\0'};
  char a10[] = {'-', 'p', '\0'};
  char a11[] = {'-', 'o', '1', '\0'};
  char *argv[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11};
  parseArguments(argc, argv);
  EXPECT_EQ(parseProgramName()[0], 'o');  // Verify program name
  EXPECT_EQ(parseProgramName()[1], 'n');
  EXPECT_EQ(parseProgramName()[2], 'e');
  EXPECT_EQ(parseProgramName()[3], 's');
  EXPECT_EQ(parseProgramName()[4], 'a');
  EXPECT_EQ(parseProgramName()[5], 'm');
  EXPECT_EQ(parseProgramName()[6], 'p');
  EXPECT_EQ(parseProgramName()[7], '\0');
  EXPECT_EQ(parseFormFlag(), TRUE);       // Verify interpreted microsats
  EXPECT_EQ(parseIterations(), 72);       // Verify number of iterations
  EXPECT_EQ(parseNLoci(), 19);            // Verify number of loci
  EXPECT_EQ(parseInputSamples(), 71);     // Verify input size
  EXPECT_EQ(parseBottleneckMin(), 5);       // Verify bottleneck size
  EXPECT_EQ(parseThetaMin(), 0.1);           // Verify theta parameter
  EXPECT_EQ(parseMRateMin(), 0.12);          // Verify mutation rate
  EXPECT_EQ(parseSyntaxCheck(), FALSE);   // Not doing a syntax check
  EXPECT_EQ(parseBottleneckLengthMin(), 10); // Verify bottleneck length
  EXPECT_EQ(parseExample(), FALSE);        // Verify not simulating an example
  EXPECT_EQ(parseExamplePop(), TRUE);
  flushArguments();
}
