#include <gtest/gtest.h>

extern "C" {
#include "../macro/refactor_macro.h"
}


// This is the pseudorandom number configuration state as determined by
// the source (Lewis and Payne 1973) mentioned elsewhere.
const int resetgfsr_standard[] =
{-1787141026,   692513452,  1183199546,  -408704453, -1948498849,
  -741268548, -2056134124,   344472089,  1971513546, -1186404342,
 -1288692714,  1300795239, -1058175985,  1279878134, -1398448201,
  2093706256,  1318340072,  2069868445,   559626742,   653860201,
   734005281,  1296960365,  -475499606,  1237127688,  1691062534,
   584524938,   598826735,   -15324583,  2011606675,   780278841,
 -1972910449,  -226245823,   668140974,  1131266630,   249592506,
   -86628624,  -275463608,  -295592482,   166446056, -1112309364,
  1293403676,  -424242630,  1590026273,  1360713837,  -751543611,
 -1645096292,  -555286679,   713491269,  2123841324,  1229902981,
   523752922,  1407975601,  1594927896,   356479372, -1011550732,
 -1215576184, -1625113692,  1881095498,  -379674164,  -537984567,
  -227157412,  -872755286,   -18090145,  1295468477, -1182986326,
 -1874750317,  -708581624,  -655308140,  1503686128,   691242801,
  1150891949, -1013129773, -1536584373, -1059302984,  1085932207,
 -1683257391, -1341523941, -1361343899,  -610447094,  -961687630,
   434014805,  1370456709,  1804175579,    64864485,  1579425988,
  1405582889, -2132743785, -1150733498,   643025249,  1288826228,
  1727979288, -1597603817,   168759895,  -383327804,  1882949212,
  1969381119,  -706548769, -885815584};

// Tests all the macro functions.
// Got this working correctly.
TEST(macro, macro_true){
  EXPECT_TRUE(TRUE);
}

TEST(macro, macro_false){
  EXPECT_TRUE(!FALSE);
}

TEST(macro, macro_ascii_zero){
  EXPECT_TRUE(ASCII_ZERO == 48);
}

TEST(macro, struct_gtype_size){
  EXPECT_TRUE(STRUCT_GTYPE_SIZE == (2 * sizeof(int *)));
}

TEST(macro, macro_p){
  EXPECT_TRUE(P() == 98);
}

TEST(macro, macro_q){
  EXPECT_TRUE(Q() == 27);
}

// These next EXPECT_EQ statements enforce a particular sequence
// of 2^98 - 1 (10 ^ 29) pseudorandom numbers.
// See (Lewis and Payne 1973).
// Generalized Feedback Shift Register Pseudorandom Number Algorithm
TEST(macro, macro_gfsr_file_access){
  int i;

  // Enable GFSR
  char *arg = (char *) malloc(3 * sizeof(char));
  arg[0] = '-';
  arg[1] = 'g';
  arg[2] = '\0';
  char **argv = (char **) malloc(10 * sizeof(char *));
  argv[9] = arg;

  resetgfsr();
  for(i = 0; i < P(); i++) EXPECT_EQ(resetgfsr_standard[i], rand_table[i]);

  EXPECT_EQ(GFSR_RESET_PERIOD(), 9800);
  EXPECT_EQ(GFSR_FLUSH_ITERATIONS(), 490000);
  EXPECT_EQ(GFSR_SBITS(), 32);

  EXPECT_EQ(intrand(), 1587561535);
  EXPECT_EQ(disrand(4, 8), 5);
  EXPECT_FLOAT_EQ(gfsr4(), 0.42877844);

  closegfsr();
  opengfsr();
  EXPECT_EQ(intrand(), 2036303646);
  EXPECT_EQ(disrand(4, 8), 6);
  EXPECT_FLOAT_EQ(gfsr4(), 0.77386665);

  resetgfsr();

  EXPECT_EQ(intrand(), 1587561535);
  EXPECT_EQ(disrand(4, 8), 5);
  EXPECT_FLOAT_EQ(gfsr4(), 0.42877844);

  resetgfsr();
  closegfsr();

  free(arg);
  free(argv);
}
