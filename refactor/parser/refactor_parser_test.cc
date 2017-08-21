#include <stdlib.h>
#include <stdio.h>
#include <gtest/gtest.h>

extern "C"{
#include "../macro/refactor_macro.h"
#include "refactor_parser.h"
}

#define test_file

// Test parser utility functions.
TEST(parser, parseGenotypeToken){
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  int r1;
  int r2;
  int *returngene1ptr = &r1;
  int *returngene2ptr = &r2;
  initializeParserQueue();
  enqueueParserToken((int) '\n');
  enqueueParserToken((int) '0');
  enqueueParserToken((int) '0');
  enqueueParserToken((int) '0');
  enqueueParserToken((int) '0');
  enqueueParserToken((int) '\n');
  parseGenotypeToken(returngene1ptr, returngene2ptr);
  EXPECT_EQ(*returngene1ptr, 0);
  EXPECT_EQ(*returngene2ptr, 0);

  enqueueParserToken((int) '0');
  enqueueParserToken((int) '1');
  enqueueParserToken((int) '0');
  enqueueParserToken((int) '1');
  enqueueParserToken((int) '\n');
  parseGenotypeToken(returngene1ptr, returngene2ptr);
  EXPECT_EQ(*returngene1ptr, 1);
  EXPECT_EQ(*returngene2ptr, 1);

  enqueueParserToken((int) '9');
  enqueueParserToken((int) '9');
  enqueueParserToken((int) '9');
  enqueueParserToken((int) '9');
  enqueueParserToken((int) '9');
  enqueueParserToken((int) '9');
  enqueueParserToken((int) '\n');
  parseGenotypeToken(returngene1ptr, returngene2ptr);
  EXPECT_EQ(*returngene1ptr, 999);
  EXPECT_EQ(*returngene2ptr, 999);

  enqueueParserToken((int) '1');
  enqueueParserToken((int) '4');
  enqueueParserToken((int) '2');
  enqueueParserToken((int) '8');
  enqueueParserToken((int) '5');
  enqueueParserToken((int) '7');
  enqueueParserToken((int) '\n');
  parseGenotypeToken(returngene1ptr, returngene2ptr);
  EXPECT_EQ(*returngene1ptr, 142);
  EXPECT_EQ(*returngene2ptr, 857);

}

// Test parser utility functions
TEST(parser, matchPOP){
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  initializeParserQueue();
  enqueueParserToken((int) '\n');
  enqueueParserToken((int) 'P');
  enqueueParserToken((int) 'O');
  enqueueParserToken((int) 'P');
  enqueueParserToken((int) '\n');
  EXPECT_EQ(matchPop(), zeroOne(TRUE));
  enqueueParserToken((int) 'p');
  EXPECT_EQ(matchPop(), zeroOne(FALSE));
  enqueueParserToken((int) 'O');
  enqueueParserToken((int) 'P');
  enqueueParserToken((int) '\n');
  EXPECT_EQ(matchPop(), zeroOne(TRUE));
}

TEST(parser, matchCharacter){

  // To test the match functions, use a test string.

  // All of the test in the loop correspond to the characters in that string.
  char testString[] = "PpOo0123456789 \t,\r\n";
  int length = strlen(testString);
  int i;
  int j;
  for(i = 0; i < length; i++){
    int j = (int) testString[i];
    EXPECT_EQ(matchSlashN(&j), zeroOne(i == 18));
    EXPECT_EQ(matchSlashR(&j), zeroOne(i == 17));
    EXPECT_EQ(matchComma(&j), zeroOne(i == 16));
    EXPECT_EQ(matchP(&j), zeroOne(i == 0 || i == 1));
    EXPECT_EQ(matchO(&j), zeroOne(i == 2 || i == 3));
    EXPECT_EQ(matchDigit(&j), zeroOne(i > 3 && i < 14));
    EXPECT_EQ(matchEndOfLine(&j), zeroOne(i == 17 || i == 18));
    EXPECT_EQ(matchWhitespace(&j),
              zeroOne(i == 14 || i == 15 || i == 17 || i == 18));
    EXPECT_EQ(matchHardDelimiter(&j), zeroOne(i > 15));
    EXPECT_EQ(matchSoftDelimiter(&j), zeroOne(i > 13));
    EXPECT_EQ(matchEOF(&j), zeroOne(FALSE));
  }

  // Test EOF
  j = EOF;
  EXPECT_EQ(matchSlashN(&j), zeroOne(FALSE));
  EXPECT_EQ(matchSlashR(&j), zeroOne(FALSE));
  EXPECT_EQ(matchComma(&j), zeroOne(FALSE));
  EXPECT_EQ(matchP(&j), zeroOne(FALSE));
  EXPECT_EQ(matchO(&j), zeroOne(FALSE));
  EXPECT_EQ(matchDigit(&j), zeroOne(FALSE));
  EXPECT_EQ(matchEndOfLine(&j), zeroOne(FALSE));
  EXPECT_EQ(matchWhitespace(&j), zeroOne(FALSE));
  EXPECT_EQ(matchHardDelimiter(&j), zeroOne(TRUE));
  EXPECT_EQ(matchSoftDelimiter(&j), zeroOne(TRUE));
  EXPECT_EQ(matchEOF(&j), zeroOne(TRUE));
}

TEST(parser, zeroOne){
  EXPECT_EQ(zeroOne(FALSE), 0);
  EXPECT_EQ(zeroOne(TRUE), 1);
}

// Test parser on entire file, using newlines to delineate text.
TEST(parser, parseA){
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  int syntax_results[2] = {0, 0};
  parse(TRUE, "../parser/testA.txt", syntax_results);
  ASSERT_EQ(syntax_results[0], 3);
  ASSERT_EQ(syntax_results[1], 3);
}

// Test parser on entire file, using newlines to delineate text.
TEST(parser, parseB){
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  int syntax_results[2] = {0, 0};
  parse(TRUE, "../parser/testB.txt", syntax_results);
  ASSERT_EQ(syntax_results[0], 3);
  ASSERT_EQ(syntax_results[1], 3);
}

// Test parser on entire file with syntax error.
TEST(parser, parseC){
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  int syntax_results[2] = {0, 0};
  ASSERT_DEATH(parse(TRUE, "../parser/testC.txt", syntax_results), "ONESAMP PARSE ERROR, line 6, column 14 \nInvalid SNP: must be a number that is 0, 1, 2, 3, or 4.\nExiting...\n");
}

// Test parser on entire file.
// TODO Still possibly an error with blank lines? Why doesn't it show up on the first line.
TEST(parser, parseD){
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  int syntax_results[2] = {0, 0};
  ASSERT_DEATH(parse(TRUE, "../parser/testD.txt", syntax_results), "ONESAMP PARSE ERROR, line 8, column 0 \nGenePop 4.0 format does not allow blank lines.\nExiting...\n");
}

// Test parser on entire file.
TEST(parser, parseE){
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  int syntax_results[2] = {0, 0};
  parse(TRUE, "../parser/testE.txt", syntax_results);
  ASSERT_EQ(syntax_results[0], 4);
  ASSERT_EQ(syntax_results[1], 3);
}

// Test parser on entire file.
TEST(parser, parseF){
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  char *argv_name = (char *) "parser_test";
  char **argv = &argv_name;
  int syntax_results[2] = {0, 0};
  parse(TRUE, "../parser/testF.txt", syntax_results);
  ASSERT_EQ(syntax_results[0], 3);
  ASSERT_EQ(syntax_results[1], 3);
}
