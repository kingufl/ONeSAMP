#include <gtest/gtest.h>
#include <stdio.h>

#ifndef GTEST_MAIN
#define GTEST_MAIN

/*! \brief Main method for unit tests for OneSamp.
 */
int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif
