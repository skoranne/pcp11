////////////////////////////////////////////////////////////////////////////////
// File   : unit_test.cpp
// Author : Sandeep Koranne (C)
// Purpose: Notes on how to integrate C++ unit testing framework
// g++ -Wall -std=c++17 -I /scratch1/skoranne/CAD/include/ unit_test.cpp -L /scratch1/skoranne/CAD/lib64/ -lgtest -pthread
////////////////////////////////////////////////////////////////////////////////

#include <gtest/gtest.h>

// Function to be tested
int Add(int a, int b) {
    return a + b;
}

// Basic test case: TEST(TestSuiteName, TestName)
TEST(AdditionTest, PositiveNumbers) {
    EXPECT_EQ(Add(1, 2), 3); // Continues if failure occurs
    ASSERT_EQ(Add(10, 5), 15); // Stops if failure occurs
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
