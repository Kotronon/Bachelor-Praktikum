
#include <gtest/gtest.h>


//@TODO write tests with TEST()

TEST(SIMPLETESTS,Multiplication){
    EXPECT_EQ(7*6,42);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
