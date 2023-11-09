
#include <gtest/gtest.h>


//@TODO write tests with TEST()
//To compile tests write cmake --build . in terminal and afterwarts ctest works

TEST(HELLOTEST, qualitycheck){
        EXPECT_EQ(7 * 6, 42);
}

int main(){
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}