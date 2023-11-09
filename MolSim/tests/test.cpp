
#include <gtest/gtest.h>


//@TODO write tests with TEST()

TEST(HELLOTEST, qualitycheck){
        EXPECT_EQ(7 * 6, 42);
}

int main(){
    return RUN_ALL_TESTS();
}