#https://www.jetbrains.com/help/clion/unit-testing-tutorial.html

project(Test)



enable_testing()

add_executable(
        tests
        tests.cpp
)
target_link_libraries(
        tests
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(tests)
