#https://www.jetbrains.com/help/clion/unit-testing-tutorial.html

project(Test)
include(CTest)
enable_testing()
add_subdirectory(lib)
add_subdirectory(${CMAKE_BINARY_DIR}/googletest)
add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# adding the Google_Tests_run target
add_executable(Google_Tests_run tests.cpp)

# linking Google_Tests_run with DateConverter_lib which will be tested
target_link_libraries(Google_Tests_run src)

target_link_libraries(Google_Tests_run gtest gtest_main)
