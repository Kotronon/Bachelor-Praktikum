#https://www.jetbrains.com/help/clion/unit-testing-tutorial.html
#http://google.github.io/googletest/quickstart-cmake.html
#https://enccs.github.io/cmake-workshop/hello-ctest/
cmake_minimum_required(VERSION 3.17)

project(Tests)
# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${GTEST_INCLUDE_DIRS})
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
FETCHCONTENT_GETPROPERTIES(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BUILD_DIR} ${googletest_BINARY_DIR})
endif ()
find_package(GTest REQUIRED)
include_directories(SYSTEM ${GMOCK_DIR}/include ${GMOCK_DIR}/gtest/include)
enable_testing()
set(gmock_build_include_dirs
        "${gmock_SOURCE_DIR}/include"
        "${gmock_SOURCE_DIR}"
        "${gtest_SOURCE_DIR}/include"
        # This directory is needed to build directly from Google Test sources.
        "${gtest_SOURCE_DIR}")
include_directories(${gmock_build_include_dirs})
add_executable(tests tests.cpp ../src/ParticleContainer.cpp ../src/Particle.cpp ../src/calculations/PositionCalculator.cpp ../src/calculations/VelocityCalculator.cpp ../src/calculations/ForceCalculator.cpp )
target_link_libraries(tests GTest::GTest GTest::Main ${GTEST_BOTH_LIBRARIES})
target_link_libraries(tests gtest gtest_main)
gtest_add_tests(TARGET tests)
gtest_discover_tests(tests)
#add_test(tests tests)
