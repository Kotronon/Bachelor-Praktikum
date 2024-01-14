#https://www.jetbrains.com/help/clion/unit-testing-tutorial.html
#http://google.github.io/googletest/quickstart-cmake.html
#https://enccs.github.io/cmake-workshop/hello-ctest/
cmake_minimum_required(VERSION 3.14)

project(Tests)
# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
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

include(FetchContent)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.x
)
FetchContent_GetProperties(spdlog)
if (NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
    add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR})
endif ()

#find_package(spdlog REQUIRED)





enable_testing()
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)
#find_package(GTest REQUIRED)
include(GoogleTest)
include(CTest)
add_executable(tests tests.cpp ../src/ParticleContainer.cpp ../src/ParticleGenerator.cpp ../src/Particle.cpp ../src/calculations/PositionCalculator.cpp ../src/calculations/VelocityCalculator.cpp ../src/calculations/ForceCalculator.cpp ../src/LinkedCellContainer.cpp ../src/Thermostat.cpp)

FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master) # or whatever tag you want
FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
target_link_libraries(tests PUBLIC matplot)

target_link_libraries(tests PRIVATE GTest::gtest_main gtest_main spdlog::spdlog )
gtest_add_tests(TARGET tests)
gtest_discover_tests(tests)
#add_test((tests tests))
