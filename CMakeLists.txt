cmake_minimum_required(VERSION 3.14)

# define project name, version
project(PSEMolDyn_GroupH VERSION 0.0.1)



# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)


# collect all cpp files
file(GLOB_RECURSE MY_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    # header don't need to be included but this might be necessary for some IDEs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

# create make target
add_executable(MolSim ${MY_SRC}
        src/ParticleContainer.cpp
        src/ParticleContainer.h
        src/calculations/ForceCalculator.cpp
        src/calculations/ForceCalculator.h
        src/calculations/VelocityCalculator.cpp
        src/calculations/VelocityCalculator.h
        src/calculations/PositionCalculator.cpp
        src/calculations/PositionCalculator.h
        src/ParticleGenerator.cpp
        src/ParticleGenerator.h
        src/LinkedCellContainer.cpp
        src/LinkedCellContainer.h

        XML_input/inputnew-pskel.hxx
        XML_input/inputnew-pskel.cxx
        XML_input/inputnew-driver.cxx
        XML_input/inputnew-pimpl.cxx
        XML_input/inputnew-pimpl.hxx XML_input/inputnew-pimpl.cxx)


# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
            cxx_std_17
)

target_include_directories(MolSim
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
            xerces-c
)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:-Wall>
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
# Doxygen

# look for Doxygen package
OPTION(BUILD_DOC "Build documentation" on)
if(BUILD_DOC)
    find_package(Doxygen)

    if (DOXYGEN_FOUND)
        # configure doxygen configuration file
        set(doxyfile ${PROJECT_SOURCE_DIR}/Doxyfile)

        configure_file(${PROJECT_SOURCE_DIR}/Doxyfile @ONLY)
        message("Doxygen build started")

        add_custom_target(doc_doxygen
                COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                COMMAND firefox ${CMAKE_CURRENT_SOURCE_DIR}/doxys_documentation/html/index.html
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM )

    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
else(BUILD_DOC)
    message("Doxygen won't start")
endif (BUILD_DOC)



## logging
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

find_package(spdlog REQUIRED)


target_link_libraries(MolSim PRIVATE spdlog::spdlog)

# cmake --log-level=<level> or --loglevel=<level> with level = ERROR, WARNING, NOTICE, STATUS (default), VERBOSE, DEBUG, or TRACE.
#trys of setting log level
set( CMAKE_MESSAGE_LOG_LEVEL "STATUS" )
target_compile_definitions(MolSim PUBLIC OFF)
set(SPDLOG_INCLUDES_LEVEL "")
if(SPDLOG_SYSTEM_INCLUDES)
    set(SPDLOG_INCLUDES_LEVEL "SYSTEM")
    set(SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF)
endif()
target_include_directories(MolSim ${SPDLOG_INCLUDES_LEVEL} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

include(CTest)
add_subdirectory(tests)
