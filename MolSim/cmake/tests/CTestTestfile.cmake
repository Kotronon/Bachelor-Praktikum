# CMake generated Testfile for 
# Source directory: /home/kathi/Bachelor-Praktikum/MolSim/tests
# Build directory: /home/kathi/Bachelor-Praktikum/MolSim/cmake/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
include("/home/kathi/Bachelor-Praktikum/MolSim/cmake/tests/tests[1]_include.cmake")
add_test(old:ParticleContainerTest.ParticleContainer "/home/kathi/Bachelor-Praktikum/MolSim/cmake/tests/tests" "--gtest_filter=ParticleContainerTest.ParticleContainer")
set_tests_properties(old:ParticleContainerTest.ParticleContainer PROPERTIES  _BACKTRACE_TRIPLES "/usr/share/cmake-3.22/Modules/GoogleTest.cmake;400;add_test;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;39;gtest_add_tests;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;0;")
add_test(old:PositionTest.stroemerVelvet "/home/kathi/Bachelor-Praktikum/MolSim/cmake/tests/tests" "--gtest_filter=PositionTest.stroemerVelvet")
set_tests_properties(old:PositionTest.stroemerVelvet PROPERTIES  _BACKTRACE_TRIPLES "/usr/share/cmake-3.22/Modules/GoogleTest.cmake;400;add_test;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;39;gtest_add_tests;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;0;")
add_test(old:VelocityTest.BrownianMotionInitialization "/home/kathi/Bachelor-Praktikum/MolSim/cmake/tests/tests" "--gtest_filter=VelocityTest.BrownianMotionInitialization")
set_tests_properties(old:VelocityTest.BrownianMotionInitialization PROPERTIES  _BACKTRACE_TRIPLES "/usr/share/cmake-3.22/Modules/GoogleTest.cmake;400;add_test;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;39;gtest_add_tests;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;0;")
add_test(old:VelocityTest.stroemerVelvet "/home/kathi/Bachelor-Praktikum/MolSim/cmake/tests/tests" "--gtest_filter=VelocityTest.stroemerVelvet")
set_tests_properties(old:VelocityTest.stroemerVelvet PROPERTIES  _BACKTRACE_TRIPLES "/usr/share/cmake-3.22/Modules/GoogleTest.cmake;400;add_test;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;39;gtest_add_tests;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;0;")
add_test(old:ForceTest.SimpleForceCalculation "/home/kathi/Bachelor-Praktikum/MolSim/cmake/tests/tests" "--gtest_filter=ForceTest.SimpleForceCalculation")
set_tests_properties(old:ForceTest.SimpleForceCalculation PROPERTIES  _BACKTRACE_TRIPLES "/usr/share/cmake-3.22/Modules/GoogleTest.cmake;400;add_test;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;39;gtest_add_tests;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;0;")
add_test(old:ForceTest.LennardJonesForce "/home/kathi/Bachelor-Praktikum/MolSim/cmake/tests/tests" "--gtest_filter=ForceTest.LennardJonesForce")
set_tests_properties(old:ForceTest.LennardJonesForce PROPERTIES  _BACKTRACE_TRIPLES "/usr/share/cmake-3.22/Modules/GoogleTest.cmake;400;add_test;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;39;gtest_add_tests;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;0;")
add_test(tests "tests")
set_tests_properties(tests PROPERTIES  _BACKTRACE_TRIPLES "/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;41;add_test;/home/kathi/Bachelor-Praktikum/MolSim/tests/CMakeLists.txt;0;")
subdirs("../_deps/googletest-build")
