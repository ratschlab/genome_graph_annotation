cmake_minimum_required (VERSION 3.6.1)
project(AnnotationSchemes)

enable_testing()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif()

string(APPEND CMAKE_CXX_FLAGS " \
  -Wall -Wextra -Werror \
  -msse4.2 \
  -fopenmp -D_THREAD_SAFE -pthread \
") # -DDBGDEBUG -Wconversion -Wsign-conversion -Werror=shadow=compatible-local

if(NOT DEFINED WITH_AVX)
set(WITH_AVX TRUE)
endif()
if(WITH_AVX)
string(APPEND CMAKE_CXX_FLAGS " -mavx -mavx2 -mfma -mbmi -mbmi2")
endif()

if(WITH_MMAP)
message(STATUS "Enabling boost iostreams for mmapped reader")
string(APPEND CMAKE_CXX_FLAGS " -D_USE_MMAP")
set(METALIBS ${METALIBS} -lboost_iostreams)
endif()

if(WITH_DEATH_TESTS)
message(STATUS "Enabling death tests in gtest")
string(APPEND CMAKE_CXX_FLAGS " -D_DEATH_TEST")
endif()

#-------------------
# KMC k-mer counter
#-------------------
set(KMC_MAIN_DIR "${PROJECT_SOURCE_DIR}/external-libraries/KMC/kmer_counter")
if(NOT DEFINED BUILD_KMC)
set(BUILD_KMC TRUE)
endif()
configure_file(CMakeListsKMC.txt.in KMC/CMakeLists.txt @ONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/KMC
)
if(result)
  message(FATAL_ERROR "CMake step for KMC failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release -- -j
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/KMC
)
if(result)
  message(FATAL_ERROR "Build step for KMC failed: ${result}")
endif()


# Profile build type
set(CMAKE_CXX_FLAGS_PROFILE "-pg -DNDEBUG -O3 -g")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg -g")

# Debug build type
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")

set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_policy(SET CMP0015 NEW)

link_directories(
  ${CMAKE_BINARY_DIR}/KMC
  ${CMAKE_BINARY_DIR}/lib
  external-libraries/libmaus2/lib
  external-libraries/sdsl-lite/lib
)

include_directories(
  external-libraries/KMC/kmc_api
  external-libraries/libmaus2/include
  external-libraries/sdsl-lite/include
  external-libraries/eigen-git-mirror
  external-libraries/caches/include
  external-libraries/sparsepp
  external-libraries/hopscotch-map/include
  external-libraries/ordered-map/include
  external-libraries/binrel_wt/include
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/external-libraries/binrel_wt/cmake/modules")
add_subdirectory(
  external-libraries/binrel_wt
)

set(FOLLY_LIBRARIES FOLLY_LIBRARIES-NOTFOUND)
find_library(FOLLY_LIBRARIES
  NAMES folly
  HINTS
  ENV LD_LIBRARY_PATH
  PATHS
  ~/.linuxbrew/lib/
)
if(FOLLY_LIBRARIES)
  string(APPEND CMAKE_CXX_FLAGS " -D_USE_FOLLY")
  set(METALIBS -lfolly -lglog)
else()
  message(WARNING "Folly was not found.\n"
  "Install folly to reduce RAM required by the row-compressed annotator.")
endif()

file(GLOB src_files "*.cpp"
  "common/*.cpp"
  "common/kmc_counter/*.cpp"
  "dbg_hash/*.cpp"
  "annotation/*.cpp"
  "annotation/*/*.cpp"
)

list(FILTER src_files EXCLUDE REGEX ".*\\._.*")
list(FILTER src_files EXCLUDE REGEX ".*main.cpp")

add_library(annographlibs STATIC ${src_files})
add_executable(annograph "main.cpp")

set(METALIBS ${METALIBS}
  -lKMC
  -lhts -lz -lbz2
  -lmaus2 -lmaus2_simd_align_128 -lmaus2_simd_align_256
  -lsdsl -ldivsufsort -ldivsufsort64
  -lssl -lcrypto -llzma
  -lboost_filesystem -lboost_system
  -lbrwt
)

if(BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -static)
  set(METALIBS ${METALIBS} -lcurl -lssl)
endif()

target_link_libraries(
  annographlibs ${METALIBS}
)

target_link_libraries(
  annograph annographlibs ${METALIBS}
)


#-------------------
# Experiments
#-------------------
file(GLOB experiment_files "experiments/*.cpp")
list(FILTER experiment_files EXCLUDE REGEX ".*\\._.*")

include_directories(
  .
  common
  common/kmc_counter
  dbg_hash
  ${ANNO_DIRS}
  experiments
  external-libraries/tclap/include
)

link_directories(
  external-libraries/libmaus2/lib
  external-libraries/sdsl-lite/lib
  external-libraries/binrel_wt/lib
)
add_executable(run_experiments ${experiment_files})

target_link_libraries(
  run_experiments annographlibs ${METALIBS}
)


#-------------------
# Unit Tests
#-------------------
# Download and unpack googletest at configure time
if(FALSE) # There is one gtest version attached with sdsl
if(NOT EXISTS ${CMAKE_BINARY_DIR}/googletest-download)
  configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(
  ${CMAKE_BINARY_DIR}/googletest-src
  ${CMAKE_BINARY_DIR}/googletest-build
  EXCLUDE_FROM_ALL
)
endif()

file(GLOB unit_tests_files "tests/*.cpp")
list(FILTER unit_tests_files EXCLUDE REGEX ".*\\._.*")

if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "${curdir}")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(ANNO_DIRS "${PROJECT_SOURCE_DIR}/annotation")

include_directories(
  .
  common
  common/kmc_counter
  dbg_hash
  ${ANNO_DIRS}
  experiments
)

link_directories(
  ${CMAKE_BINARY_DIR}
)

add_executable(unit_tests ${unit_tests_files})

target_link_libraries(unit_tests gtest_main gtest
  annographlibs
  ${METALIBS}
)

add_test(NAME unit_tests COMMAND unit_tests)