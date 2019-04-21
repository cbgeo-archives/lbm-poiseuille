project(lbm LANGUAGES CXX)

# Require C++14-compliant compiler; only available for CMAke v. 3.1 and up
set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.1)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Debug")
  #SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GPU specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# Intel specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Clang specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

option(LBM_BUILD_TESTING "enable testing for lbm" ON)


# Include directories
include_directories(BEFORE
  ${lbm_SOURCE_DIR}/include/
  ${lbm_SOURCE_DIR}/external/
)

# lbm executable
SET(lbm_src
  ${lbm_SOURCE_DIR}/src/lbgrid.cc
)

add_library(llbm SHARED ${lbm_src})

add_executable(lbm ${lbm_SOURCE_DIR}/src/main.cc)
target_link_libraries(lbm llbm)

# Unit test
if(LBM_BUILD_TESTING)
  SET(test_src
    ${lbm_SOURCE_DIR}/tests/test_main.cc
    ${lbm_SOURCE_DIR}/tests/lbgrid_test.cc
    )
  add_executable(lbmtest ${test_src})
  target_link_libraries(lbmtest llbm)
  add_test(NAME lbmtest COMMAND $<TARGET_FILE:lbmtest>)
  enable_testing()
endif()
