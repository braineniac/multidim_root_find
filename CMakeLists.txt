project (newton-raphson)

SET(CMAKE_C_COMPILER /bin/clang)

cmake_minimum_required (VERSION 2.6)
find_package(GSL REQUIRED)

SET( GSL_C_FLAGS $(gsl-config --cflags))
SET( GSL_LIB_FLAGS $(gsl-config --libs))

SET( CMAKE_C_FLAGS ${GSL_C_FLAGS})
SET( CMAKE_VERBOSE_MAKEFILE on)
SET( CMAKE_EXECUTABLE_SUFFIX ${GSL_LIB_FLAGS})

add_executable(newton-raphson newton-raphson.c)
