project (multidim-root-find)

SET(CMAKE_C_COMPILER /bin/clang)

cmake_minimum_required (VERSION 2.6)
find_package(GSL REQUIRED)

#SET( GSL_C_FLAGS $(gsl-config --cflags))
#SET( GSL_LIB_FLAGS $(gsl-config --libs))

SET( DEFAULT_FLAGS "-Wall")

SET( CMAKE_C_FLAGS ${DEFAULT_FLAGS} ${GSL_C_FLAGS})
SET( CMAKE_VERBOSE_MAKEFILE on)
#SET( CMAKE_EXECUTABLE_SUFFIX ${GSL_LIB_FLAGS})

SET (SOURCE_FILES multidim-root-find.c)

add_executable(multidim-root-find ${SOURCE_FILES})
target_link_libraries(multidim-root-find m)
#target_link_libraries(newton_raphson dl)
target_link_libraries(multidim-root-find GSL::gsl GSL::gslcblas)
