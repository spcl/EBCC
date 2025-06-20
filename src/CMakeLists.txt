cmake_minimum_required(VERSION 3.18)
project(h5z_j2k)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_DEBUG "-g -O2 -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O2")
set(BUILD_THIRDPARTY True)

set(CMAKE_POSITION_INDEPENDENT_CODE True)


add_library(h5z_j2k SHARED h5z_j2k.c spiht/spiht_re.c log/log.c)
set_property(TARGET h5z_j2k PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

option(ZSTD_BUILD_PROGRAMS "" OFF)
option(BUILD_CODEC "Build OPENJPEG CLI APPS" OFF)
add_subdirectory(zstd/build/cmake)
add_subdirectory(openjpeg)

find_package(HDF5 REQUIRED)

include_directories(${HDF5_INCLUDE_DIRS} zstd/lib openjpeg/src/lib/openjp2 ${CMAKE_CURRENT_BINARY_DIR}/openjpeg/src/lib/openjp2 spiht log)
target_link_libraries(h5z_j2k PRIVATE libzstd_static openjp2_static ${HDF5_LIBRARIES})

install(
    TARGETS h5z_j2k
    LIBRARY DESTINATION lib
)