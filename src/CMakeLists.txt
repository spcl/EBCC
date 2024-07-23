cmake_minimum_required(VERSION 3.18)
project(h5z_j2k)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS_RELEASE "-O3")
set(BUILD_THIRDPARTY True)

set(CMAKE_POSITION_INDEPENDENT_CODE True)

set(SOURCES h5z_j2k.c)
add_library(h5z_j2k SHARED ${SOURCES})

add_subdirectory(hdf5)
add_subdirectory(zstd/build/cmake)
add_subdirectory(wavelib/src)
add_subdirectory(openjpeg)


include_directories(hdf5/src zstd/lib openjpeg/src/lib/openjp2 build/openjpeg/src/lib/openjp2 wavelib/header)
target_link_libraries(h5z_j2k libzstd_static openjp2_static wavelib hdf5)

# install(TARGETS h5z_j2k LIBRARY)