cmake_minimum_required(VERSION 3.18)
project(h5z_ebcc)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  set(CMAKE_C_FLAGS_RELEASE "/O2")
  set(CMAKE_CXX_FLAGS_RELEASE "/O2")
  set(CMAKE_C_FLAGS_DEBUG "/Od /Zi /DDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi")
else()
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_DEBUG "-g -O2 -DDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O2")
endif()
set(BUILD_THIRDPARTY True)

set(CMAKE_POSITION_INDEPENDENT_CODE True)


add_library(h5z_ebcc SHARED h5z_ebcc.c ebcc_codec.c spiht/spiht_re.c log/log.c)
set_property(TARGET h5z_ebcc PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

# Add static EBCC library for Rust bindings
add_library(ebcc STATIC ebcc_codec.c spiht/spiht_re.c log/log.c)
set_property(TARGET ebcc PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

option(ZSTD_BUILD_PROGRAMS "" OFF)
option(BUILD_CODEC "Build OPENJPEG CLI APPS" OFF)
add_subdirectory(zstd/build/cmake zstd)
add_subdirectory(openjpeg)

include_directories(zstd/lib openjpeg/src/lib/openjp2 ${CMAKE_CURRENT_BINARY_DIR}/openjpeg/src/lib/openjp2 spiht log)
target_link_libraries(h5z_ebcc PRIVATE libzstd_static openjp2_static)
target_link_libraries(ebcc PRIVATE libzstd_static openjp2_static)

install(
    TARGETS h5z_ebcc ebcc
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)