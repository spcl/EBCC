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

option(ENABLE_PERF "Enable Linux perf event collection around compression" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE True)

if(WIN32)
  # OpenJPEG exposes only the "openjp2" target on Windows; force a static build
  # so EBCC can link it without pulling in the shared variant.
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
  set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)

  # zstd builds both shared and static variants by default; disable the shared
  # target on Windows to avoid unnecessary shared-library build failures.
  set(ZSTD_BUILD_SHARED OFF CACHE BOOL "Build zstd shared libraries" FORCE)
  set(ZSTD_BUILD_STATIC ON CACHE BOOL "Build zstd static libraries" FORCE)
  set(OPENJPEG_LINK_TARGET openjp2)
else()
  set(OPENJPEG_LINK_TARGET openjp2_static)
endif()

add_library(h5z_ebcc SHARED h5z_ebcc.c ebcc_codec.c spiht/spiht_re.c log/log.c)
set_property(TARGET h5z_ebcc PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
if(WIN32)
  set_property(TARGET h5z_ebcc PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

# Add static EBCC library for Rust bindings
add_library(ebcc STATIC ebcc_codec.c spiht/spiht_re.c log/log.c)

option(ZSTD_BUILD_PROGRAMS "" OFF)
option(BUILD_CODEC "Build OPENJPEG CLI APPS" OFF)
add_subdirectory(zstd/build/cmake zstd)
add_subdirectory(openjpeg)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/zstd/lib
  ${CMAKE_CURRENT_SOURCE_DIR}/zstd/build/cmake/lib
  ${CMAKE_CURRENT_SOURCE_DIR}/openjpeg/src/lib/openjp2
  ${CMAKE_CURRENT_BINARY_DIR}/openjpeg/src/lib/openjp2
  ${CMAKE_CURRENT_SOURCE_DIR}/spiht
  ${CMAKE_CURRENT_SOURCE_DIR}/log
)

target_link_libraries(h5z_ebcc PRIVATE libzstd_static ${OPENJPEG_LINK_TARGET})
target_link_libraries(ebcc PRIVATE libzstd_static ${OPENJPEG_LINK_TARGET})

if(ENABLE_PERF)
    target_compile_definitions(h5z_ebcc PRIVATE ENABLE_PERF)
    target_compile_definitions(ebcc PRIVATE ENABLE_PERF)
endif()

install(
    TARGETS h5z_ebcc ebcc
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
