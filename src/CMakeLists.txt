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

option(ENABLE_PERF "Enable Linux perf event collection around compression" OFF)
option(ENABLE_JXL "Enable JPEG XL backend" OFF)

option(ZSTD_BUILD_PROGRAMS "" OFF)
option(BUILD_CODEC "Build OPENJPEG CLI APPS" OFF)
add_subdirectory(zstd/build/cmake zstd)
add_subdirectory(openjpeg)

if(ENABLE_JXL)
  # Build libjxl from the bundled submodule with minimal options
  set(JPEGXL_ENABLE_TOOLS OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_MANPAGES OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_BENCHMARK OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_JNI OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_SJPEG OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_OPENEXR OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_PLUGINS OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_VIEWERS OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_FUZZERS OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_DEVTOOLS OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_JPEGLI OFF CACHE BOOL "" FORCE)
  set(JPEGXL_ENABLE_JPEGLI_LIBJPEG OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  add_subdirectory(jxl)
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/zstd/lib
  ${CMAKE_CURRENT_SOURCE_DIR}/zstd/build/cmake/lib
  ${CMAKE_CURRENT_SOURCE_DIR}/openjpeg/src/lib/openjp2
  ${CMAKE_CURRENT_BINARY_DIR}/openjpeg/src/lib/openjp2
  ${CMAKE_CURRENT_SOURCE_DIR}/spiht
  ${CMAKE_CURRENT_SOURCE_DIR}/log
)

set(EBCC_CORE_SOURCES
  ebcc_codec.c
  spiht/spiht_re.c
  log/log.c
)

add_library(h5z_ebcc SHARED h5z_ebcc.c ${EBCC_CORE_SOURCES})
add_library(h5z_ebcc_jxl SHARED h5z_ebcc.c ${EBCC_CORE_SOURCES})
add_library(ebcc STATIC ${EBCC_CORE_SOURCES})

set_property(TARGET h5z_ebcc PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
set_property(TARGET h5z_ebcc_jxl PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

# HDF5 filter 308 -> J2K default backend
# HDF5 filter 309 -> JXL default backend
# Both plugins can decode either backend by reading the serialized EBCC header.
target_compile_definitions(h5z_ebcc PRIVATE
  EBCC_HDF5_PLUGIN_J2K
)
target_compile_definitions(h5z_ebcc_jxl PRIVATE
  EBCC_HDF5_PLUGIN_JXL
)

target_link_libraries(h5z_ebcc PRIVATE libzstd_static openjp2_static)
target_link_libraries(h5z_ebcc_jxl PRIVATE libzstd_static openjp2_static)
target_link_libraries(ebcc PRIVATE libzstd_static openjp2_static)

if(ENABLE_JXL)
  message(STATUS "JXL support: ON")
  target_compile_definitions(h5z_ebcc PRIVATE ENABLE_JXL)
  target_compile_definitions(h5z_ebcc_jxl PRIVATE ENABLE_JXL)
  target_compile_definitions(ebcc PRIVATE ENABLE_JXL)

  target_include_directories(h5z_ebcc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/jxl/lib/include)
  target_include_directories(h5z_ebcc_jxl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/jxl/lib/include)
  target_include_directories(ebcc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/jxl/lib/include)

  target_link_libraries(h5z_ebcc PRIVATE jxl)
  target_link_libraries(h5z_ebcc_jxl PRIVATE jxl)
  target_link_libraries(ebcc PRIVATE jxl)

  if(TARGET jxl_threads)
    target_link_libraries(h5z_ebcc PRIVATE jxl_threads)
    target_link_libraries(h5z_ebcc_jxl PRIVATE jxl_threads)
    target_link_libraries(ebcc PRIVATE jxl_threads)
  endif()
else()
  message(STATUS "JXL support: OFF")
endif()

if(ENABLE_PERF)
  target_compile_definitions(h5z_ebcc PRIVATE ENABLE_PERF)
  target_compile_definitions(h5z_ebcc_jxl PRIVATE ENABLE_PERF)
  target_compile_definitions(ebcc PRIVATE ENABLE_PERF)
endif()

install(
  TARGETS h5z_ebcc h5z_ebcc_jxl ebcc
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
