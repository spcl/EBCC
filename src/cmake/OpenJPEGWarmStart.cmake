# EBCC's private OpenJPEG extension. Keep the upstream checkout pristine and
# compile every core translation unit against the same extended private header.
set(_opj_source "${CMAKE_CURRENT_SOURCE_DIR}/openjpeg/src/lib/openjp2")
set(_opj_copy "${CMAKE_CURRENT_BINARY_DIR}/ebcc_openjpeg")
file(MAKE_DIRECTORY "${_opj_copy}")
file(GLOB _opj_files CONFIGURE_DEPENDS "${_opj_source}/*.c" "${_opj_source}/*.h")
foreach(_file IN LISTS _opj_files)
  get_filename_component(_name "${_file}" NAME)
  configure_file("${_file}" "${_opj_copy}/${_name}" COPYONLY)
endforeach()

# Fail clearly if an upstream update moves a hook, rather than silently building
# a partially patched encoder. Each hook must occur exactly once.
function(ebcc_opj_replace file before after)
  file(READ "${file}" _text)
  string(FIND "${_text}" "${before}" _position)
  if(_position LESS 0)
    message(FATAL_ERROR "OpenJPEG warm-start hook missing in ${file}")
  endif()
  string(LENGTH "${before}" _length)
  math(EXPR _next "${_position} + ${_length}")
  string(SUBSTRING "${_text}" ${_next} -1 _tail)
  string(FIND "${_tail}" "${before}" _duplicate)
  if(NOT _duplicate LESS 0)
    message(FATAL_ERROR "OpenJPEG warm-start hook is ambiguous in ${file}")
  endif()
  string(REPLACE "${before}" "${after}" _text "${_text}")
  file(WRITE "${file}" "${_text}")
endfunction()

ebcc_opj_replace("${_opj_copy}/tcd.h"
  "} opj_tcd_t;"
  "    /* EBCC: prepared Tier-1 state, owned by one encoder session. */\n    OPJ_BOOL ebcc_reuse_tier1;\n} opj_tcd_t;")
ebcc_opj_replace("${_opj_copy}/tcd.c"
  "        /* FIXME _ProfStart(PGROUP_DC_SHIFT); */"
  "        if (!p_tcd->ebcc_reuse_tier1) {\n        /* FIXME _ProfStart(PGROUP_DC_SHIFT); */")
ebcc_opj_replace("${_opj_copy}/tcd.c"
  "        /* FIXME _ProfStart(PGROUP_RATE); */"
  "        }\n\n        /* FIXME _ProfStart(PGROUP_RATE); */")
file(APPEND "${_opj_copy}/j2k.c" "\n#include \"ebcc_openjpeg_bridge.inc\"\n")

get_target_property(_opj_sources ${OPENJPEG_LINK_TARGET} SOURCES)
set(_opj_patched_sources)
foreach(_file IN LISTS _opj_sources)
  string(REPLACE "${_opj_source}/" "${_opj_copy}/" _patched "${_file}")
  list(APPEND _opj_patched_sources "${_patched}")
endforeach()
set_property(TARGET ${OPENJPEG_LINK_TARGET} PROPERTY SOURCES "${_opj_patched_sources}")
target_include_directories(${OPENJPEG_LINK_TARGET} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
# The extended private layout must never bind to another OpenJPEG loaded by
# the host process. OPJ_STATIC hides its public entry points on ELF/Mach-O;
# the visibility preset also hides internal functions and the EBCC bridge.
target_compile_definitions(${OPENJPEG_LINK_TARGET} PUBLIC OPJ_STATIC)
set_target_properties(${OPENJPEG_LINK_TARGET} PROPERTIES C_VISIBILITY_PRESET hidden)
