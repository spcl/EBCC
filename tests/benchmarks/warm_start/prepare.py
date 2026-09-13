"""Build an isolated OpenJPEG warm-start probe; leave the submodule untouched."""

import argparse
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True,
                        help="Existing OpenJPEG build directory with opj_config*.h")
    parser.add_argument("--zstd-library", type=Path,
                        help="Build an isolated EBCC probe using this existing static zstd library")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[2]
    source = root / "src/openjpeg/src/lib/openjp2"
    dest = args.work_dir.resolve() / "openjp2"
    dest.mkdir(parents=True, exist_ok=True)
    for p in source.iterdir():
        if p.suffix in (".c", ".h"):
            shutil.copy2(p, dest / p.name)
    for name in ("opj_config.h", "opj_config_private.h"):
        shutil.copy2(args.config_dir / name, dest / name)

    p = dest / "tcd.c"
    text = p.read_text()
    start = "        /* FIXME _ProfStart(PGROUP_DC_SHIFT); */"
    end = "        /* FIXME _ProfStart(PGROUP_RATE); */"
    assert text.count(start) == text.count(end) == 1
    text = text.replace(start, "        if (p_tcd != ebcc_warm_tile) {\n" + start)
    text = text.replace(end, "        }\n\n" + end)
    include = '#include "opj_includes.h"'
    assert text.count(include) == 1
    text = text.replace(include, include + "\nextern _Thread_local opj_tcd_t *ebcc_warm_tile;")
    p.write_text(text)
    p = dest / "j2k.c"
    p.write_text(p.read_text() + "\n" + (here / "j2k_warm.inc").read_text())
    shutil.copy2(here / "probe.c", dest.parent / "probe.c")
    # Match the upstream library's core source list, excluding optional JPIP.
    names = "thread bio cio dwt event ht_dec image invert j2k jp2 mct mqc openjpeg opj_clock pi t1 t2 tcd tgt function_list opj_malloc sparse_array".split()
    cmake = ["cmake_minimum_required(VERSION 3.18)", "project(ebcc_warm_probe C)",
             "set(CMAKE_C_STANDARD 11)", "set(CMAKE_C_FLAGS_RELEASE \"-O3 -DNDEBUG\")",
             "add_library(warm_probe SHARED probe.c"]
    cmake += [f"  openjp2/{name}.c" for name in names]
    cmake += [")", "target_include_directories(warm_probe PRIVATE openjp2)",
              "find_package(Threads REQUIRED)",
              "target_link_libraries(warm_probe PRIVATE Threads::Threads m)"]
    if args.zstd_library:
        # Keep this historical prototype reproducible after production warm
        # start changes the current codec. Its reference source is fixed.
        text = subprocess.run(["git", "-C", str(root), "show", "399ade6:src/ebcc_codec.c"],
                              check=True, capture_output=True, text=True).stdout
        include = '#include "openjpeg.h"'
        declarations = '''
/* Temporary warm-start experiment; cache lifetime is one ebcc_encode call. */
extern void *probe_create(const uint16_t *, unsigned, unsigned);
extern void probe_destroy(void *);
extern int probe_encode(void *, float, float, int);
extern const void *probe_data(void *);
extern size_t probe_size(void *);
static _Thread_local void *ebcc_warm_session;
'''
        text = text.replace(include, include + declarations)
        marker = "    size_t n_tiles = image_dims[0] / tile_dims[0];"
        start = text.index("void j2k_encode_internal_mode(")
        pos = text.index(marker, start)
        warm = '''    const char *warm = getenv("EBCC_WARM_START");
    if (!lossless && image_dims[0] == tile_dims[0] && warm && strcmp(warm, "1") == 0) {
        if (!ebcc_warm_session)
            ebcc_warm_session = probe_create(data, image_dims[1], image_dims[0]);
        if (!ebcc_warm_session || !probe_encode(ebcc_warm_session, psnr,
                base_cr / EBCC_J2K_RATE_DIVISOR, TRUE)) abort();
        write_to_buffer_stream((void *)probe_data(ebcc_warm_session),
                               probe_size(ebcc_warm_session), codec_data_buffer);
        return;
    }
'''
        text = text[:pos] + warm + text[pos:]
        cleanup = "    if (!const_field) free(scaled_data);"
        assert text.count(cleanup) == 1
        text = text.replace(cleanup, "    probe_destroy(ebcc_warm_session);\n"
                            "    ebcc_warm_session = NULL;\n" + cleanup)
        (dest.parent / "ebcc_codec_probe.c").write_text(text)
        cmake += [f'add_library(ebcc_probe SHARED ebcc_codec_probe.c "{root}/src/spiht/spiht_re.c" "{root}/src/log/log.c")',
                  f'target_include_directories(ebcc_probe PRIVATE openjp2 "{root}/src" "{root}/src/spiht" "{root}/src/log" "{root}/src/zstd/lib")',
                  f'target_link_libraries(ebcc_probe PRIVATE warm_probe "{args.zstd_library.resolve()}")']
    (dest.parent / "CMakeLists.txt").write_text("\n".join(cmake) + "\n")
    subprocess.run(["cmake", "-S", str(dest.parent), "-B", str(dest.parent / "build"),
                    "-DCMAKE_BUILD_TYPE=Release"], check=True)
    subprocess.run(["cmake", "--build", str(dest.parent / "build"), "-j", "4"], check=True)


if __name__ == "__main__":
    main()
