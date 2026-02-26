import os
import re
import subprocess
import sys
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        ext_path = Path(self.get_ext_fullpath(ext.name)).resolve()

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = 'Debug' if debug else 'Release'
        build_temp = (Path(self.build_temp) / ext.name).resolve()
        build_temp.mkdir(parents=True, exist_ok=True)
        lib_output_dir = build_temp / "lib"
        lib_output_dir.mkdir(parents=True, exist_ok=True)
        lib_output_dir_str = str(lib_output_dir)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={lib_output_dir_str}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]

        # Always build with JXL enabled.
        enable_jxl = True
        cmake_args.append("-DENABLE_JXL=ON")
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        if self.compiler.compiler_type == "msvc":
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", "x64"]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={lib_output_dir_str}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        # Check if git submodules are initialized
        src_dir = Path(ext.sourcedir) / "src"
        openjpeg_dir = src_dir / "openjpeg"
        zstd_dir = src_dir / "zstd"
        jxl_dir = src_dir / "jxl"

        need_submodule_init = (
            not (openjpeg_dir / "CMakeLists.txt").exists()
            or not (zstd_dir / "CMakeLists.txt").exists()
        )
        if enable_jxl and not (jxl_dir / "CMakeLists.txt").exists():
            need_submodule_init = True

        def _dir_has_entries(path: Path) -> bool:
            return path.exists() and any(path.iterdir())

        if need_submodule_init:
            print("Initializing git submodules...")
            subprocess.run(
                ["git", "submodule", "update", "--init", "--recursive"],
                cwd=ext.sourcedir,
                check=True
            )
        elif enable_jxl:
            # Ensure jxl's own nested submodules (highway, brotli, etc.) are initialized
            jxl_highway = jxl_dir / "third_party" / "highway"
            jxl_brotli = jxl_dir / "third_party" / "brotli"
            if not _dir_has_entries(jxl_highway) or not _dir_has_entries(jxl_brotli):
                print("Initializing jxl nested submodules...")
                subprocess.run(
                    ["git", "submodule", "update", "--init", "--recursive"],
                    cwd=ext.sourcedir,
                    check=True
                )

        # Run CMake configure
        subprocess.run(
            ["cmake", str(src_dir)] + cmake_args, cwd=build_temp, check=True
        )
        
        # Run CMake build
        subprocess.run(
            ["cmake", "--build", "."] + build_args,
            cwd=build_temp,
            check=True,
        )

        src_lib_dir = Path(ext.sourcedir) / "ebcc"
        src_lib_dir.mkdir(parents=True, exist_ok=True)
        
        # Find and copy the built library to where setuptools expects it
        lib_name = "libh5z_ebcc.so"
        extra_plugin_name = "libh5z_ebcc_jxl.so"
        if sys.platform == "win32":
            lib_name = "h5z_ebcc.dll"
            extra_plugin_name = "h5z_ebcc_jxl.dll"
        elif sys.platform == "darwin":
            lib_name = "libh5z_ebcc.dylib"
            extra_plugin_name = "libh5z_ebcc_jxl.dylib"

        # Look for the library produced by CMake
        search_roots = [
            build_temp / "lib",
            build_temp / "openjpeg" / "bin",
            build_temp,
        ]
        candidates = [root / lib_name for root in search_roots] + [build_temp / "lib" / ext_path.name]
        built_lib = next((p for p in candidates if p.exists()), None)

        print(f"Built library candidates: {[str(p) for p in candidates]}")
        
        if built_lib:
            import shutil
            # Copy to where setuptools expects the built extension
            ext_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(built_lib, ext_path)
            print(f"Copied {built_lib.name} to {ext_path}")

            # Always copy both HDF5 plugin binaries into ebcc/
            plugin_names = (lib_name, extra_plugin_name)
            for plugin_name in plugin_names:
                plugin_candidates = [root / plugin_name for root in search_roots]
                plugin_path = next((p for p in plugin_candidates if p.exists()), None)
                if plugin_path is None:
                    print(f"Warning: Could not find plugin binary {plugin_name}")
                    continue
                dst = src_lib_dir / plugin_path.name
                shutil.copy2(plugin_path, dst)
                print(f"Copied {plugin_path.name} to {dst}")
        else:
            print(f"Warning: Could not find built library {lib_name}")
            print(f"Searched in: {build_temp / 'lib'}, {build_temp}")
            # List files to debug
            if build_temp.exists():
                print(f"Build temp contents: {list(build_temp.iterdir())}")
                lib_dir = build_temp / "lib"
                if lib_dir.exists():
                    print(f"Lib dir contents: {list(lib_dir.iterdir())}")


setup(
    name="ebcc",
    packages=["ebcc"],
    ext_modules=[CMakeExtension("ebcc.libh5z_ebcc", sourcedir=os.path.dirname(__file__))],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    python_requires=">=3.8",
)
