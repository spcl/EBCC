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
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        
        # Required for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = 'Debug' if debug else 'Release'

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]
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
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"
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

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        # Check if git submodules are initialized
        src_dir = Path(ext.sourcedir) / "src"
        openjpeg_dir = src_dir / "openjpeg"
        zstd_dir = src_dir / "zstd"
        
        if not (openjpeg_dir / "CMakeLists.txt").exists() or not (zstd_dir / "CMakeLists.txt").exists():
            print("Initializing git submodules...")
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

        src_lib_dir = Path(ext.sourcedir) / "src" / "build" 
        src_lib_dir.mkdir(parents=True, exist_ok=True)
        
        # Run CMake install to copy the library to the right location
        subprocess.run(
            ["cmake", "--install", ".", f"--prefix={src_lib_dir}"],
            cwd=build_temp,
            check=True,
        )
        
        # Find and copy the built library to where setuptools expects it
        lib_name = "libh5z_ebcc.so"
        if sys.platform == "win32":
            lib_name = "h5z_ebcc.dll"
        elif sys.platform == "darwin":
            lib_name = "libh5z_ebcc.dylib"
        
        # Look for the library in the install location
        installed_lib = src_lib_dir / "lib" / lib_name
        if not installed_lib.exists():
            # Try the build directory
            installed_lib = build_temp / "lib" / lib_name
            if not installed_lib.exists():
                installed_lib = build_temp / lib_name
        
        # Copy to the location setuptools expects (where the .so should be)
        if installed_lib.exists():
            import shutil
            target_lib = Path(self.get_ext_fullpath(ext.name))
            target_lib.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(installed_lib, target_lib)
            print(f"Copied {lib_name} to {target_lib}")
        else:
            print(f"Warning: Could not find built library {lib_name}")
            print(f"Searched in: {src_lib_dir / 'lib'}, {build_temp / 'lib'}, {build_temp}")
            # List files to debug
            if build_temp.exists():
                print(f"Build temp contents: {list(build_temp.iterdir())}")
                lib_dir = build_temp / "lib"
                if lib_dir.exists():
                    print(f"Lib dir contents: {list(lib_dir.iterdir())}")


setup(
    name="ebcc",
    ext_modules=[CMakeExtension("libh5z_ebcc", sourcedir=".")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    python_requires=">=3.8",
)