import os
from glob import glob

from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup

from cuda_build import locate_CUDA, CudaBuildExt

# Load environment variables

use_GPU = int(os.environ.get("USE_GPU", 0))
GPU_increments = int(os.environ.get("GPU_INCREMENTS", 100))
float_type = str(os.environ.get("TYPE", "float"))


# Set compilation flags

if os.name == "nt":
    os.environ["CL"] = "/std:c++17 /O2 /arch:AVX2 /fp:fast"
else:
    os.environ["CPPFLAGS"] = "-std=c++17 -Ofast -mavx2 -ffast-math -flto=auto"

if use_GPU:
    cuda_compile_args = ["--use_fast_math", "-O3", "-std=c++17"]


# Set macros

macros = []

if GPU_increments < 1 or GPU_increments > 100:
    raise ValueError("GPU_INCREMENTS must be set to an integer in the range [1, 100]!")

if os.name == "nt":
    os.environ["CL"] += f" /DUSE_GPU#{use_GPU} /DGPU_INCREMENTS#{GPU_increments} /DTYPE#{float_type}"
else:
    macros += [("USE_GPU", use_GPU), ("GPU_INCREMENTS", GPU_increments), ("TYPE", float_type)]

if use_GPU:
    cuda_macros = [("GPU_INCREMENTS", GPU_increments), ("TYPE", float_type)]


# Define source paths

header_include_dirs = [
    "extern/CTPL/",
    "extern/pybind11/include/",
    "extern/HardwareAcceleration/",
    "src/",
    "src/Benchmark/",
    "src/Coil/",
    "src/Coil/CUDAKernels",
    "src/CoilGroup/",
    "src/CoilGroup/CUDAKernels",
    "src/Compare/",
    "src/CUDAUtils",
    "src/LegendreMatrix/",
    "src/Math/",
    "src/Tensor/",
    "src/Test/",
    "src/ThreadPool/",
    "src/Utils/",
]

sources = [
    *glob("src/*.cxx"),

    *glob("src/Benchmark/*.cxx"),
    *glob("src/Benchmark/Coil/Fields/*.cxx"),
    *glob("src/Benchmark/Coil/MInductanceAndForce/*.cxx"),
    *glob("src/Benchmark/CoilGroup/*.cxx"),
    *glob("src/Benchmark/Other/*.cxx"),

    *glob("src/Compare/*.cxx"),
    *glob("src/Compare/MInductanceAndForce/*.cxx"),
    *glob("src/Compare/Other/*.cxx"),
    *glob("src/Test/*.cxx"),
    *glob("src/Test/MInductanceAndForce/*.cxx"),
    *glob("src/Test/Other/*.cxx"),

    *glob("src/CoilGroup/*.cxx"),
    *glob("src/CoilGroup/Fields/*.cxx"),
    *glob("src/CoilGroup/ForceAndTorque/*.cxx"),
    *glob("src/CoilGroup/MInductance/*.cxx"),
    *glob("src/CoilGroup/Utils/*.cxx"),

    *glob("src/Coil/*.cxx"),
    *glob("src/Coil/Fields/*.cxx"),
    *glob("src/Coil/ForceAndTorque/*.cxx"),
    *glob("src/Coil/MInductance/*.cxx"),
    *glob("src/Coil/PrecisionArguments/*.cxx"),
    *glob("src/Coil/Utils/*.cxx"),

    *glob("src/Tensor/*.cxx"),
    *glob("src/Tensor/Matrix/*.cxx"),
    *glob("src/Tensor/Vector/*.cxx"),

    *glob("src/LegendreMatrix/*.cxx"),
    *glob("src/ThreadPool/*.cxx"),
    *glob("src/Utils/*.cxx"),
]


# Set up GPU building

extra_kwargs = {}

if use_GPU:
    CUDA = locate_CUDA()
    print(f"Using CUDA installation: {CUDA}")
    sources = sources + [
        *glob("src/Coil/CUDAKernels/*.cu"),
        *glob("src/CoilGroup/CUDAKernels/*.cu"),
        *glob("src/CUDAUtils/ErrorCheck/*.cu"),
        *glob("src/CUDAUtils/MemoryManagement/*.cu"),
    ]
    header_include_dirs += CUDA['include']

    library_dirs = []
    library_dirs += [CUDA['lib'] + '/x64'] if os.name == "nt" else []
    library_dirs += [CUDA['lib64']] if os.name != "nt" else []

    extra_kwargs = {
        'library_dirs': library_dirs,
        'libraries': ['cudart', 'cudadevrt'],
    }

    if os.name != "nt":
        extra_kwargs['runtime_library_dirs'] = library_dirs

    if os.name == "nt":
        cuda_compile_args += ["--compiler-options=/MD"]

    CudaBuildExt.setup(CUDA, cuda_compile_args, header_include_dirs, cuda_macros)
    cmdclass = {'build_ext': CudaBuildExt}
else:
    cmdclass = {}


# Finish setting up the building environment and run the build

ext_modules = [
    Pybind11Extension(
        name="c_coil",
        sources=sources,
        define_macros=macros,
        **extra_kwargs
    ),
]

with open("README.md", "r") as f:
    long_description = f.read().strip()

with open("LICENSE", "r") as f:
    project_license = f.read().strip()

with open("VERSION", "r") as f:
    project_version = f.read().strip()

setup(
    cmdclass=cmdclass,
    name="c-coil",
    version=project_version,
    ext_modules=ext_modules,
    author="Davor Dobrota, Nikola Socec, Lara Vrabac",
    author_email="nikola.socec@gmail.com",
    url="https://github.com/DavorDobrota/C-Coil",
    description="A library for performing calculations with coils.",
    long_description=long_description,
    zip_safe=False,
    python_requires=">=3.6",
    define_macros=macros,
    include_dirs=header_include_dirs,
    license=project_license,
)
