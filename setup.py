import os
from os.path import join
from glob import glob
from setuptools import setup
from distutils.command.build_ext import build_ext

from pybind11.setup_helpers import Pybind11Extension


use_GPU = int(os.environ.get("USE_GPU", 0))
GPU_increments = int(os.environ.get("GPU_INCREMENTS", 100))

if os.name == "nt":
    os.environ["CL"] = "/std:c++17 /O2 /arch:AVX2 /fp:fast"
else:
    os.environ["CPPFLAGS"] = "-std=c++17 -Ofast -mavx2 -ffast-math -flto=auto"

macros = []

if GPU_increments < 1 or GPU_increments > 100:
    raise ValueError("GPU_INCREMENTS must be set to an integer in the range [1, 100]!")

if os.name == "nt":
    os.environ["CL"] += f" /DUSE_GPU#{use_GPU} /DGPU_INCREMENTS#{GPU_increments}"
else:
    macros += [("USE_GPU", use_GPU), ("GPU_INCREMENTS", GPU_increments)]

if use_GPU:
    cuda_compile_args = ["--use_fast_math", "-O3"]

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
    *glob("src/Compare/*.cxx"),
    *glob("src/Test/*.cxx"),

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
    *glob("src/Math/*.cxx"),
    *glob("src/ThreadPool/*.cxx"),
    *glob("src/Utils/*.cxx"),
]


def find_in_path(names, path):
    for dir in path.split(os.pathsep):
        print(f"Looking for NVCC in: {dir}")
        for name in names:
            binpath = join(dir, name)
            if os.path.exists(binpath):
                return os.path.abspath(binpath)
    return None


def locate_CUDA():
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = join(home, 'bin', 'nvcc')
    else:
        nvcc = find_in_path(['nvcc', 'nvcc.exe'], os.environ['PATH'])
        if nvcc is None:
            raise EnvironmentError('The nvcc binary could not be located in your $PATH. Either add it to your path, or set $CUDAHOME')
        home = os.path.dirname(os.path.dirname(nvcc))
    cudaconfig = {'home': home, 'nvcc': nvcc, 'include': join(home, 'include'), 'lib64': join(home, 'lib64'), 'lib': join(home, 'lib')}
    for k, v in cudaconfig.items():
        if not os.path.exists(v):
            print('Warning: The CUDA %s path could not be located in %s' % (k, v))

    return cudaconfig


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
    extra_kwargs = {
        'library_dirs': [CUDA['lib64']],
        'libraries': ['cudart'],
        'runtime_library_dirs': [CUDA['lib64']],
        'extra_link_args': ['-lcudadevrt', '-lcudart'],
    }

ext_modules = [
    Pybind11Extension(
        name="coil_evolution",
        sources=sources,
        define_macros=macros,
        **extra_kwargs
    ),
]

with open("README.md", "r") as f:
    long_description = f.read()


def customize_compiler_for_nvcc(self):
    self.src_extensions.append('.cu')

    default_compiler_so = self.compiler_so
    super = self._compile

    nvcc_location = CUDA["nvcc"]

    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            self.set_executable('compiler_so', CUDA['nvcc'])
            extra_postargs = cuda_compile_args + ["--compiler-options=-fpic"]

        super(obj, src, ext, cc_args, extra_postargs, pp_opts)
        self.compiler_so = default_compiler_so

    self._compile = _compile


class cuda_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


cmdclass = {'build_ext': cuda_build_ext} if use_GPU else {}

setup(
    cmdclass=cmdclass,
    name="coil-evolution",
    version="1.0.0a0",
    ext_modules=ext_modules,
    author="Davor Dobrota, Nikola Socec",
    author_email="nikola.socec@gmail.com",
    url="https://github.com/DavorDobrota/Coil-Evolution",
    description="A library for performing calculations with coils.",
    long_description=long_description,
    zip_safe=False,
    python_requires=">=3.6",
    define_macros=macros,
    include_dirs=header_include_dirs,
)
