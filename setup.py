from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension


macros = [
    ("USE_GPU", 0)
]

header_include_dirs = [
    "extern/CTPL/",
    "extern/pybind11/include",
    "extern/HardwareAcceleration",
    "src/",
    "src/Coil/",
    "src/CoilData/",
    "src/CoilGroup/",
    "src/LegendreMatrix/",
    "src/Math/",
    "src/Tensor/",
    "src/Test/",
    "src/ThreadPool/",
]

sources = sorted([
    *glob("src/*.cxx"),
    *glob("src/Coil/*.cxx"),
    *glob("src/CoilData/*.cxx"),
    *glob("src/CoilGroup/*.cxx"),
    *glob("src/LegendreMatrix/*.cxx"),
    *glob("src/Math/*.cxx"),
    *glob("src/Tensor/*.cxx"),
    *glob("src/Test/*.cxx"),
    *glob("src/ThreadPool/*.cxx"),
])

ext_modules = [
    Pybind11Extension(
        "coil_evolution",
        sources,
    ),
]

with open("README.md", "r") as f:
    long_description=f.read()

setup(
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
