# C-Coil (Circular Coil Object-oriented Interaction Library)

## Overview
- This library is used for fast and precise calculation
  of magnetostatic fields, as well as interactions between circular
  coils in the form of mutual inductance, force, and torque.
- The model is satisfactory when the coil is made of a homogeneous 
  material, with no radial current, or when the field near the coil 
  is not of great interest, in the case of a compactly wound solenoid 
  with a large number of turns. 

### Coil Definition
- The main class Coil is characterized by 4 immutable attributes: inner 
  radius (R), thickness (a), length (b), and number of turns of 
  wire (N). Additional parameters are optional.
- Depending on the relationship of thickness and length, 4
  types of coils are distinguished with Enum CoilType:
  - Loop of wire (FILAMENT) - negligible length and thickness (a = 0.0, b = 0.0)
  - Thin solenoid (THIN) - negligible thickness (a = 0.0)
  - Pancake coil (FLAT) - negligible length (b = 0.0)
  - Circular coil with a rectangular cross-section (RECTANGULAR)
- Calculations with multi-coil systems, implemented in CoilGroup, are
  especially efficient when using the GPU

### Getting Started
- The library can be used natively from C++ and via Python bindings, which
  can also be used from MATLAB
- Documentation is available, and it is generated with Doxygen for C++. 
  Python methods have the same signature, but their names are written in 
  snake_case. It can be found in folder docs.
- Examples for 3 Python and 3 MATLAB use cases are provided.
- Performance benchmarks are available in module Benchmark and precision tests
  are available in module Comparison.

## Installation
The library can be built on Windows and Linux platforms, 
either as a C++ library (using CMake) or as a Python module (using setuptools and pybind11).

### Choosing build options
You can choose to build the project with different options:
- USE_GPU (bool/int) - enables or disables CUDA modules for GPU hardware acceleration
- GPU_INCREMENTS (int) 
  - determines the amount of increments used in GPU calculations
  - should be a value between 1 and 80, bounds included
  - increasing the increment count increases the precision of GPU calculations (it also increases compute time and VRAM usage), and vice versa
- TYPE (str)
  - determines the floating point type used in GPU calculations
  - should be either 'float' or 'double'
  - choosing 'double' will significantly decrease the calculation speed on most modern GPUs 

### Common requirements
#### CMake
- CMake version 3.8 or higher

#### Python 3
- Python version 3.8 or higher

#### CUDA
- CUDA toolkit version 11 or higher

#### Project dependencies
Make sure that all the project's dependencies (git submodules located in the 'extern' directory)
are cloned properly. This can be done by cloning the project with the command:  
`git clone https://github.com/DavorDobrota/C-COIL.git --recursive`  
If you have already cloned the repository, you can clone all the submodules by running:  
`git submodule update --init`

### Linux requirements
#### GCC
- GCC version 10 or higher

### Windows requirements
#### MSVC
- MSVC version 1920 or newer (Visual Studio 2019)

### Building with CMake

#### Running CMake
To build the project with CMake, run the following commands:
```shell
mkdir build
cd build
cmake ../
cmake --build . --config release
```

#### Choosing build options
To choose build options while building with CMake, simply append them to the `cmake ../` call,  
e.g. `cmake ../ -DUSE_GPU=1 -DGPU_INCREMENTS=80`

### Building the Python Module

#### Running the setup.py build
To build the Python module, run the following commands 
(using a [virtual environment](https://docs.python.org/3/library/venv.html#:~:text=A%20virtual%20environment%20is%20a,part%20of%20your%20operating%20system.) is recommended):
```shell
pip install -U -r requirements.txt
pip install .
```

#### Choosing build options
To choose build options while build the Python module on Linux, you need to pass them to
the build script as environment variables,  
e.g. `USE_GPU=1 GPU_INCREMENTS=80 pip install .`

The process is very similar on Windows and this example will work in Powershell:
```shell
$env:USE_GPU = 1
$env:GPU_INCREMENTS = 80
pip install .
```

#### Using the Python module in MATLAB 
After building the Python module, no further steps should be required to use it from MATLAB.
Make sure you open MATLAB from the virtual environment if you are using it.
See the MATLAB examples for extra information on using the Python module within MATLAB.


#### Important notes
1. On Windows, to use the CUDA modules with the Python code, 
   the CUDA binaries need to be added in Python as a DLL directory.
   See the Python examples for further reference.
   
2. When testing the Python module within MATLAB, we encountered an error caused by a
   mismatch of GCC's native libstdc++ version and MATLAB's bundeled libstdc++ version.
   We fixed the problem by overwriting MATLAB's version with the native GCC's version.
   This procedure comes with a risk, however, and we therefore advise caution if you are trying
   to replicate the hack on your machine.
