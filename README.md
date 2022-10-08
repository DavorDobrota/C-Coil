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

### Getting started
- The library can be used natively from C++ and via Python bindings, which
  can also be used from MATLAB
- Documentation is available, and it is generated with Doxygen for C++. 
  Python methods have the same signature, but their names are written in 
  snake_case. It can be found in folder docs.
- Examples for 3 Python and 3 MATLAB use cases are provided.
- Performance benchmarks are available in module Benchmark and precision tests
  are available in module Comparison.

## Installation
- bla
