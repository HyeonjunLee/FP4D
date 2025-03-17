# Compilation Instructions

This document provides step-by-step instructions for compiling the project.

## Prerequisites
Ensure you have the following dependencies installed:
- **HDF5** (compiled with the same compiler as the project)
- **NetCDF4** (compiled with the same compiler as the project)
- **Make** (build system)
- [Optional Tools] (debuggers, profilers, or additional utilities)

## Supported Build Environments
The project currently supports compilation using `Makefile` but does **not** provide CMake support.
By default, the following environments are provided:
- **HNUC**
- **KAIROS**

For other environments, users must create their own `Makefile`.

## Recommended Compiler
- **Intel Compiler (recommended)**
- **GNU Compiler (supported but requires HDF5 and NetCDF4 to be compiled with GNU as well)**

## Compilation Steps

### Using Makefile
To compile the project, navigate to the top-level directory and run:
```sh
make
```
For a clean build, use:
```sh
make clean && make
```
To install the compiled binaries:
```sh
make install
```

## Running the Executable
After compilation, run the program with:
```sh
./output_file [arguments]
```
For Windows:
```sh
output_file.exe [arguments]
```

## Additional Notes
- If you encounter errors, check the dependencies and paths.
- Ensure that HDF5 and NetCDF4 are compiled with the same compiler as the project.
- For optimized builds with Intel Compiler, additional flags may be required.
- For debugging, enable verbose output:
  ```sh
  make VERBOSE=1
  ```

For further details, refer to the project documentation or contact the maintainer.



///eqmodule/// courtesy to J. Song

  1) BlaBlaBla

///Miscellaneous///
  1) "blacs_mod.F90" in src

    Not currently used.

  2) "interpol.f" in src

    Not currently used.

  3) "blacs_mod.F90"

    Not currently uesd.
