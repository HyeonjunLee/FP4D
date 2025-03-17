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

## Environment Variables

Before compiling, set the following environment variables:

```sh
export FP4D_PLATFORM=HNUC   # Set this to the target environment (e.g., HNUC, KAIROS)
export FP4D_ROOT=/path/to/project   # Set this to the root directory of the project
```

## Creating a Custom Makefile

If you need to create a custom `Makefile` for a new environment, place it in the following directory:

```sh
${FP4D_ROOT}/platforms
```

The `Makefile` should be named as follows:

```sh
makefile.${FP4D_PLATFORMS}
```

where `${FP4D_PLATFORMS}` corresponds to the name of the target environment (e.g., `HNUC`, `KAIROS`).

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

## Running the Executable

After compilation, run the program with:

```sh
./FP4D [arguments]
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


