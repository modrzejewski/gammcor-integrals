# gammcor-cholesky
Authors: Marcin Modrzejewski, Michał Hapka
---

Molecular integrals library for the Gammcor code of Pernal et al.:
* Cholesky decomposition of Coulomb integrals,
* AO->MO transformation,
* interfaces for handling density matrices from external programs.

The build script of gammcor-cholesky generates static library `cholesky.a`, which should be linked to the main gammcor program.

## Installation
#### 1. Clone the repository
```
git clone git@github.com:modrzejewski/gammcor-cholesky.git
```

#### 2. Compile static library file
The compilation is done by python script `Build.py`. The user needs to define the `CompilerFlags` argument, which is specifies the compiler's command line. The value of `CompilerFlags` is the name of a subdirectory in `./src/CompilerFlags` where the command line for the compiler and linker are defined. For example, the subdirectory `./src/CompilerFlags/ifort-gammcor` stores two text files: `linker` and `compiler`. Those are raw text files which specify the linker and compiler command lines. The directory `./src/CompilerFlags` includes some predefined compiler options, but the user can freely add new subdirectories with customized parameters. The compilation should be done with parallelization to reduce the build time. The number of concurrent processes used during compilation is `-np`. Example build command using ifort and four concurrent processes:

```
cd <repository_name>
./Build.py -np 4 ifort-gammcor
```

The build script generates static library file `./lib/cholesky.a`. This file should be linked with the main gammcor program. Fortran header files (`.mod`) are stored in `./include`. The mod files directory should be included in gammcor compilation in order to import the subroutine interfaces from gammcor-cholesky.
