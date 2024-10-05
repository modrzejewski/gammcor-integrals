# gammcor-integrals
Authors: Marcin Modrzejewski, Michał Hapka
---

Molecular integrals library for the Gammcor code of Pernal et al.:
* Cholesky decomposition of Coulomb integrals,
* AO->MO transformation,
* interfaces for converting AO indices between Gammcor and external programs.

The build script of gammcor-cholesky generates static library `cholesky.a`, which should be linked to Gammcor's main binary file.

## Installation
#### 1. Clone the repository
```
git clone git@github.com:modrzejewski/gammcor-integrals.git
```

#### 2. Compile static library file
The compilation is done by running `Build.py` followed by two arguments:
* `-np` The number of concurrent processes used during compilation. Parallelization saves lots of build time if automatically generated ERI subroutines are compiled.
* `CompilerFlags` Name of the compiler and linker commands set. There are ready to use predefined sets in `./src/CompilerFlags/`. Just pick one of the subdirectory names, e.g., `ifort-gammcor`. User can create new subdirectories with custom params.

#### 3. Example: build using Intel and four concurrent compilation processes

```
cd <repository_name>
./Build.py -np 4 ifort-gammcor
```

The build script generates static library file `./lib/cholesky.a`. This file should be linked with the main gammcor program. Fortran header files (`.mod`) are stored in `./include`. The mod files directory should be added to the include path of Gammcor.
