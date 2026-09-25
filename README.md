# gammcor-integrals
Molecular integrals library for the Gammcor code of Pernal et al.:
* Cholesky decomposition of Coulomb integrals,
* AO->MO transformation,
* interfaces for converting AO indices between Gammcor and external programs.

The build system generates static library `lib/cholesky.a` and Fortran header files (`.mod`) in `./include`, which should be linked and included in Gammcor's build.

## Installation

### 1. Clone the repository
```bash
git clone git@github.com:modrzejewski/gammcor-integrals.git
cd gammcor-integrals
```

### 2. Setup the build directory with a chosen compiler profile
The compiler flags are configured using profile files located in `.meson/profiles/`. For example:
```bash
meson setup build --native-file .meson/profiles/ifx.ini
```

Profiles: `ifx.ini` (release) and `ifx-debug.ini`. The CPU target is detected at setup (`-Dcpu=auto`); it can be set with `-Dcpu=avx512`, `avx512-amd`, `avx2`, `native` or `dispatch`.

### 3. Compile the project
```bash
cd build
meson compile -j 4
```
Using a parallel build (e.g., `-j 4`) is highly recommended. It significantly speeds up compilation of the automatically generated ERI subroutines.

The build process generates:
* Static library file: `./lib/cholesky.a`
* Fortran module files (`.mod`): `./include/`
* Test binary: `./test`


---

## Contributors
* Marcin Modrzejewski
* Michał Hapka
* Aleksandra Tucholska
