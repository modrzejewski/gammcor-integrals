#!/usr/bin/env python3
import sys
import shutil
from pathlib import Path


def main():
    if len(sys.argv) < 6:
        print("Usage: post_build.py <lib_file> <test_exe> <build_root> <source_root> <dummy_out>")
        sys.exit(1)

    lib_file = Path(sys.argv[1]).resolve()
    test_exe = Path(sys.argv[2]).resolve()
    build_root = Path(sys.argv[3]).resolve()
    source_root = Path(sys.argv[4]).resolve()
    dummy_out = Path(sys.argv[5])

    lib_dir = source_root / 'lib'
    include_dir = source_root / 'include'

    lib_dir.mkdir(parents=True, exist_ok=True)
    include_dir.mkdir(parents=True, exist_ok=True)

    # 1. Copy static library to lib/cholesky.a
    lib_dest = lib_dir / 'cholesky.a'
    shutil.copy2(lib_file, lib_dest)
    print(f"Copied static library to {lib_dest}")

    # 2. Copy test executable to gammcor-integrals/test
    test_dest = source_root / 'test'
    shutil.copy2(test_exe, test_dest)
    print(f"Copied test executable to {test_dest}")

    # 3. Copy all Fortran module files (.mod / .MOD) to include/
    mod_count = 0
    for mod_path in build_root.rglob('*.[mM][oO][dD]'):
        if mod_path.is_file():
            dest_file = include_dir / mod_path.name
            shutil.copy2(mod_path, dest_file)
            # Also create lowercase name if different
            dest_lower = include_dir / mod_path.name.lower()
            if dest_lower != dest_file:
                shutil.copy2(mod_path, dest_lower)
            mod_count += 1

    print(f"Synchronized {mod_count} module files to {include_dir}")

    # Touch dummy output for Meson custom_target
    dummy_out.write_text("done")


if __name__ == "__main__":
    main()
