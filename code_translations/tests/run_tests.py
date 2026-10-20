"""Build all implementations, run every test case, then check the results.

Usage:  py translations/tests/run_tests.py [case ...]
Needs gfortran, g++, MATLAB (matlab on the PATH) and a Python with numpy (the one running
this script). Outputs go to translations/tests/build/.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRANSLATIONS = HERE.parent
REPO = TRANSLATIONS.parent
BUILD = HERE / "build"
RUNS = BUILD / "runs"

CASES = ["units", "eps10_dt001", "eps10", "cone25", "eps25", "nocorr", "bent", "pitchstep",
         "noconv", "multi", "rootzero"]
IMPLS = ["original", "fortran", "cpp", "python", "matlab"]

# Windows: link the runtime libraries statically, so that another MinGW's DLLs on the PATH
# cannot be picked up at run time
STATIC = ["-static"] if os.name == "nt" else []


def build():
    """Compile the original (with host stubs), the Fortran and the C++ drivers."""
    (BUILD / "mod_original").mkdir(parents=True, exist_ok=True)
    (BUILD / "mod_fortran").mkdir(parents=True, exist_ok=True)
    original = ["gfortran", "-O2", *STATIC, "-J", BUILD / "mod_original",
                "-o", BUILD / "driver_original.exe",
                HERE / "host_stubs.f90", REPO / "al_nw_smearing_correction.f",
                HERE / "test_common.f90", HERE / "driver_original.f90"]
    fortran = ["gfortran", "-O2", *STATIC, "-J", BUILD / "mod_fortran",
               "-o", BUILD / "driver_fortran.exe",
               TRANSLATIONS / "fortran" / "al_nw_smearing_correction.f90",
               HERE / "test_common.f90", HERE / "driver_fortran.f90"]
    cpp = ["g++", "-O2", "-std=c++17", *STATIC, "-I", TRANSLATIONS / "cpp",
           "-o", BUILD / "driver_cpp.exe",
           TRANSLATIONS / "cpp" / "al_nw_smearing_correction.cpp", HERE / "driver_cpp.cpp"]
    for name, cmd in [("original", original), ("fortran", fortran), ("cpp", cpp)]:
        print(f"Building {name} ...")
        subprocess.run([str(c) for c in cmd], check=True)


def run(impl, cases):
    """Run all cases of one implementation, each in its own folder; return the wall time [s]."""
    for case in cases:
        work = RUNS / impl / case
        work.mkdir(parents=True, exist_ok=True)
        for old in work.iterdir():
            old.unlink()
    start = time.perf_counter()
    if impl == "matlab":
        # One MATLAB session for all cases, since MATLAB takes long to start
        script = f"addpath('{TRANSLATIONS / 'matlab'}'); addpath('{HERE}');"
        for case in cases:
            script = script + f" cd('{RUNS / 'matlab' / case}'); driver_matlab('{case}');"
        with open(RUNS / "matlab" / "log.txt", "w") as log:
            subprocess.run(["matlab", "-batch", script], stdout=log, stderr=subprocess.STDOUT,
                           check=True)
    else:
        for case in cases:
            if impl == "python":
                cmd = [sys.executable, str(HERE / "driver_python.py"), case]
            else:
                cmd = [str(BUILD / f"driver_{impl}.exe"), case]
            with open(RUNS / impl / case / "log.txt", "w") as log:
                subprocess.run(cmd, cwd=RUNS / impl / case, stdout=log, stderr=subprocess.STDOUT,
                               check=True)
    return time.perf_counter() - start


def main():
    cases = sys.argv[1:] or CASES
    build()
    print(f"\nRunning {len(cases)} case(s), wall time per implementation:")
    for impl in IMPLS:
        print(f"  {impl:10s} {run(impl, cases):6.1f} s")
    sys.path.insert(0, str(HERE))
    import check_results
    return check_results.main(cases)


if __name__ == "__main__":
    sys.exit(main())
