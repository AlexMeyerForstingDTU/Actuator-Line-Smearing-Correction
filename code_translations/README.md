# Translations of `al_nw_smearing_correction.f`

The fixed-form Fortran smearing correction v3 ([../al_nw_smearing_correction.f](../al_nw_smearing_correction.f)) translated to Python, MATLAB, free-form Fortran and C++. Each translation keeps the algorithm, the routine structure and the order of the floating-point operations. The tests show that the Fortran and C++ translations reproduce the original bit for bit, and Python and MATLAB to within 3e-15.

The code is written to be read next to the papers: plain loops and array operations, one routine per model step, and the equation numbers in the comments. It has no error handling beyond that of the original.

## Files

| File | Content |
|---|---|
| [python/al_nw_smearing_correction.py](python/al_nw_smearing_correction.py) | Python 3 with numpy; the inner loops over the sections are numpy slices |
| [matlab/](matlab/) | MATLAB, one function per file |
| [fortran/al_nw_smearing_correction.f90](fortran/al_nw_smearing_correction.f90) | Free-form Fortran module |
| [cpp/al_nw_smearing_correction.hpp](cpp/al_nw_smearing_correction.hpp), [cpp/al_nw_smearing_correction.cpp](cpp/al_nw_smearing_correction.cpp) | C++17, standard library only |
| [tests/](tests/) | Test drivers for all five implementations, test runner and checks |

## Using a translation

The steps are the same in every language:

1. Create the state for `n_rotors` x `n_blades` actuator lines with `n_sections` control points each. The fourth argument, `project`, only names the output: it is the prefix of the diagnostics files (see below). `''` writes no files.
2. Every time step, call the correction once for each rotor `i` and blade `j` (both counted from 1). It returns the corrected section forces (spanwise, tangential, normal).

```python
from al_nw_smearing_correction import init_nw_smearing_correction, nw_smearing_correction

state = init_nw_smearing_correction(n_rotors, n_blades, n_sections, "my_run")
Fstn = nw_smearing_correction(state, i, j, time, dt, density, cyl_cp, cyl_vtp, vstn_cfd,
                              vcyl_cfd, eps, omega, chord, twist, pitch, cl_2d, cd_2d,
                              alpha_2d, cone, fbend, sweep, dsweep, l1point)
```

```matlab
state = init_nw_smearing_correction(n_rotors, n_blades, n_sections, 'my_run');
[Fstn, state] = nw_smearing_correction(state, i, j, time, dt, density, cyl_cp, cyl_vtp, ...
    vstn_cfd, vcyl_cfd, eps, omega, chord, twist, pitch, cl_2d, cd_2d, alpha_2d, cone, ...
    fbend, sweep, dsweep, l1point);
```

```fortran
call init_nw_smearing_correction(n_rotors, n_blades, n_sections, 'my_run')
call nw_smearing_correction(i, j, time, dt, density, cyl_cp, cyl_vtp, vstn_cfd, vcyl_cfd, &
                            eps, omega, chord, twist, pitch, cl_2d, cd_2d, alpha_2d, &
                            cone, fbend, sweep, dsweep, l1point, Fstn)
```

```cpp
NWState state = init_nw_smearing_correction(n_rotors, n_blades, n_sections, "my_run");
nw_smearing_correction(state, i, j, time, dt, density, cyl_cp, cyl_vtp, vstn_cfd, vcyl_cfd,
                       eps, omega, chord, twist, pitch, cl_2d, cd_2d, alpha_2d,
                       cone, fbend, sweep, dsweep, l1point, Fstn);
```

Inputs, units and array shapes are documented at the top of each file. Positions are in the cylindrical rotor frame (r, theta, z), angles in degrees, `omega` in rad/s. In MATLAB, the updated state is the second output and must be passed to the next call.

**Diagnostics files.** With `project = 'my_run'`, the correction writes one text file per rotor blade, named after the rotor and blade numbers: `my_run.nw01.bl01`, `my_run.nw01.bl02`, ... Each call of `nw_smearing_correction` adds one line with:

- the time [s], the number of sub-iterations and the final residual [m/s];
- at each section: the angle of attack from the CFD inflow and the corrected angle of attack [deg];
- at each trailing point: the helix angle [rad];
- at each section: the correction velocities `vt` and `vn` [m/s];
- at each trailing point: the trailed vortex strength `dGamma` [m^2/s];
- the CPU time of the call [s].

The format is the same as the original's. The tests use these files to compare each translation with the original line by line.

**Interface differences to the original** (the algorithm is unchanged):

- **Host values are arguments.** The original reads the time step, air density and time from the host modules `ns_contr` (`timestep0(grlev)`, `density`, `time`).
- **One state record per blade.** The original stores all blades in module arrays indexed `(i, j, ...)`. The translations keep one record per blade: `state.blade[i-1][j-1]` (Python), `state.blade(i, j)` (MATLAB), `blade(i, j)` (Fortran module variable) and `state.blade[(i-1)*n_blades + j-1]` (C++).
- **On by default.** `kaya` starts `true`. The original starts `.false.`, and the host switches it on before calling `init_nw_smearing_correction`.
- **File prefix is an argument.** The original takes the project name for the diagnostics files from the host code (`ACLproject` in the module `TurbineData`) and always writes the files on processor 1. The translations take it as the argument `project` and write no files when it is `''`. In Python and MATLAB, close the files at the end of the run: `state.blade[i][j].fid.close()`, `fclose(state.blade(i, j).fid)`.
- **Polar layout.** Python, MATLAB and Fortran use `(len_2d, n_s)` as the original. C++ stores one table per section, `cl_2d[l][m]`.

The settings `nw_max_resid` (convergence limit, 1e-6 m/s) and `x_aoa` (0.5) and the per-blade `refresh` flag (the original `kaya_refresh`) work as in the original.

## Routine map

| `al_nw_smearing_correction.f` | Python | MATLAB | Fortran | C++ |
|---|---|---|---|---|
| `init_nw_smearing_correction` | `init_nw_smearing_correction` | `init_nw_smearing_correction` | `init_nw_smearing_correction` | `init_nw_smearing_correction` |
| `nw_smearing_correction` | `nw_smearing_correction` | `nw_smearing_correction` | `nw_smearing_correction` | `nw_smearing_correction` |
| `smearing_factor` | `smearing_factor` | `smearing_factor` | `smearing_factor` | `smearing_factor` |
| `nw_hr_components` | `nw_hr_components` | `nw_hr_components` | `nw_hr_components` | `nw_hr_components` |
| `nw_phi_components` | `nw_phi_components` | `nw_phi_components` | `nw_phi_components` | `nw_phi_components` |
| `nw_XY_components` | `nw_xy_components` | `nw_xy_components` | `nw_xy_components` | `nw_xy_components` |
| `nw_circulation` | `nw_circulation` | `nw_circulation` | `nw_circulation` | `nw_circulation` |
| `nw_calc_relax` | `nw_calc_relax` | `nw_calc_relax` | `nw_calc_relax` | `nw_calc_relax` |
| `nw_interp1d` | `nw_interp1d` | `interp1(..., 'linear', 'extrap')` | `nw_interp1d` | `nw_interp1d` |
| `nw2vstn`, `Vxyz_Vstn` and rotations | `nw2vstn`, `vxyz_vstn`, `rot_x/y/z` | `nw2vstn`, `vxyz_vstn` (with `rot_x/y/z`) | `nw2vstn`, `vxyz_vstn`, `rot_x/y/z` | `nw2vstn`, `vxyz_vstn`, `rot_x/y/z` |

The interpolation extrapolates linearly beyond both ends of the table, as the original does. The three rotations about the x, y and z axes replace the original's rotation matrices and give the same numbers.

## Tests

```
py translations/tests/run_tests.py             # all cases, about 5 minutes
py translations/tests/run_tests.py eps10 units # selected cases
```

The runner needs `gfortran`, `g++`, MATLAB (`matlab` on the `PATH`) and a Python with numpy (the one that runs the script). It:

1. Builds three drivers: the original v3 with stubs for its host modules ([tests/host_stubs.f90](tests/host_stubs.f90)), the Fortran translation and the C++ translation. On Windows it links them with `-static`, so that another MinGW's runtime DLLs on the `PATH` cannot crash them.
2. Runs every case in every implementation, each in its own folder under `tests/build/runs/`. Each run writes `results.txt` and the diagnostics files. MATLAB runs all cases in one session, since it is slow to start.
3. Calls [tests/check_results.py](tests/check_results.py), which performs two checks:
   - **Translation:** every translation reproduces the original v3. `results.txt` must agree to `rtol = 1e-9`. Every line of the diagnostics files must agree to `rtol = 1e-6` (they print 8 digits), and the sub-iteration counts must agree exactly.
   - **Implementation:** all five implementations reproduce the numbers of the review note and exact properties of the method.

The test blade and inflow are those of the note's standalone test (Sect. 5): an NREL 5MW-like blade with 64 trailing points, 8 m/s at 9.2 rpm, and a frozen "CFD" inflow with induction 0.3. [tests/test_common.f90](tests/test_common.f90) defines them; the C++, Python and MATLAB drivers repeat the same definitions.

| Case | Setup | What it tests | Expected (source) |
|---|---|---|---|
| `units` | Single routines | Smearing factor, near-wake constants, root factor, Phi*, indicial functions, interpolation, rotations, circulation filter, relaxation factor | D_Y < 0 (note 3.1); C = 1.000, 0.992, 0.89, 0.69, 0.56 at h/r = 0.1 to 0.8 (note 3.2); `nw2vstn` gives -0.999 on every call (note 6.1); filter settles at Gamma_QS (note 4.2); relaxation factor 0.20 at 0.01 s (note 3.5) |
| `eps10_dt001` | eps = 10 m, dt = 0.01 s, 3000 steps | Main case of the note | P = 1.3307 MW, T = 287.82 kN, sub-iterations 3.2 / 9 (note 6.9) |
| `eps10` | eps = 10 m, dt = 0.1 s, 300 steps | Time-step independence | Same P and T as `eps10_dt001` (note 5) |
| `cone25` | 2.5 deg cone | 3D kinematics, `nw2vstn` | P = 1.3281 MW, T = 287.58 kN (note 6.9) |
| `eps25` | eps = 2.5 m | Smaller loop ranges | P between `eps10` and `nocorr` |
| `nocorr` | eps = 0.01 m | Correction vanishes | P = 1.4837 MW, T = 298.24 kN (note 5) |
| `bent` | Pre-bend, sweep, out-of-plane offsets, cone, Ang Li correction | All v2.0 geometry inputs | Converges; otherwise translation check only |
| `pitchstep` | Pitch +2 deg after 300 steps | Unsteady response: filter, histories, relaxation | Power drops |
| `noconv` | Convergence limit 0 | 500-iteration fallback | Correction stays zero, as the note describes (note 4.5) |
| `multi` | 2 rotors x 2 blades, one blade pitched | Per-blade state | Identical blades give identical loads |
| `rootzero` | First trailing point at r = 0 | Built-in switch-off | Correction switches itself off |

Result of the run on 13 October 2026 (gfortran and g++ 13.2, Python 3.14 with numpy 2.4.6, MATLAB R2024b): the Fortran and C++ translations are identical to the original in all 11 cases; Python and MATLAB agree to within 2.5e-15, with the same sub-iteration counts; all 21 implementation checks pass in all five implementations. The full run takes about 5 minutes, most of it in Python and MATLAB.
