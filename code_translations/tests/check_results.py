"""Check the test results written by run_tests.py.

1. Translation: each translation (fortran, cpp, python, matlab) reproduces the original v3
   code, both results.txt and every line of the diagnostics files test.nwRR.blBB
   (sub-iteration counts must match exactly).
2. Implementation: every implementation reproduces the numbers in
   Note_review_of_PhD_comments.md ("note") and simple exact properties.

Usage: py translations/tests/check_results.py [case ...]
"""
import math
import sys
from pathlib import Path

import numpy as np

RUNS = Path(__file__).resolve().parent / "build" / "runs"
REFERENCE = "original"
TRANSLATIONS = ["fortran", "cpp", "python", "matlab"]

# Translation tolerances: results.txt (17 digits) and diagnostics files (8 digits)
RTOL, ATOL = 1e-9, 1e-12
RTOL_DIAG, ATOL_DIAG = 1e-6, 1e-12

# h/r of the unit test of nw_hr_components, and the fit constants P(1:4, 1:4) of the
# mixing factor for h/r > 0 (Pirrung et al., 2017b, Eq. 14)
HR_LIST = [0.1, 0.3, 0.5, 0.7, 0.8, -0.5, -3.0]
FIT_P = [[-1.64636754184988, 8.14821474772595, -12.1784861715161, 5.02652654794305],
         [-0.499006672527007, 6.08465268675599, -15.1712005446736, 14.8254135311781],
         [3.90835606394205, -18.7662330478906, 39.1243282362892, -29.4870086106427],
         [-1.60622843516794, 7.42952711674470, -15.8594777074799, 11.6870163960092]]


def read_results(folder):
    """results.txt as a dictionary: name -> array of values."""
    res = {}
    for line in (folder / "results.txt").read_text().splitlines():
        name, *numbers = line.split()
        res[name] = np.array([float(x) for x in numbers])
    return res


def read_diagnostics(folder):
    """Diagnostics files as a dictionary: file name -> (sub-iterations, values).

    Each line holds the time (15 characters), the sub-iterations (5), then 15 characters per
    value. The last value, the CPU time, is left out. Empty files are skipped: the original
    opens its files in the initialisation, also when the correction never writes to them.
    """
    diag = {}
    for path in sorted(folder.glob("test.nw*.bl*")):
        lines = path.read_text().splitlines()
        if len(lines) == 0:
            continue
        sub_iter = []
        values = []
        for line in lines:
            sub_iter.append(int(line[15:20]))
            values.append([float(line[0:15])]
                          + [float(line[k:k + 15]) for k in range(20, len(line) - 15, 15)])
        diag[path.name] = (np.array(sub_iter), np.array(values))
    return diag


def compare(res, diag, res_ref, diag_ref):
    """Compare one run with the original: (ok, identical, maximum relative difference)."""
    if res.keys() != res_ref.keys() or diag.keys() != diag_ref.keys():
        return False, False, math.nan
    ok = True
    pairs = []
    for name in res_ref:
        pairs.append((res[name], res_ref[name], RTOL, ATOL))
    for name in diag_ref:
        ok = ok and np.array_equal(diag[name][0], diag_ref[name][0])   # sub-iterations: exact
        pairs.append((diag[name][1], diag_ref[name][1], RTOL_DIAG, ATOL_DIAG))
    identical = ok
    worst = 0.0
    for a, b, rtol, atol in pairs:
        if a.shape != b.shape:
            return False, False, math.nan
        ok = ok and np.all(np.abs(a - b) <= atol + rtol*np.abs(b))
        identical = identical and np.array_equal(a, b)
        if b.size > 0:
            worst = max(worst, np.max(np.abs(a - b))/max(np.max(np.abs(b)), 1e-300))
    return ok, identical, worst


def wang_coton(h_r):
    """Phi of Wang and Coton (Pirrung et al., 2016, Eq. 4)."""
    if h_r > 0:
        return math.pi/4*abs((1 + 0.5*h_r)*math.log(1 - h_r))
    return math.log(1 - h_r)/(1.5 + math.log(1 - 0.5*h_r))


def root_factor(phi):
    """Root correction C for D_X = 1.359 and D_Y = -0.359/4 (2016 Eq. 27-29)."""
    return ((1.359*(1 - math.exp(-math.pi/(2*phi))) - 0.359/4*(1 - math.exp(-2*math.pi/phi)))
            / (1.359 - 0.359/4))


def k_fit(h_r, phi):
    """Mixing factor k_Phi for h/r > 0 before clipping (2017b Eq. 10 and 14)."""
    a = [p[0] + p[1]*h_r + p[2]*h_r**2 + p[3]*h_r**3 for p in FIT_P]
    return a[0] + a[1]*phi + a[2]*phi**2 + a[3]*phi**3


def rel(a, b):
    return abs(a - b)/abs(b)


def sub_iterations(diag):
    """Mean and maximum sub-iterations of rotor 1, blade 1, and the number of calls with 500."""
    s = diag["test.nw01.bl01"][0]
    return s.mean(), s.max(), np.sum(s == 500)


def implementation_checks(res, diag):
    """Checks of one implementation against the note and exact properties.

    res[case] and diag[case] hold the results and diagnostics of each case that was run.
    Returns a list of (description, passed, detail).
    """
    checks = []

    if "units" in res:
        u = res["units"]

        # Smearing factor
        f = u["smearing_factor"]
        expected = [math.exp(-36/100), math.exp(-(36 + (2*math.cos(0.3))**2)/100),
                    math.exp(-(25 + math.cos(0.5)**2)/6.25)]
        checks.append(("Gaussian smearing factor, with and without out-of-plane offset "
                       "[MF2019 Eq. 21-22]",
                       np.allclose(f, expected, rtol=1e-14, atol=0),
                       f"f = {f[0]:.6f}, {f[1]:.6f}, {f[2]:.6f}"))

        # Amplitudes D_X and D_Y, and Phi_s. Rows of v: phi, phi_s, dx, dy, a(1:4) per h/r
        v = u["hr_components"].reshape(7, 8)
        ok = True
        for h_r, row in zip(HR_LIST, v):
            h = 10*h_r
            ok = ok and rel(row[2]*4*math.pi*abs(h)*h_r, 1.359) < 1e-12
            ok = ok and rel(row[3]/row[2], -0.359/4/1.359) < 1e-12
            ok = ok and rel(row[1], 0.788*abs(h_r)) < 1e-14     # Phi_s, 2017b Eq. 7
        checks.append(("D_X = 1.359 r/(4 pi h|h|), D_Y = -(0.359/4) r/(4 pi h|h|), i.e. D_Y < 0 "
                       "[note Sect. 3.1]",
                       ok, f"D_Y/D_X = {v[0, 3]/v[0, 2]:.6f}"))

        # Root factor C = Phi_C/Phi
        C = {h_r: row[0]/wang_coton(h_r) for h_r, row in zip(HR_LIST, v)}
        note = {0.1: 1.000, 0.3: 0.992, 0.5: 0.89, 0.7: 0.69, 0.8: 0.56}
        ok = all(abs(C[h_r] - note[h_r]) <= 0.005 for h_r in note)
        ok = ok and all(abs(C[h_r] - root_factor(wang_coton(h_r))) < 1e-12 for h_r in HR_LIST)
        checks.append(("Root factor C = Phi_C/Phi applied once: 1.000, 0.992, 0.89, 0.69, 0.56 at "
                       "h/r = 0.1, 0.3, 0.5, 0.7, 0.8 [note Sect. 3.2]",
                       ok, "C = " + ", ".join(f"{C[h_r]:.4f}" for h_r in note)))

        # Phi*: first pair clipped to k = 0, second pair inside [0, 1]
        ps = u["phi_components"]
        k1 = k_fit(0.1, math.radians(5))
        k2 = k_fit(0.3, math.radians(30))
        ok = k1 < 0 and ps[0] == v[0, 0]
        ok = ok and 0 < k2 < 1 and rel(ps[1], k2*v[1, 1] + (1 - k2)*v[1, 0]) < 1e-12
        checks.append(("Phi* = k Phi_s + (1 - k) Phi_C, k clipped to [0, 1] "
                       "[2017b Eq. 8-10, corrigendum]",
                       ok, f"k = {k1:.3f} -> 0 at (0.1, 5 deg), k = {k2:.3f} at (0.3, 30 deg)"))

        # Indicial functions in steady state
        X, Y = u["xy_steady"]
        ok = rel(X, 2*0.0123*0.09) < 1e-12 and rel(Y, 2*-0.0008125*0.09) < 1e-12
        checks.append(("Indicial functions settle at dGamma D_X Phi* and dGamma D_Y Phi* "
                       "[2016 Eq. 25-26]",
                       ok, f"X = {X:.6e}, Y = {Y:.6e}"))

        # Interpolation
        y = u["interp1d"]
        checks.append(("Interpolation inside and linear extrapolation beyond both ends",
                       np.allclose(y, [-10, 0, 5, 15, 22.5, 30, 45], rtol=0, atol=1e-12),
                       "y = " + ", ".join(f"{value:g}" for value in y)))

        # nw2vstn called repeatedly
        v_rep = u["nw2vstn_repeat"]
        checks.append(("nw2vstn returns -0.999 on every call for a 2.5 deg cone "
                       "[note Sect. 4.1, 6.1]",
                       np.allclose(v_rep, -math.cos(math.radians(2.5)), rtol=1e-14, atol=0),
                       f"{v_rep[0]:.6f} ... {v_rep[-1]:.6f}"))

        # Rotations keep the vector length
        g = u["nw2vstn_general"].reshape(2, 3)
        rotated = [g[0], g[1], u["vxyz_vstn"]]
        original = [[0, 1.5, -0.7], [0, -0.3, 2.0], [1, -2, 0.5]]
        ok = True
        for a, b in zip(rotated, original):
            ok = ok and rel(np.linalg.norm(a), np.linalg.norm(b)) < 1e-14
        checks.append(("Rotations keep vector lengths",
                       ok, f"|v| = {np.linalg.norm(u['vxyz_vstn']):.12f}"))

        # Trailed vortex strengths
        G = u["circulation_gamma"]
        dG = u["circulation_dgamma"]
        ok = dG[0] == G[0] and dG[-1] == -G[-1] and abs(dG.sum()) < 1e-12*np.max(np.abs(G))
        checks.append(("Trailed strengths: dGamma(1) = Gamma(1), dGamma(n_v) = -Gamma(n_s), "
                       "sum = 0",
                       ok, f"max Gamma = {np.max(G):.4f} m^2/s"))

        # Circulation filter in steady state
        ratio = u["circulation_ratio"][0]
        checks.append(("Circulation filter settles at Gamma_QS (A_1 + A_2 + A_3 = 1) "
                       "[note Sect. 4.2]",
                       ratio < 1e-9, f"max |Gamma/Gamma_QS - 1| = {ratio:.1e}"))

        # Relaxation factor
        f = u["calc_relax"]
        ok = abs(f[0] - 0.20) <= 0.005 and f[1] > f[0] and all(0.2 <= x <= 0.99 for x in f)
        checks.append(("Relaxation factor 0.20 at dt = 0.01 s, larger at 0.1 s, within [0.2, 0.99] "
                       "[note Sect. 3.5]",
                       ok, f"factor = {f[0]:.3f}, {f[1]:.3f}"))

    if "eps10_dt001" in res:
        P = res["eps10_dt001"]["P"][0]/1e6
        T = res["eps10_dt001"]["T"][0]/1e3
        mean, mx, n500 = sub_iterations(diag["eps10_dt001"])
        ok = abs(P - 1.3307) <= 5e-5 and abs(T - 287.82) <= 5e-3
        ok = ok and abs(mean - 3.2) <= 0.05 and mx == 9 and n500 == 0
        checks.append(("eps = 10 m, dt = 0.01 s: P = 1.3307 MW, T = 287.82 kN, sub-iterations "
                       "3.2 / 9 [note Sect. 6.9]",
                       ok, f"P = {P:.5f} MW, T = {T:.3f} kN, sub-iterations {mean:.2f} / {mx}"))

    if "eps10" in res and "eps10_dt001" in res:
        a = res["eps10"]
        b = res["eps10_dt001"]
        mean, mx, n500 = sub_iterations(diag["eps10"])
        ok = rel(a["P"][0], b["P"][0]) < 1e-5 and rel(a["T"][0], b["T"][0]) < 1e-5 and n500 == 0
        checks.append(("Converged loads independent of the time step (0.1 s vs 0.01 s) "
                       "[note Sect. 5]",
                       ok, f"P differs by {rel(a['P'][0], b['P'][0]):.1e}, "
                           f"sub-iterations {mean:.2f} / {mx}"))

    if "cone25" in res:
        P = res["cone25"]["P"][0]/1e6
        T = res["cone25"]["T"][0]/1e3
        checks.append(("2.5 deg cone: P = 1.3281 MW, T = 287.58 kN [note Sect. 6.9]",
                       abs(P - 1.3281) <= 5e-5 and abs(T - 287.58) <= 5e-3,
                       f"P = {P:.5f} MW, T = {T:.3f} kN"))

    if "nocorr" in res:
        r = res["nocorr"]
        P = r["P"][0]/1e6
        T = r["T"][0]/1e3
        zero = not np.any(r["vn_history"]) and not np.any(r["vt_history"])
        checks.append(("No correction (eps = 0.01 m): P = 1.4837 MW, T = 298.24 kN, correction "
                       "zero [note Sect. 5]",
                       abs(P - 1.4837) <= 5e-5 and abs(T - 298.24) <= 5e-3 and zero,
                       f"P = {P:.5f} MW, T = {T:.3f} kN"))

    if "eps10" in res and "eps25" in res and "nocorr" in res:
        P = [res[case]["P"][0]/1e6 for case in ["eps10", "eps25", "nocorr"]]
        checks.append(("Smaller eps, smaller correction: P(10 m) < P(2.5 m) < P(none) "
                       "[note Sect. 5]",
                       P[0] < P[1] < P[2], f"P = {P[0]:.5f} < {P[1]:.5f} < {P[2]:.5f} MW"))

    if "bent" in res:
        mean, mx, n500 = sub_iterations(diag["bent"])
        P = res["bent"]["P"][0]/1e6
        checks.append(("Pre-bent, swept, coned blade with the Ang Li correction converges",
                       n500 == 0 and math.isfinite(P),
                       f"P = {P:.5f} MW, sub-iterations {mean:.2f} / {mx}"))

    if "pitchstep" in res:
        P = res["pitchstep"]["P_series"]
        T = res["pitchstep"]["T_series"]
        n500 = sub_iterations(diag["pitchstep"])[2]
        checks.append(("Pitching 2 deg to feather lowers power and thrust",
                       P[-1] < P[0] and T[-1] < T[0] and n500 == 0,
                       f"P: {P[0]/1e6:.5f} -> {P[-1]/1e6:.5f} MW"))

    if "noconv" in res and "nocorr" in res:
        r = res["noconv"]
        s = diag["noconv"]["test.nw01.bl01"][0]
        zero = not np.any(r["vn_history"]) and not np.any(r["vt_history"])
        same_P = r["P"][0] == res["nocorr"]["P"][0]
        checks.append(("No convergence: 500 sub-iterations, correction silently stays zero "
                       "[note Sect. 4.5]",
                       np.all(s == 500) and zero and same_P,
                       "sub-iterations " + ", ".join(str(n) for n in s)
                       + "; P equals the uncorrected P"))

    if "multi" in res:
        P = res["multi"]["P_blade"]
        T = res["multi"]["T_blade"]
        ok = P[0] == P[1] == P[2] and T[0] == T[1] == T[2] and P[3] != P[0] and T[3] != T[0]
        checks.append(("Each rotor blade keeps its own history (2 rotors x 2 blades)",
                       ok, f"P per blade = {P[0]/1e6:.5f} (x3), {P[3]/1e6:.5f} MW (pitch 1 deg)"))

    if "rootzero" in res:
        checks.append(("First trailing point at r = 0: the correction switches itself off",
                       res["rootzero"]["kaya_after"][0] == 0, "kaya = off"))

    return checks


def main(cases=None):
    cases = cases or sorted(p.name for p in (RUNS / REFERENCE).iterdir())
    impls = [REFERENCE] + TRANSLATIONS
    res = {}
    diag = {}
    for impl in impls:
        res[impl] = {case: read_results(RUNS / impl / case) for case in cases}
        diag[impl] = {case: read_diagnostics(RUNS / impl / case) for case in cases}
    failures = 0

    # 1. Translation: compare each translation with the original
    print("\nTranslation check: same results as the original v3 code")
    print(f"  {'case':12s}" + "".join(f"{impl:>20s}" for impl in TRANSLATIONS))
    for case in cases:
        line = f"  {case:12s}"
        for impl in TRANSLATIONS:
            ok, identical, worst = compare(res[impl][case], diag[impl][case],
                                           res[REFERENCE][case], diag[REFERENCE][case])
            if identical:
                text = "identical"
            else:
                text = f"max diff {worst:.1e}"
            if not ok:
                text = "FAIL " + text
                failures += 1
            line += f"{text:>20s}"
        print(line)

    # 2. Implementation: the same checks for every implementation
    print("\nImplementation check: all five implementations "
          "(note = Note_review_of_PhD_comments.md)")
    checks = {impl: implementation_checks(res[impl], diag[impl]) for impl in impls}
    for n, (description, _, detail) in enumerate(checks[REFERENCE]):
        failed = [impl for impl in impls if not checks[impl][n][1]]
        if failed:
            label = "FAIL " + ", ".join(failed)
            failures += 1
        else:
            label = "PASS"
        print(f"  {label:4s}  {description}")
        print(f"        {detail}")

    if failures == 0:
        print("\nAll checks passed.")
        return 0
    print(f"\n{failures} check(s) failed.")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or None))
