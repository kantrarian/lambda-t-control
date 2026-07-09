# -*- coding: utf-8 -*-
"""lzsm_sweeps.py -- Two-level TDSE simulation of LZSM sweeps for
"Predictive Geometric Control of Quantum Dynamics via Spectral Curvature
Functional Lambda(t)" (R. J. Mathews).

System (paper Eq. 2, hbar-normalized generator B(t) = H(t)/hbar):

    B(t) = (1/2) * (Delta(t) * sigma_z + Omega * sigma_x)      [rad/s]

Curvature functional (paper Eq. 3):  Lambda(t) = ||[B(t), Bdot(t)]||_op.
Lemma 1 (paper Eq. 4):               Lambda(t) = (1/2) * Omega * |Delta_dot(t)|.

Units: nondimensionalized with Omega_ref = 1 rad/s; all frequencies in
rad/s, times in seconds of that unit system. hbar-normalization B = H/hbar
is declared as a telemetry OperationRecord.

THEORY FORMS UNDER TEST (both simulated; algebra check done before coding):
  * Standard Landau-Zener for this Hamiltonian:
        P_std   = exp(-pi * Omega^2 / (2 * nu*)),   nu* = |Delta_dot(t*)|
    which via Lemma 1 (nu* = 2*Lambda*/Omega) reads
        P_std   = exp(-pi * Omega^3 / (4 * Lambda(t*))).
  * The paper's printed form (its Eq. 1 uses exp(-pi Omega^2/nu); its Eq. 5
    follows consistently from that):
        P_paper = exp(-pi * Omega^3 / (2 * Lambda(t*)))
             = exp(-pi * Omega^2 / nu*).
  These differ by a factor 2 in the exponent. The paper's Eq. 5 IS
  algebraically consistent with its own Eq. 1, but its Eq. 1 differs from
  the standard LZ result exp(-pi Omega^2/(2 nu)) [Zener 1932] by that
  factor 2. The TDSE below is the arbiter; both forms are overlaid and the
  discrepancy is reported, not papered over.

DETUNING SCHEDULES (three; the .tex names none, so the task's defaults),
all crossing Delta = 0 at t* = 0, sweeping -Delta0 -> +Delta0,
parametrized by the crossing slope nu* = |Delta_dot(0)|:
  linear:     Delta(t) = nu* t,                      t in [-Delta0/nu*, +Delta0/nu*]
  tanh:       Delta(t) = Delta0 tanh(t/tau),         tau = Delta0/nu*, t in [-6 tau, 6 tau]
  gauss_flat: Delta_dot(t) = A (1 - a exp(-t^2/(2 sigma^2))), a = 0.75,
              A = nu*/(1-a), sigma = 6 * max(1/sqrt(nu*), Omega/nu*),
              Delta(t) = A [t - a sigma sqrt(pi/2) erf(t/(sigma sqrt(2)))],
              domain solved so Delta(+-T/2) = +-Delta0.
              (Locally curvature-flattened: 4x faster than nu* away from
              the crossing, slope nu* at the crossing.)

SWEEP GRID: Omega in {0.5, 1.0, 2.0} rad/s, Delta0 = 30*Omega (fixed
ratio), 10 target exponents x = pi*Omega^2/(2 nu*) log-spaced in [0.2, 8]
=> P_std in [3.4e-4, 0.82]; 3 schedules x 3 Omegas x 10 rates = 90 runs.
MATCHED-T COMPARISON (the abstract's claim): Omega = 1, x_lin in
{0.5, 1, 2}; linear vs gauss_flat (a = 0.85, sigma = 0.1 T) at identical
sweep range AND identical total time T; measured leakage ratio reported.

INTEGRATOR: scipy solve_ivp DOP853, rtol = 1e-11, atol = 1e-13. Initial
state = adiabatic ground state of B(-T/2); leakage P_sim = final
population of the excited adiabatic state of B(+T/2).

============================ PRE-REGISTRATION ============================
Acceptance criteria fixed BEFORE the first full run (task discipline):

SC1 (GATED): Lemma-1 identity check. On a 2001-point time grid per run,
  Lambda_num(t) = ||[B(t), Bdot_fd(t)]||_2 with Bdot_fd a 4th-order
  central finite difference (step h = 1e-3 * schedule timescale) must
  match Lambda_analytic(t) = (1/2) Omega |Delta_dot(t)| with relative
  error < 1e-6 at every grid point excluding the outermost 2% of points
  at each end (endpoint exclusion pre-registered because tanh's crossing
  slope decays to ~2.4e-5 of peak there, amplifying rounding). PASS/FAIL.

SC2 (REPORTED, NOT GATED -- per task instruction): closed-form agreement.
  Report max |ln P_sim - ln P_theory| over all sweep points with
  P_sim > 1e-6, separately for P_std and P_paper, per schedule and
  pooled; plus the least-squares slope b of -ln P_sim vs
  pi*Omega^2/(2 nu*) (b ~ 1 favors the standard form, b ~ 2 the paper
  form). No pass bar is set here; numbers are reported as measured.
  Directional expectation stated for honesty: the standard form should
  fit; the paper form should show a systematic factor ~2.

SC3 (GATED): unitarity. | ||psi(T)||_2 - 1 | < 1e-8 for every run.

SC4 (GATED): every telemetry record emitted (capsules.jsonl,
  operations.jsonl, defects.jsonl) validates against the three schemas
  in ../schemas/ (jsonschema).

SC5 (REPORTED): solver convergence. Re-run 6 runs (each schedule at
  Omega in {0.5, 2.0}, slowest sweep x = 8) at rtol = 1e-9/atol = 1e-11;
  report max |Delta ln P_sim| against the tight-tolerance result.

Seed: 0 (numpy default_rng(0)); the pipeline is deterministic -- the seed
is recorded so any future stochastic extension inherits it.
==========================================================================

Outputs (paths relative to the build dir, i.e. this file's parent's parent):
  results/lzsm_sweeps.json           all sweep numbers, seeds, tolerances
  figures/fig_curvature_sweeps.{png,pdf}
  telemetry/{capsules,operations,defects}.jsonl

Telemetry tiers (minimal-sufficiency declaration, SKILL.md rule 4):
  Tier 3 for the representative runs (one per schedule, Omega = 1, mid
  ramp rate): B(t) capsules at 33 samples (units rad/s, Pauli frame,
  spectral signature with declared pairing sorted_by_real -- gap is
  bounded below by Omega at the avoided crossing, so the pairing is
  well-conditioned), Bdot capsules, commutator DefectRecords per sample,
  eigenvalue-drift DefectRecords between consecutive samples,
  OperationRecords for hbar-normalization, finite-difference derivative,
  TDSE integration, and the adiabatic-projection measurement (a quotient:
  the relative phase of the adiabatic amplitudes is destroyed and NOT
  side-channeled -- declared, with a claim_effects invalidation).
  Tier 0/1 for scalar sweep outputs (P values in results JSON).

Run:  python lzsm_sweeps.py          (from the sim/ directory)
"""

import json
import math
import os
import platform
import sys
import time
import uuid

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import erf

import st_emit

# --------------------------------------------------------------- paths
SIM_DIR = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.dirname(SIM_DIR)
RESULTS_DIR = os.path.join(BUILD_DIR, "results")
FIG_DIR = os.path.join(BUILD_DIR, "figures")
TEL_DIR = os.path.join(BUILD_DIR, "telemetry")
SCHEMA_DIR = os.path.join(BUILD_DIR, "schemas")

RUN_ID = "lzsm_sweeps_v1_seed0"
SEED = 0
RTOL, ATOL = 1e-11, 1e-13
RTOL_LOOSE, ATOL_LOOSE = 1e-9, 1e-11
DELTA0_OVER_OMEGA = 30.0
OMEGAS = (0.5, 1.0, 2.0)
X_GRID = np.geomspace(0.2, 8.0, 10)          # target standard-form exponents
GAUSS_A = 0.75
GAUSS_SIGMA_FACTOR = 6.0
TANH_HALF_WIDTHS = 6.0                        # domain +- 6 tau
N_LAMBDA_GRID = 2001
FD_H_FACTOR = 1e-3
SC1_EDGE_EXCLUDE = 0.02                       # outermost 2% of grid points
SC1_GATE = 1e-6
SC3_GATE = 1e-8
N_TELEMETRY_SAMPLES = 33

SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]])
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]])


# ----------------------------------------------------------- schedules

class Schedule:
    """Closed-form detuning schedule Delta(t) with analytic derivative.

    Delta and dDelta are defined on all of R (they extend smoothly past
    the integration domain), so central finite differences are valid at
    the domain endpoints too.
    """

    def __init__(self, name, Delta, dDelta, t_half, timescale, params):
        self.name = name
        self.Delta = Delta
        self.dDelta = dDelta
        self.t_half = float(t_half)       # integrate t in [-t_half, +t_half]
        self.timescale = float(timescale)  # characteristic time for FD step
        self.params = params


def make_linear(nu, Delta0):
    return Schedule("linear",
                    lambda t: nu * t,
                    lambda t: nu * np.ones_like(np.asarray(t, dtype=float)),
                    Delta0 / nu, 2.0 * Delta0 / nu,
                    {"nu_star": nu, "Delta0": Delta0})


def make_tanh(nu, Delta0):
    tau = Delta0 / nu
    return Schedule("tanh",
                    lambda t: Delta0 * np.tanh(t / tau),
                    lambda t: (Delta0 / tau) / np.cosh(t / tau) ** 2,
                    TANH_HALF_WIDTHS * tau, tau,
                    {"nu_star": nu, "Delta0": Delta0, "tau": tau})


def make_gauss_flat(nu, Delta0, Omega, a=GAUSS_A, sigma=None):
    if sigma is None:
        sigma = GAUSS_SIGMA_FACTOR * max(1.0 / math.sqrt(nu), Omega / nu)
    A = nu / (1.0 - a)

    def Delta(t):
        t = np.asarray(t, dtype=float)
        return A * (t - a * sigma * math.sqrt(math.pi / 2.0)
                    * erf(t / (sigma * math.sqrt(2.0))))

    def dDelta(t):
        t = np.asarray(t, dtype=float)
        return A * (1.0 - a * np.exp(-t ** 2 / (2.0 * sigma ** 2)))

    # Solve Delta(t_half) = Delta0 exactly (Delta is strictly increasing).
    upper = 2.0 * (Delta0 / A + a * sigma) + 10.0 * sigma
    t_half = brentq(lambda t: float(Delta(t)) - Delta0, 1e-12, upper,
                    xtol=1e-14, rtol=1e-15)
    return Schedule("gauss_flat", Delta, dDelta, t_half, sigma,
                    {"nu_star": nu, "Delta0": Delta0, "a": a,
                     "sigma": sigma, "A": A})


def make_gauss_flat_matched_T(T, Delta0, a, sigma_frac):
    """gauss_flat schedule with prescribed total time T and range +-Delta0
    (for the matched-T comparison); returns schedule; nu* = A(1-a)."""
    sigma = sigma_frac * T
    th = T / 2.0
    # A from Delta(T/2) = Delta0:
    denom = th - a * sigma * math.sqrt(math.pi / 2.0) * erf(
        th / (sigma * math.sqrt(2.0)))
    A = Delta0 / denom
    nu = A * (1.0 - a)

    def Delta(t):
        t = np.asarray(t, dtype=float)
        return A * (t - a * sigma * math.sqrt(math.pi / 2.0)
                    * erf(t / (sigma * math.sqrt(2.0))))

    def dDelta(t):
        t = np.asarray(t, dtype=float)
        return A * (1.0 - a * np.exp(-t ** 2 / (2.0 * sigma ** 2)))

    return Schedule("gauss_flat_matchedT", Delta, dDelta, th, sigma,
                    {"nu_star": nu, "Delta0": Delta0, "a": a,
                     "sigma": sigma, "A": A, "T": T})


def build_schedule(kind, nu, Delta0, Omega):
    if kind == "linear":
        return make_linear(nu, Delta0)
    if kind == "tanh":
        return make_tanh(nu, Delta0)
    if kind == "gauss_flat":
        return make_gauss_flat(nu, Delta0, Omega)
    raise ValueError(kind)


# ------------------------------------------------------------- physics

def B_matrix(Delta_val, Omega):
    return 0.5 * (Delta_val * SIGMA_Z + Omega * SIGMA_X)


def run_tdse(sched, Omega, rtol=RTOL, atol=ATOL):
    """Integrate i dpsi/dt = B(t) psi from -t_half to +t_half starting in
    the adiabatic ground state; return leakage = final excited-adiabatic
    population, plus the norm error."""
    th = sched.t_half

    def rhs(t, y):
        d = float(sched.Delta(t))
        return (-0.5j) * np.array([d * y[0] + Omega * y[1],
                                   Omega * y[0] - d * y[1]])

    w0, v0 = np.linalg.eigh(B_matrix(float(sched.Delta(-th)), Omega))
    psi0 = v0[:, 0].astype(complex)          # adiabatic ground state
    sol = solve_ivp(rhs, (-th, th), psi0, method="DOP853",
                    rtol=rtol, atol=atol, dense_output=False)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    psiT = sol.y[:, -1]
    norm_err = abs(np.linalg.norm(psiT) - 1.0)
    wT, vT = np.linalg.eigh(B_matrix(float(sched.Delta(th)), Omega))
    P = float(abs(np.vdot(vT[:, 1], psiT)) ** 2)   # excited-adiabatic pop.
    return P, norm_err, sol.nfev


def lambda_profile(sched, Omega, n=N_LAMBDA_GRID):
    """Lambda(t) on a grid: analytic (Lemma 1) and numeric
    (4th-order central FD of the B matrix -> commutator -> 2-norm)."""
    th = sched.t_half
    t = np.linspace(-th, th, n)
    lam_analytic = 0.5 * Omega * np.abs(sched.dDelta(t))

    h = FD_H_FACTOR * sched.timescale
    # 4th-order central difference of Delta(t) (B is affine in Delta)
    dD_fd = (-sched.Delta(t + 2 * h) + 8.0 * sched.Delta(t + h)
             - 8.0 * sched.Delta(t - h) + sched.Delta(t - 2 * h)) / (12.0 * h)
    # Build full matrices, take the commutator, take the operator norm --
    # tests the whole Lemma-1 chain, not just the scalar derivative.
    D = np.asarray(sched.Delta(t))
    Bt = 0.5 * (D[:, None, None] * SIGMA_Z + Omega * SIGMA_X)
    Bdot = 0.5 * dD_fd[:, None, None] * SIGMA_Z
    C = Bt @ Bdot - Bdot @ Bt          # real 2x2, batched
    lam_num = np.linalg.svd(C, compute_uv=False)[:, 0]
    return t, lam_analytic, lam_num, h


def sc1_check(lam_analytic, lam_num, edge_exclude=SC1_EDGE_EXCLUDE):
    n = lam_analytic.size
    k = int(round(edge_exclude * n))
    sl = slice(k, n - k)
    rel = np.abs(lam_num - lam_analytic) / lam_analytic
    return float(np.max(rel[sl])), float(np.max(rel))


# ----------------------------------------------------------- telemetry

class Telemetry:
    def __init__(self, tel_dir):
        self.paths = {k: os.path.join(tel_dir, f"{k}.jsonl")
                      for k in ("capsules", "operations", "defects")}
        for p in self.paths.values():
            open(p, "w").close()          # fresh sidecars for this run

    def cap(self, rec):
        st_emit.emit(st_emit.capsule(**rec), self.paths["capsules"])
        return rec["capsule_id"]

    def op(self, rec):
        st_emit.emit(st_emit.operation(**rec), self.paths["operations"])
        return rec["operation_id"]

    def defect(self, rec):
        st_emit.emit(rec, self.paths["defects"])
        return rec["defect_id"]


PAULI_FRAME = {"basis_id": "pauli:{I,sigma_x,sigma_y,sigma_z}",
               "coordinate_chart": "diabatic (sigma_z eigenbasis)",
               "version": "1"}


def emit_run_telemetry(tel, sched, Omega, run_tag, P_sim, norm_err):
    """Tier-3 telemetry for one representative run."""
    th = sched.t_half
    ts = np.linspace(-th, th, N_TELEMETRY_SAMPLES)
    h = FD_H_FACTOR * sched.timescale

    # --- hbar-normalization operation (declared unit change J -> rad/s)
    capH = tel.cap(dict(
        capsule_id=f"cap-H-{run_tag}",
        value_ref=f"symbolic:H(t)=(hbar/2)*(Delta(t)*sigma_z+Omega*sigma_x);"
                  f" schedule={sched.name}; Omega={Omega}; params={sched.params}",
        object_kind="operator_family",
        carrier_space="Herm(C^2), t in [-T/2, T/2]",
        units="J",
        frame=PAULI_FRAME,
        law_profile={"hermiticity": {"status": "proven_true",
                                     "evidence": "theorem",
                                     "scope": "class"}},
        provenance={"run_id": RUN_ID, "step_index": 0,
                    "created_at": "2026-07-09"}))
    capBfam = tel.cap(dict(
        capsule_id=f"cap-Bfam-{run_tag}",
        value_ref=f"symbolic:B(t)=H(t)/hbar; schedule={sched.name};"
                  f" Omega={Omega}",
        object_kind="operator_family",
        carrier_space="Herm(C^2), t in [-T/2, T/2]",
        units="rad/s (nondimensionalized, Omega_ref = 1 rad/s)",
        frame=PAULI_FRAME,
        law_profile={"hermiticity": {"status": "proven_true",
                                     "evidence": "theorem",
                                     "scope": "class"}},
        provenance={"run_id": RUN_ID, "step_index": 0}))
    tel.op(dict(
        operation_id=f"op-hbarnorm-{run_tag}",
        operation_name="hbar_normalization: B(t) = H(t)/hbar",
        run_id=RUN_ID, step_index=0,
        input_capsules=[capH], output_capsules=[capBfam],
        laws_used=["hermiticity preserved under scalar division",
                   "spectrum scales by 1/hbar (J -> rad/s)"],
        preconditions=["hbar > 0 (unit conversion, invertible)"],
        side_channels={"hbar": "1 in the nondimensional unit system; "
                               "inverse operation is multiplication by hbar",
                       "units_before": "J", "units_after": "rad/s"},
        severity="observe"))

    # --- per-sample capsules + defects
    defect_ids = []
    prev_ev, prev_cap = None, None
    for k, tk in enumerate(ts):
        d = float(sched.Delta(tk))
        dd_fd = float((-sched.Delta(tk + 2 * h) + 8.0 * sched.Delta(tk + h)
                       - 8.0 * sched.Delta(tk - h)
                       + sched.Delta(tk - 2 * h)) / (12.0 * h))
        B = B_matrix(d, Omega)
        Bdot = 0.5 * dd_fd * SIGMA_Z
        ev = np.linalg.eigvalsh(B)
        gap = float(ev[1] - ev[0])
        capB = tel.cap(dict(
            capsule_id=f"cap-B-{run_tag}-s{k:03d}",
            value_ref=f"inline:B=0.5*(Delta*sigma_z+Omega*sigma_x);"
                      f" t={tk:.9g}; Delta={d:.12g}; Omega={Omega}",
            object_kind="operator",
            carrier_space="Herm(C^2)",
            shape="2x2", units="rad/s",
            frame=PAULI_FRAME,
            spectral_signature={
                "eigenvalues": [float(ev[0]), float(ev[1])],
                "spectral_gaps": [gap],
                "pairing_rule": "sorted_by_real",
                "normality": "normal (Hermitian)",
                "reference_frame": PAULI_FRAME["basis_id"],
                "pairing_conditioning": "gap >= Omega = %g rad/s at the "
                                        "avoided crossing; sorted_by_real "
                                        "pairing well-conditioned" % Omega},
            law_profile={"hermiticity": {"status": "proven_true",
                                         "evidence": "theorem",
                                         "scope": "object"}},
            non_assumptions=["commutativity of B(t) with B(t')"],
            provenance={"run_id": RUN_ID, "step_index": k,
                        "t": float(tk)}))
        capBdot = tel.cap(dict(
            capsule_id=f"cap-Bdot-{run_tag}-s{k:03d}",
            value_ref=f"inline:Bdot_fd=0.5*dDelta_fd*sigma_z;"
                      f" t={tk:.9g}; dDelta_fd={dd_fd:.12g}; fd_order=4;"
                      f" fd_step={h:.9g}",
            object_kind="operator",
            carrier_space="Herm(C^2)",
            shape="2x2", units="rad/s^2",
            frame=PAULI_FRAME,
            provenance={"run_id": RUN_ID, "step_index": k, "t": float(tk),
                        "derived_by": f"op-fd-{run_tag}"}))
        # commutator defect [B, Bdot]: both endomorphisms on Herm(C^2),
        # normalized value dimensionless (units rad/s * rad/s^2 cancel).
        dc = st_emit.commutator_defect(B, Bdot,
                                       object_ids=(capB, capBdot))
        lam_op = float(np.linalg.norm(B @ Bdot - Bdot @ B, 2))
        dc["payload"] = {
            "t": float(tk), "run_id": RUN_ID, "step_index": k,
            "lambda_operator_norm_rad2_per_s3": lam_op,
            "lambda_analytic_lemma1": 0.5 * Omega * abs(
                float(sched.dDelta(tk))),
            "note": "Lambda(t) = ||[B,Bdot]||_op; raw_value is Frobenius "
                    "(= sqrt(2)*operator norm for this commutator); "
                    "normalized_value is the units-invariant "
                    "noncommutativity fraction"}
        defect_ids.append(tel.defect(dc))
        if prev_ev is not None:
            dr = st_emit.eigenvalue_drift_defect(
                ev, prev_ev, object_ids=(capB, prev_cap),
                pairing="sorted_by_real")
            dr["payload"] = {"run_id": RUN_ID, "step_index": k,
                             "gap_lower_bound_rad_per_s": Omega,
                             "note": "avoided crossing: gap = "
                                     "sqrt(Delta^2+Omega^2) >= Omega, "
                                     "pairing well-conditioned"}
            defect_ids.append(tel.defect(dr))
        prev_ev, prev_cap = ev, capB

    # --- finite-difference derivative operation
    tel.op(dict(
        operation_id=f"op-fd-{run_tag}",
        operation_name="finite_difference_derivative (4th-order central)",
        run_id=RUN_ID, step_index=0,
        input_capsules=[capBfam],
        output_capsules=[f"cap-Bdot-{run_tag}-s{k:03d}"
                         for k in range(N_TELEMETRY_SAMPLES)],
        preconditions=[f"Delta(t) smooth; fd step h = {h:.6g} s = 1e-3 * "
                       f"schedule timescale {sched.timescale:.6g} s"],
        laws_used=["4th-order truncation error O(h^4 * d5Delta/dt5)"],
        defects_measured=defect_ids,
        severity="observe"))

    # --- TDSE integration operation + endpoint state capsules
    w0, v0 = np.linalg.eigh(B_matrix(float(sched.Delta(-th)), Omega))
    cap_psi0 = tel.cap(dict(
        capsule_id=f"cap-psi0-{run_tag}",
        value_ref=f"inline:psi0 = adiabatic ground state of B(-T/2); "
                  f"coeffs(diabatic)={[float(x) for x in v0[:, 0]]}",
        object_kind="state_vector",
        carrier_space="C^2 (unit sphere)",
        units="dimensionless amplitude",
        frame=PAULI_FRAME,
        provenance={"run_id": RUN_ID, "step_index": 0}))
    cap_psiT = tel.cap(dict(
        capsule_id=f"cap-psiT-{run_tag}",
        value_ref=f"results/lzsm_sweeps.json#representative_runs/{run_tag}",
        object_kind="state_vector",
        carrier_space="C^2 (unit sphere)",
        units="dimensionless amplitude",
        frame=PAULI_FRAME,
        provenance={"run_id": RUN_ID,
                    "step_index": N_TELEMETRY_SAMPLES - 1,
                    "norm_error": norm_err}))
    tel.op(dict(
        operation_id=f"op-tdse-{run_tag}",
        operation_name="tdse_integration (solve_ivp DOP853)",
        run_id=RUN_ID, step_index=N_TELEMETRY_SAMPLES - 1,
        input_capsules=[capBfam, cap_psi0],
        output_capsules=[cap_psiT],
        preconditions=[f"rtol={RTOL}, atol={ATOL}"],
        laws_used=["unitarity of exp(-i integral B dt) (tested: SC3)"],
        laws_tested=[f"norm preservation: |norm(psiT)-1| = {norm_err:.3e} "
                     f"< {SC3_GATE}"],
        severity="observe"))

    # --- adiabatic projection / measurement quotient
    tel.op(dict(
        operation_id=f"op-adproj-{run_tag}",
        operation_name="adiabatic_projection_measurement: "
                       "P = |<excited_adiabatic(T/2)|psi(T/2)>|^2",
        run_id=RUN_ID, step_index=N_TELEMETRY_SAMPLES - 1,
        input_capsules=[cap_psiT],
        output_capsules=[],
        information_lost={
            "lost_components": "relative phase between the two adiabatic "
                               "amplitudes (Stueckelberg phase)",
            "preserved_components": "excited-adiabatic population P_LZ",
            "loss_fraction": 0.5},
        side_channels={"note": "phase NOT side-channeled -- this is a "
                               "declared destructive quotient; repair "
                               "requires re-running the integration"},
        claim_effects=[{"claim_type": "phase_coherence",
                        "effect": "invalidate",
                        "severity": "observe",
                        "downstream_warning":
                            "P_LZ is a population; interference/"
                            "Stueckelberg phases are not recoverable "
                            "from it"}],
        severity="observe"))


def validate_telemetry(tel):
    import jsonschema
    schema_map = {"capsules": "capsule.schema.json",
                  "operations": "operation_record.schema.json",
                  "defects": "defect_record.schema.json"}
    counts, errors = {}, []
    for kind, schema_file in schema_map.items():
        with open(os.path.join(SCHEMA_DIR, schema_file),
                  encoding="utf-8") as f:
            schema = json.load(f)
        n = 0
        with open(tel.paths[kind], encoding="utf-8") as f:
            for line_no, line in enumerate(f, 1):
                if not line.strip():
                    continue
                n += 1
                try:
                    jsonschema.validate(json.loads(line), schema)
                except jsonschema.ValidationError as e:
                    errors.append(f"{kind}:{line_no}: {e.message}")
        counts[kind] = n
    return counts, errors


# ---------------------------------------------------------------- main

def theory_P(Omega, nu):
    x_std = math.pi * Omega ** 2 / (2.0 * nu)
    return math.exp(-x_std), math.exp(-2.0 * x_std)   # (P_std, P_paper)


def main():
    rng = np.random.default_rng(SEED)  # deterministic pipeline; recorded
    _ = rng.random()                   # touch it so the seed is meaningful
    t_start = time.time()
    tel = Telemetry(TEL_DIR)

    rows = []
    sc1_worst_interior, sc1_worst_all = 0.0, 0.0
    sc3_worst = 0.0
    rep_profiles = {}       # schedule -> profile arrays for the figure
    rep_x_index = 5         # mid ramp rate for representative/telemetry runs
    total_runs = 0

    for kind in ("linear", "tanh", "gauss_flat"):
        for Omega in OMEGAS:
            Delta0 = DELTA0_OVER_OMEGA * Omega
            for ix, x in enumerate(X_GRID):
                nu = math.pi * Omega ** 2 / (2.0 * x)
                sched = build_schedule(kind, nu, Delta0, Omega)
                t_grid, lam_a, lam_n, h_fd = lambda_profile(sched, Omega)
                e_int, e_all = sc1_check(lam_a, lam_n)
                sc1_worst_interior = max(sc1_worst_interior, e_int)
                sc1_worst_all = max(sc1_worst_all, e_all)

                P_sim, norm_err, nfev = run_tdse(sched, Omega)
                sc3_worst = max(sc3_worst, norm_err)
                P_std, P_paper = theory_P(Omega, nu)
                Lambda_star = 0.5 * Omega * nu
                rows.append({
                    "schedule": kind, "Omega": Omega, "Delta0": Delta0,
                    "x_target_std_exponent": float(x), "nu_star": nu,
                    "Lambda_star": Lambda_star,
                    "Lambda_star_over_Omega3": Lambda_star / Omega ** 3,
                    "T_total": 2.0 * sched.t_half,
                    "schedule_params": {k2: float(v2) for k2, v2
                                        in sched.params.items()},
                    "P_sim": P_sim, "P_theory_std": P_std,
                    "P_theory_paper": P_paper,
                    "ln_P_sim": math.log(P_sim),
                    "dlnP_vs_std": math.log(P_sim) - math.log(P_std),
                    "dlnP_vs_paper": math.log(P_sim) - math.log(P_paper),
                    "norm_error": norm_err, "solver_nfev": int(nfev),
                    "sc1_max_rel_err_interior": e_int,
                    "sc1_max_rel_err_all": e_all,
                    "fd_step": h_fd})
                total_runs += 1
                if Omega == 1.0 and ix == rep_x_index:
                    ss = max(1, N_LAMBDA_GRID // 20)
                    rep_profiles[kind] = {
                        "t_over_T": (t_grid / (2.0 * sched.t_half)).tolist(),
                        "Delta_over_Delta0":
                            (np.asarray(sched.Delta(t_grid))
                             / Delta0).tolist(),
                        "lam_analytic_over_peak":
                            (lam_a / (0.5 * Omega * nu)).tolist(),
                        "lam_fd_over_peak_subsampled":
                            (lam_n[::ss] / (0.5 * Omega * nu)).tolist(),
                        "t_over_T_subsampled":
                            (t_grid[::ss] / (2.0 * sched.t_half)).tolist(),
                        "nu_star": nu, "P_sim": P_sim}
                    emit_run_telemetry(tel, sched, Omega,
                                       f"{kind}-rep", P_sim, norm_err)
                print(f"  {kind:11s} Omega={Omega:3.1f} x={x:6.3f} "
                      f"nu*={nu:8.4f} P_sim={P_sim:.6e} "
                      f"P_std={P_std:.6e} dlnP_std="
                      f"{math.log(P_sim) - math.log(P_std):+.4f}",
                      flush=True)

    # ---------------- SC2 summary (reported, not gated)
    def sc2_block(rs):
        sel = [r for r in rs if r["P_sim"] > 1e-6]
        return {
            "n_points": len(sel),
            "max_abs_dlnP_vs_std": max(abs(r["dlnP_vs_std"]) for r in sel),
            "max_abs_dlnP_vs_paper": max(abs(r["dlnP_vs_paper"])
                                         for r in sel),
            "max_abs_dlog10P_vs_std": max(abs(r["dlnP_vs_std"])
                                          for r in sel) / math.log(10),
            "max_abs_dlog10P_vs_paper": max(abs(r["dlnP_vs_paper"])
                                            for r in sel) / math.log(10)}

    def fitted_slope(rs):
        xs = np.array([r["x_target_std_exponent"] for r in rs])
        ys = np.array([-r["ln_P_sim"] for r in rs])
        return float(np.sum(xs * ys) / np.sum(xs * xs))

    sc2 = {"pooled": sc2_block(rows),
           "per_schedule": {k: sc2_block([r for r in rows
                                          if r["schedule"] == k])
                            for k in ("linear", "tanh", "gauss_flat")},
           "fitted_exponent_slope_b": {
               "pooled": fitted_slope(rows),
               **{k: fitted_slope([r for r in rows if r["schedule"] == k])
                  for k in ("linear", "tanh", "gauss_flat")}},
           "interpretation": "b ~ 1 => standard form "
                             "exp(-pi Omega^3/(4 Lambda*)); b ~ 2 => "
                             "paper form exp(-pi Omega^3/(2 Lambda*))"}

    # ---------------- matched-T comparison (abstract's claim)
    matched = []
    Omega = 1.0
    Delta0 = DELTA0_OVER_OMEGA * Omega
    for x_lin in (0.5, 1.0, 2.0):
        nu_lin = math.pi * Omega ** 2 / (2.0 * x_lin)
        T = 2.0 * Delta0 / nu_lin
        lin = make_linear(nu_lin, Delta0)
        flat = make_gauss_flat_matched_T(T, Delta0, a=0.85, sigma_frac=0.10)
        P_lin, ne1, _ = run_tdse(lin, Omega)
        P_flat, ne2, _ = run_tdse(flat, Omega)
        sc3_worst = max(sc3_worst, ne1, ne2)
        nu_flat = flat.params["nu_star"]
        Pl_std, _ = theory_P(Omega, nu_lin)
        Pf_std, _ = theory_P(Omega, nu_flat)
        tau_lz = max(1.0 / math.sqrt(nu_flat), Omega / nu_flat)
        matched.append({
            "Omega": Omega, "Delta0": Delta0, "T_total": T,
            "x_linear": x_lin,
            "linear": {"nu_star": nu_lin,
                       "Lambda_star": 0.5 * Omega * nu_lin,
                       "P_sim": P_lin, "P_theory_std": Pl_std},
            "gauss_flat": {"nu_star": nu_flat,
                           "Lambda_star": 0.5 * Omega * nu_flat,
                           "a": 0.85, "sigma": flat.params["sigma"],
                           "sigma_over_tau_LZ": flat.params["sigma"] / tau_lz,
                           "P_sim": P_flat, "P_theory_std": Pf_std},
            "leakage_reduction_measured": P_lin / P_flat,
            "leakage_reduction_predicted_std": Pl_std / Pf_std,
            "locality_caveat": "LZ point prediction degrades when "
                               "sigma ~ tau_LZ; measured P is the "
                               "authoritative number"})
        print(f"  matched-T x_lin={x_lin}: P_lin={P_lin:.4e} "
              f"P_flat={P_flat:.4e} reduction={P_lin / P_flat:.1f}x",
              flush=True)

    # ---------------- SC5 solver convergence (reported)
    sc5_rows = []
    for kind in ("linear", "tanh", "gauss_flat"):
        for Omega in (0.5, 2.0):
            Delta0 = DELTA0_OVER_OMEGA * Omega
            x = float(X_GRID[-1])
            nu = math.pi * Omega ** 2 / (2.0 * x)
            sched = build_schedule(kind, nu, Delta0, Omega)
            P_tight, _, _ = run_tdse(sched, Omega)
            P_loose, _, _ = run_tdse(sched, Omega,
                                     rtol=RTOL_LOOSE, atol=ATOL_LOOSE)
            sc5_rows.append({"schedule": kind, "Omega": Omega,
                             "x": x, "P_tight": P_tight,
                             "P_loose": P_loose,
                             "abs_dlnP": abs(math.log(P_tight)
                                             - math.log(P_loose))})
    sc5_max = max(r["abs_dlnP"] for r in sc5_rows)

    # ---------------- SC4 telemetry validation (gated)
    tel_counts, tel_errors = validate_telemetry(tel)

    # ---------------- results JSON capsule (Tier 0/1)
    tel.cap(dict(
        capsule_id=f"cap-results-{RUN_ID}",
        value_ref="results/lzsm_sweeps.json",
        object_kind="results_table",
        carrier_space="JSON document",
        units="probabilities dimensionless; rates rad/s^2; "
              "curvature rad^2/s^3",
        provenance={"run_id": RUN_ID, "step_index": total_runs,
                    "created_at": "2026-07-09"}))
    # re-validate after the late capsule
    tel_counts, tel_errors = validate_telemetry(tel)

    sc1_pass = sc1_worst_interior < SC1_GATE
    sc3_pass = sc3_worst < SC3_GATE
    sc4_pass = not tel_errors

    results = {
        "meta": {
            "run_id": RUN_ID, "seed": SEED, "date": "2026-07-09",
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "platform": platform.platform(),
            "solver": {"method": "DOP853", "rtol": RTOL, "atol": ATOL},
            "units": "nondimensional, Omega_ref = 1 rad/s; hbar = 1 "
                     "after declared hbar-normalization B = H/hbar",
            "conventions": {
                "B": "B(t) = H(t)/hbar = 0.5*(Delta(t) sigma_z "
                     "+ Omega sigma_x)",
                "Lambda": "Lambda(t) = ||[B, Bdot]||_op "
                          "= 0.5 * Omega * |Delta_dot| (Lemma 1)",
                "P_theory_std": "exp(-pi Omega^2/(2 nu*)) "
                                "= exp(-pi Omega^3/(4 Lambda*)) "
                                "[standard LZ]",
                "P_theory_paper": "exp(-pi Omega^2/nu*) "
                                  "= exp(-pi Omega^3/(2 Lambda*)) "
                                  "[paper Eqs. 1/5 as printed]"},
            "runtime_seconds": None},
        "preregistration": {
            "SC1": f"GATED: Lemma-1 FD-vs-analytic max rel err < {SC1_GATE}"
                   f" excluding outermost {SC1_EDGE_EXCLUDE:.0%} of grid "
                   f"points per end",
            "SC2": "REPORTED (not gated): max |ln P_sim - ln P_theory| "
                   "over points with P_sim > 1e-6, both theory forms, "
                   "plus fitted exponent slope b",
            "SC3": f"GATED: | ||psi(T)|| - 1 | < {SC3_GATE} every run",
            "SC4": "GATED: all telemetry records validate against the "
                   "three schemas",
            "SC5": "REPORTED: max |dlnP| tight vs loose solver tolerance "
                   "on 6-run subset"},
        "parameters": {
            "Delta0_over_Omega": DELTA0_OVER_OMEGA,
            "Omegas": list(OMEGAS),
            "x_grid_std_exponents": X_GRID.tolist(),
            "gauss_flat": {"a": GAUSS_A,
                           "sigma_rule": f"{GAUSS_SIGMA_FACTOR} * "
                                         f"max(1/sqrt(nu*), Omega/nu*)"},
            "tanh_domain_half_widths": TANH_HALF_WIDTHS,
            "lambda_grid_points": N_LAMBDA_GRID,
            "fd": {"order": 4, "step_rule": f"{FD_H_FACTOR} * schedule "
                                            f"timescale"}},
        "sweep_rows": rows,
        "matched_T_comparison": matched,
        "self_checks": {
            "SC1": {"max_rel_err_interior": sc1_worst_interior,
                    "max_rel_err_including_endpoints": sc1_worst_all,
                    "gate": SC1_GATE, "pass": bool(sc1_pass)},
            "SC2": sc2,
            "SC3": {"max_norm_error": sc3_worst, "gate": SC3_GATE,
                    "pass": bool(sc3_pass)},
            "SC4": {"record_counts": tel_counts,
                    "validation_errors": tel_errors,
                    "pass": bool(sc4_pass)},
            "SC5": {"rows": sc5_rows, "max_abs_dlnP": sc5_max}},
        "representative_profiles_note":
            "Lambda(t)/Delta(t) profiles for the figure are the "
            "Omega = 1, mid-rate runs (x = %.4f)" % float(X_GRID[rep_x_index]),
    }
    results["meta"]["runtime_seconds"] = round(time.time() - t_start, 1)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "lzsm_sweeps.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    with open(os.path.join(RESULTS_DIR, "rep_profiles.json"), "w",
              encoding="utf-8") as f:
        json.dump(rep_profiles, f)

    print("\n================ SELF-CHECK SUMMARY ================")
    print(f"SC1 (Lemma 1, GATED)  max rel err interior = "
          f"{sc1_worst_interior:.3e}  (gate {SC1_GATE})  "
          f"{'PASS' if sc1_pass else 'FAIL'}")
    print(f"     incl. endpoints  = {sc1_worst_all:.3e}")
    print(f"SC2 (reported)  pooled max|dlnP| vs STANDARD form = "
          f"{sc2['pooled']['max_abs_dlnP_vs_std']:.4f}")
    print(f"                pooled max|dlnP| vs PAPER form    = "
          f"{sc2['pooled']['max_abs_dlnP_vs_paper']:.4f}")
    print(f"                fitted exponent slope b (pooled)  = "
          f"{sc2['fitted_exponent_slope_b']['pooled']:.4f}  "
          f"(1 = standard, 2 = paper)")
    print(f"SC3 (unitarity, GATED)  max norm err = {sc3_worst:.3e}  "
          f"{'PASS' if sc3_pass else 'FAIL'}")
    print(f"SC4 (telemetry, GATED)  counts = {tel_counts}  "
          f"errors = {len(tel_errors)}  "
          f"{'PASS' if sc4_pass else 'FAIL'}")
    print(f"SC5 (reported)  max |dlnP| tight-vs-loose = {sc5_max:.3e}")
    print(f"total runs: {total_runs} + {2 * len(matched)} matched-T "
          f"+ {2 * len(sc5_rows)} convergence")
    print(f"runtime: {results['meta']['runtime_seconds']} s")
    print(f"results -> {out_path}")
    return 0 if (sc1_pass and sc3_pass and sc4_pass) else 1


if __name__ == "__main__":
    sys.exit(main())
