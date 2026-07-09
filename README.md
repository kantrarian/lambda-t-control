# lambda-t-control — LZSM simulation for the Spectral Curvature paper

Numerical companion to **"Predictive Geometric Control of Quantum Dynamics
via Spectral Curvature Functional Λ(t)"** (R. J. Mathews). The published
manuscript's Section 5 ("Simulation") was a placeholder; this repo contains
the missing simulation, the figure the placeholder promised, structural-
telemetry sidecars for the numerical protocol, and the recompiled PDF with
Section 5 filled in.

## What this is

A two-level time-dependent Schrödinger (TDSE) study of
Landau–Zener–Stückelberg–Majorana sweeps for

    B(t) = H(t)/ħ = ½ (Δ(t) σ_z + Ω σ_x)        [rad/s]

testing the paper's two claims:

1. **Lemma 1 (exact identity):** Λ(t) = ‖[B, Ḃ]‖ = ½ Ω |Δ̇(t)| — verified
   numerically end-to-end (finite-difference Ḃ → matrix commutator →
   operator norm vs. the closed form) to relative error < 4×10⁻⁹ over all
   grid points of all 90 runs (pre-registered gate: 1×10⁻⁶). **PASS**
2. **Lemma 2 (leakage law):** P_LZ as a function of peak curvature Λ(t*).
   Both candidate exponent constants were simulated (see "Honest results").

Three detuning schedules (linear, tanh, Gaussian-flattened), Ω ∈
{0.5, 1, 2}, ten crossing slopes each: 90 TDSE runs (DOP853,
rtol 1e-11 / atol 1e-13), plus 6 matched-total-time comparison runs and
12 solver-convergence re-runs. Everything is deterministic; seed = 0 is
recorded in the results JSON for any future stochastic extension.

## Honest results summary

| Check | Pre-registered? | Result |
|---|---|---|
| SC1 Lemma-1 identity, FD vs analytic | gated, < 1e-6 | max rel err 2.6e-9 interior, 3.5e-9 incl. endpoints — **PASS** |
| SC2 closed-form agreement | reported, not gated | max \|ln P_sim − ln P_theory\| vs **standard** form: 6e-4 (linear), 2.2e-3 (tanh), 0.084 (gauss-flat); vs **paper** form: 8.0 |
| SC3 unitarity | gated, < 1e-8 | max norm error 1.6e-9 — **PASS** |
| SC4 telemetry schema validation | gated | 211 capsules + 12 operations + 195 defects, 0 errors — **PASS** |
| SC5 solver convergence | reported | max \|Δ ln P\| tight-vs-loose tolerance = 5.8e-7 |

**The factor-of-two finding — CORRECTED in v2 (2026-07-09).** The paper's
original Eq. (1) wrote the LZ probability as exp(−πΩ²/ν); the standard
Landau–Zener result for this Hamiltonian is exp(−πΩ²/(2ν)) (Zener 1932).
The original Eq. (5), P ≈ exp(−πΩ³/2Λ(t*)), **was** algebraically
consistent with its own Eq. (1) via Lemma 1 (ν = 2Λ/Ω), but the standard
result maps to **exp(−πΩ³/4Λ(t*))**. The TDSE data decided it: regressing
−ln P_sim on πΩ²/(2ν*) gives slope b = 1.003 (pooled; 1.000 / 1.000 /
1.010 per schedule) — the data follow the standard 4Λ constant; the
originally printed 2Λ form overstated suppression by a factor 2 in the
exponent (deviation reaches 8 in ln P at the slowest sweeps). Both forms
are overlaid in Fig. 1(c). **With the inventor's approval (2026-07-09),
the correction is propagated through the theory sections in this repo's
tex and PDF** (Eqs. 1, 5, 9–11 now carry the standard constant; the two
spectral lower bounds tighten from 2π²Ω² to π²Ω² in the numerator
accordingly), with a dated revision footnote and Section 5 documenting
the provenance. Lemma 1 and the pointwise control law were never
affected.

**Pointwise-law confirmation.** The Gaussian-flattened schedule is 4×
steeper than its crossing slope everywhere except near t*, yet its leakage
is set by Λ(t*) alone. The residual 0.084 in ln P for that schedule is the
locality limit of the pointwise formula (its flattening window is only
≈ 6 LZ jump times wide by design); linear and tanh agree to < 2.2e-3.

**Matched-time comparison (the abstract's claim).** At equal sweep range
(±30 Ω) and equal total time, replacing a linear ramp with a
curvature-flattened waveform (a = 0.85, σ = 0.1 T) reduced measured
leakage by 8.0× (ΩT = 19.1), 5.7×10³× (ΩT = 38.2), 3.5×10²× (ΩT = 76.4).
Caveat: for these strongly flattened pulses σ ≈ 1.1 τ_LZ, so the pointwise
formula is not quantitative there (predicted reductions 8.4× / 70× /
4.9×10³×) — the measured TDSE numbers are authoritative; the qualitative
design rule (flatten Δ̇ at the crossing) is what survives.

**Ω-collapse note.** The three Ω branches produce identical P at fixed
target exponent to < 6×10⁻¹³ relative. That is exact scale invariance of
the nondimensionalized problem (t → Ωt maps branches onto each other), so
the Ω sweep is a consistency check of the Λ(t*)/Ω³ scaling, not three
independent datasets.

## Scope boundary

- Two-level LZSM only. No multi-qubit blocks, no time-varying Ω(t), no
  counter-diabatic benchmark, no hardware/latency modeling — those are
  outside this simulation and remain future work in the paper's terms.
- The structural-telemetry sidecars are the *working-level* discipline of
  the Structural Telemetry concept paper (capsules / operation records /
  defect records as JSONL, schema-validated); they are not the Phase-1 SDK.
- Nothing here touches the related patent filings; this repo is the public
  manuscript plus its public numerical companion.

## Layout

```
Quantum_Dynamics_Spectral_Curvature.tex   manuscript (Section 5 filled, abstract sim sentence updated)
Quantum_Dynamics_Spectral_Curvature.pdf   recompiled PDF (pdflatex x2)
sim/lzsm_sweeps.py    TDSE harness; pre-registered acceptance criteria in the module docstring
sim/plot_fig.py       figure generator (validated palette, colorblind-safe)
sim/st_emit.py        structural-telemetry helper (pinned copy; `--selftest` passes)
schemas/              ST record schemas (capsule / operation / defect)
results/lzsm_sweeps.json   all sweep numbers, tolerances, seeds, self-check outcomes
results/rep_profiles.json  Λ(t)/Δ(t) profiles behind Fig. 1(a,b)
figures/fig_curvature_sweeps.{png,pdf}
telemetry/{capsules,operations,defects}.jsonl   Tier-3 sidecars for the representative runs
```

## How to run

Requires Python 3.12 with numpy, scipy, matplotlib, jsonschema.

```
cd sim
python st_emit.py --selftest      # confirm the telemetry helper copy
python lzsm_sweeps.py             # ~2 min; writes results/ + telemetry/, prints self-check summary, exit 0 iff gates pass
python plot_fig.py                # writes figures/fig_curvature_sweeps.{png,pdf}
cd ..
pdflatex Quantum_Dynamics_Spectral_Curvature.tex   # run twice for references
```

The pre-registration discipline: every gated threshold above was written
into `sim/lzsm_sweeps.py`'s docstring before the first full run; SC2 was
pre-declared as report-only. Environment (numpy/scipy versions, platform,
runtime) is recorded in `results/lzsm_sweeps.json` under `meta`.
