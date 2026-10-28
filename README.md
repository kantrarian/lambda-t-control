# lambda-t-control: two-level commutator telemetry

Numerical companion to R. J. Mathews' spectral-curvature research.

**Current technical manuscript:** [v3 corrected paper](Quantum_Dynamics_Spectral_Curvature_v3.pdf)
and [LaTeX source](Quantum_Dynamics_Spectral_Curvature_v3.tex), October 27, 2026.
The unversioned PDF/TeX retain the historical v2 and are superseded for claim
scope by v3. The original results and telemetry remain unchanged.

## What is established

For a fixed Pauli basis and constant positive coupling,
`B=H/hbar=(Delta sigma_z+Omega sigma_x)/2` gives the exact operator-norm identity
`Lambda=||[B,Bdot]||=Omega |Delta_dot|/2`. In this model Lambda is a rescaled
slope, not an independent predictor of transition probability.

For an asymptotic linear crossing, the standard Landau-Zener transition
probability is `exp(-pi Omega^3/(4 Lambda))`. Substituting the crossing slope
of a finite nonlinear pulse gives a **surrogate**, not a universal law.
The Fourier moment bounds Lambda from above, so it bounds that surrogate
from **above**, not below. Without an approximation-error bound this is not a
rigorous bound on actual leakage. A cubic crossing has zero Lambda at the
crossing but about 0.10923 final excited-state probability in the specified
finite-interval test. V3 states the exact assumptions and derivations.

## Archived evidence

The July 2026 experiment includes 90 TDSE sweep runs, six matched-duration
runs and twelve tolerance-check integrations. The 90-run cohort consists of
three schedules, three coupling scales and ten crossing slopes. Coupling
scale collapse is a consistency check, not three independent datasets.

| Check | Archived result |
| --- | --- |
| Exact identity, finite difference versus analytic | maximum relative error 3.5e-9 including endpoints |
| Maximum absolute log-probability residual against slope surrogate | 0.000574 linear; 0.002221 tanh; 0.083920 Gaussian-flattened |
| Norm preservation | maximum error about 1.6e-9 |
| Tolerance sensitivity | maximum absolute change in log probability about 5.8e-7 |
| Telemetry schema checks | 211 capsules, 12 operations, 195 defects; no reported schema errors |
| Three matched-duration comparisons | measured reductions about 8.0, 5713.3, 349.2 |

These observations do not demonstrate a universal pointwise control law,
multilevel leakage suppression, closed-loop hardware operation or quantum
advantage. Gaussian flattening produces a crossing slope minimum, not a peak.
The archived figure's crossing normalization is clarified in the v3 caption.

The original factor-of-two exponent error was corrected in v2. V3 corrects
additional inequality and interpretation errors without changing the old
data or pretending the new scope tests were preregistered.

## Reproduce the checks

Requires numpy and scipy; the archived full sweep also needs jsonschema and
the figure generator needs matplotlib.

```powershell
python -m unittest discover -s tests -v
python sim/claim_scope.py
pdflatex -interaction=nonstopmode -halt-on-error Quantum_Dynamics_Spectral_Curvature_v3.tex
pdflatex -interaction=nonstopmode -halt-on-error Quantum_Dynamics_Spectral_Curvature_v3.tex
```

The five new tests check the operator norm, inequality direction, cubic
counterexample with solver convergence, Gaussian crossing minimum and
consistency of the archived residuals. They do not rerun the full archive.
To rerun the original protocol, use `python sim/lzsm_sweeps.py`; it overwrites
the original generated results and telemetry, so use an isolated checkout.
Its historical names and comments are preserved for provenance, not endorsed
as the current interpretation.

No patent filing is changed by these manuscript corrections.
