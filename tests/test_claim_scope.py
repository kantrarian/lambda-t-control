import json
from pathlib import Path
import unittest

import numpy as np

from sim.claim_scope import cubic_crossing, slope_surrogate


class ClaimScopeTests(unittest.TestCase):
    def test_identity_is_operator_norm_not_frobenius(self):
        sx = np.array([[0., 1.], [1., 0.]])
        sz = np.diag([1., -1.])
        for omega, delta, slope in ((1., .2, -3.), (2., 10., .1), (.5, 0., 0.)):
            b, db = .5*(delta*sz + omega*sx), .5*slope*sz
            comm = b @ db - db @ b
            self.assertAlmostEqual(np.linalg.norm(comm, 2), .5*omega*abs(slope))
            self.assertAlmostEqual(np.linalg.norm(comm, "fro"), np.sqrt(2)*.5*omega*abs(slope))

    def test_fourier_bound_implies_surrogate_upper_not_lower(self):
        omega, moment = 1., 10.
        upper = omega * moment / (4*np.pi)
        for fraction in (0., .1, .5, 1.):
            self.assertLessEqual(slope_surrogate(omega, fraction*upper),
                                 np.exp(-np.pi**2*omega**2/moment) + 1e-15)

    def test_cubic_crossing_refutes_universal_slope_law(self):
        loose, tight = cubic_crossing(1e-10, 1e-12), cubic_crossing()
        self.assertEqual(tight["lambda_at_crossing"], 0.)
        self.assertGreater(tight["P_excited"], .10)
        self.assertLess(tight["P_excited"], .12)
        self.assertLess(abs(loose["P_excited"] - tight["P_excited"]), 1e-9)
        self.assertLess(tight["norm_error"], 1e-9)

    def test_gaussian_crossing_is_slope_minimum(self):
        slopes = 1 - .75*np.exp(-np.array([-1., 0., 1.])**2/2)
        self.assertLess(slopes[1], slopes[0])
        self.assertLess(slopes[1], slopes[2])

    def test_archived_counts_and_reported_residuals(self):
        root = Path(__file__).resolve().parents[1]
        data = json.loads((root / "results/lzsm_sweeps.json").read_text())
        rows = data["sweep_rows"]
        self.assertEqual(len(rows), 90)
        for schedule in ("linear", "tanh", "gauss_flat"):
            selected = [r for r in rows if r["schedule"] == schedule]
            self.assertEqual(len(selected), 30)
            measured = max(abs(np.log(r["P_sim"]) - np.log(r["P_theory_std"])) for r in selected)
            stored = data["self_checks"]["SC2"]["per_schedule"][schedule]["max_abs_dlnP_vs_std"]
            self.assertAlmostEqual(measured, stored, places=12)
        self.assertGreater(data["self_checks"]["SC2"]["per_schedule"]["gauss_flat"]["max_abs_dlnP_vs_std"], .08)


if __name__ == "__main__":
    unittest.main()
