"""Post-review scope diagnostic, not part of the July preregistered cohort."""

import json
import math
import numpy as np
from scipy.integrate import solve_ivp


def slope_surrogate(omega, crossing_lambda):
    if not (math.isfinite(omega) and omega > 0
            and math.isfinite(crossing_lambda) and crossing_lambda >= 0):
        raise ValueError("finite positive coupling and nonnegative Lambda required")
    return (math.exp(-math.pi * omega**3 / (4 * crossing_lambda))
            if crossing_lambda > 0 else 0.0)


def cubic_crossing(rtol=1e-12, atol=1e-14):
    sx = np.array([[0., 1.], [1., 0.]])
    sz = np.diag([1., -1.])
    end = 30.**(1/3)

    def b(t):
        return .5 * (t**3 * sz + sx)

    _, initial = np.linalg.eigh(b(-end))
    _, final = np.linalg.eigh(b(end))
    solved = solve_ivp(lambda t, psi: -1j * b(t) @ psi, (-end, end),
                       initial[:, 0].astype(complex), method="DOP853",
                       rtol=rtol, atol=atol)
    if not solved.success:
        raise RuntimeError(solved.message)
    psi = solved.y[:, -1]
    return {"schedule": "Delta=t^3", "Omega": 1., "interval": [-end, end],
            "initial_state": "instantaneous ground", "measurement": "final instantaneous excited",
            "lambda_at_crossing": 0., "P_excited": float(abs(np.vdot(final[:, 1], psi))**2),
            "norm_error": float(abs(np.linalg.norm(psi) - 1)), "rtol": rtol, "atol": atol}


if __name__ == "__main__":
    print(json.dumps({"status": "post_review_scope_check_not_preregistered",
                      "cases": [cubic_crossing(1e-10, 1e-12), cubic_crossing()]}, indent=2))
