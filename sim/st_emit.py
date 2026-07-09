# -*- coding: utf-8 -*-
"""st_emit.py -- Structural Telemetry emission helper (v0.2).

Working-level implementation of the Structural Telemetry concept paper
v1.1 (kantrarian/Structural_Telemetry): capsule / operation-record /
defect-record builders, JSONL emission, and defect calculators enforcing
the paper's well-posedness rules (norm_type, undefined_by_type, eigenvalue
pairing, principal-angle aggregation -- see SKILL.md).

Two kinds of function:
  RECORD PRODUCERS return complete DefectRecords that validate against
  ../schemas/defect_record.schema.json: commutator_defect, grade_leak,
  projection_loss, eigenvalue_drift_defect, frame_drift_defect.
  MEASUREMENT HELPERS return plain values/dicts for embedding in
  signatures or payloads: grade_norms_from_coeffs, grade_energy_fractions,
  grade_entropy, principal_angles.

Usage pattern: COPY this file into the experiment directory (pins the
version that ran), then `python st_emit.py --selftest`. The self-test
reproduces the paper's Appendix B/C worked numbers exactly and, when
`jsonschema` is installed, validates every record producer's output
(including undefined_by_type branches) against the schemas. If it fails
after an edit, the helper has diverged from the paper -- fix whichever is
wrong, never both silently.

Requires numpy. Ad-hoc extra fields belong in payload (defects),
provenance (capsules), or side_channels/frame_change/... (operations);
the top-level schemas are closed.
"""
import json
import math
import sys
import uuid

import numpy as np

SEVERITIES = ("observe", "warn", "invalidate", "reject", "quarantine")
DEFECT_STATUSES = ("measured_zero", "measured_nonzero", "undefined_by_type",
                   "skipped", "ill_conditioned")
EIGENVALUE_PAIRINGS = ("sorted_by_real", "sorted_by_abs")


# ---------------------------------------------------------------- algebra

def grade_dims(n):
    """Grade dimensions of an n-generator Clifford/exterior algebra:
    binom(n, k) for k = 0..n. Cl(3,1) -> [1, 4, 6, 4, 1]."""
    return [math.comb(n, k) for k in range(n + 1)]


def grade_slices(n):
    """(start, end) coefficient slice per grade, assuming coefficients are
    stored grade-contiguously in the declared basis-blade order. If your
    layout differs, pass explicit per-grade norms instead -- the layout is
    frame metadata and must match the capsule's declared frame."""
    dims, out, pos = grade_dims(n), [], 0
    for d in dims:
        out.append((pos, pos + d))
        pos += d
    return out


def grade_norms_from_coeffs(coeffs, n=4, norm_type="euclidean_coefficient"):
    """Per-grade norms of a flat coefficient vector (measurement helper).

    Only euclidean_coefficient is implemented: the invariant form is
    indefinite in Cl(3,1) (null blades norm to zero -- false negatives for
    leakage) and needs the algebra's metric, which this helper does not
    carry. Declare norm_type in the emitted GradeSignature either way.
    """
    if norm_type != "euclidean_coefficient":
        raise NotImplementedError(
            "only euclidean_coefficient is implemented; see SKILL.md "
            "'Grade norms need a declared norm_type'")
    coeffs = np.asarray(coeffs, dtype=float).ravel()
    expected = sum(grade_dims(n))
    if coeffs.size != expected:
        raise ValueError(f"expected {expected} coefficients for n={n}, "
                         f"got {coeffs.size}")
    return {k: float(np.linalg.norm(coeffs[a:b]))
            for k, (a, b) in enumerate(grade_slices(n))}


def grade_energy_fractions(grade_norms):
    """Squared-norm energy shares -- the paper Appendix B
    'normalized_grade_energy' estimator (measurement helper)."""
    energies = {int(k): float(v) ** 2 for k, v in grade_norms.items()}
    total = sum(energies.values())
    if total <= 0.0:
        raise ValueError("zero total energy -- entropy undefined_by_type")
    return {k: e / total for k, e in energies.items()}


def grade_entropy(grade_norms, log_base=2):
    """Shannon entropy of the grade-energy distribution (measurement
    helper). RV = grade_energy_distribution, estimator =
    normalized_grade_energy."""
    fracs = grade_energy_fractions(grade_norms)
    return -sum(p * math.log(p, log_base) for p in fracs.values() if p > 0.0)


# ------------------------------------------------------- defect calculators

def _defect_base(defect_type, object_ids, **extra):
    """Assemble a DefectRecord, dropping None-valued keys so that
    undefined/skipped branches never trip the schema's presence-based
    dependentRequired constraints."""
    rec = {"defect_id": f"def-{uuid.uuid4().hex}",
           "defect_type": defect_type,
           "object_ids": [str(o) for o in object_ids]}
    rec.update({k: v for k, v in extra.items() if v is not None})
    return rec


def commutator_defect(A, B, object_ids=("A", "B"), eps=0.0, tolerance=1e-12):
    """[A, B] defect (record producer) with the paper's typing rule: both
    operands must be endomorphisms on a common carrier (square, same shape)
    or the defect is undefined_by_type -- never reported as zero. NOTE
    (SKILL.md): the value is not invariant under per-component unit
    rescalings; only compute this after the declared normalization that
    makes A, B dimensionless."""
    A, B = np.asarray(A, dtype=float), np.asarray(B, dtype=float)
    if A.ndim != 2 or B.ndim != 2 or A.shape[0] != A.shape[1] \
            or B.shape[0] != B.shape[1] or A.shape != B.shape:
        return _defect_base("noncommutativity", object_ids,
                            defect_status="undefined_by_type",
                            payload={"reason": "operands are not "
                                     "endomorphisms on a common carrier",
                                     "shapes": [list(A.shape), list(B.shape)]})
    prod_norms = float(np.linalg.norm(A, "fro") * np.linalg.norm(B, "fro"))
    if prod_norms <= 0.0:
        # test BEFORE adding eps: eps must never turn "untestable" into
        # "law holds" (SKILL.md status-vocabulary gotcha)
        return _defect_base("noncommutativity", object_ids,
                            defect_status="undefined_by_type",
                            payload={"reason": "zero-norm operand"})
    raw = float(np.linalg.norm(A @ B - B @ A, "fro"))
    normalized = raw / (prod_norms + eps)
    status = "measured_zero" if normalized <= tolerance else "measured_nonzero"
    return _defect_base("noncommutativity", object_ids,
                        defect_status=status, raw_value=raw,
                        norm_type="frobenius",
                        normalizer="norm(A)*norm(B)+eps", eps=eps,
                        normalized_value=normalized, tolerance=tolerance)


def grade_leak(coeffs, allowed_grades, n=4, object_ids=("M",),
               tolerance=1e-12):
    """norm(proj_forbidden(M)) / norm(M) (record producer),
    euclidean_coefficient norm (declared), undefined_by_type on a zero-norm
    operand (paper 6.4)."""
    norms = grade_norms_from_coeffs(coeffs, n=n)
    allowed = sorted(int(g) for g in allowed_grades)
    total = math.sqrt(sum(v ** 2 for v in norms.values()))
    forbidden = math.sqrt(sum(v ** 2 for k, v in norms.items()
                              if k not in set(allowed)))
    if total <= 0.0:
        return _defect_base("grade_leakage", object_ids,
                            defect_status="undefined_by_type",
                            payload={"reason": "zero-norm multivector",
                                     "allowed_grades": allowed})
    normalized = forbidden / total
    status = "measured_zero" if normalized <= tolerance else "measured_nonzero"
    return _defect_base("grade_leakage", object_ids, defect_status=status,
                        raw_value=forbidden, norm_type="euclidean_coefficient",
                        normalizer="norm(M)", normalized_value=normalized,
                        tolerance=tolerance,
                        payload={"allowed_grades": allowed})


def projection_loss(x, x_recon, object_ids=("X",), tolerance=1e-12):
    """norm(x - lift(P x)) / norm(x) against a provided reconstruction
    (record producer)."""
    x = np.asarray(x, dtype=float).ravel()
    r = np.asarray(x_recon, dtype=float).ravel()
    total = float(np.linalg.norm(x))
    if total <= 0.0:
        return _defect_base("projection_loss", object_ids,
                            defect_status="undefined_by_type",
                            payload={"reason": "zero-norm input"})
    raw = float(np.linalg.norm(x - r))
    normalized = raw / total
    status = "measured_zero" if normalized <= tolerance else "measured_nonzero"
    return _defect_base("projection_loss", object_ids, defect_status=status,
                        raw_value=raw, norm_type="euclidean_coefficient",
                        normalizer="norm(X)", normalized_value=normalized,
                        tolerance=tolerance)


def eigenvalue_drift_defect(ev_t, ev_ref, object_ids=("A_t", "A_ref"),
                            pairing="sorted_by_real", tolerance=1e-12):
    """Paired eigenvalue drift under a DECLARED pairing rule (paper 6.6;
    record producer). Near crossings/degeneracy prefer subspace summaries
    (frame_drift_defect). A count mismatch is a DimensionEvent, reported
    here as undefined_by_type with the reason in payload."""
    if pairing not in EIGENVALUE_PAIRINGS:
        raise ValueError(f"pairing must be one of {EIGENVALUE_PAIRINGS}")
    key = {"sorted_by_real": lambda z: (float(np.real(z)), float(np.imag(z))),
           "sorted_by_abs": lambda z: float(np.abs(z))}[pairing]
    a = sorted(np.asarray(ev_t).ravel(), key=key)
    b = sorted(np.asarray(ev_ref).ravel(), key=key)
    if len(a) != len(b):
        return _defect_base("spectral_drift", object_ids,
                            defect_status="undefined_by_type",
                            pairing_rule=pairing,
                            payload={"reason": "eigenvalue counts differ -- "
                                     "emit a DimensionEvent, not a drift",
                                     "counts": [len(a), len(b)]})
    raw = float(np.linalg.norm(np.asarray(a) - np.asarray(b)))
    ref_scale = float(np.linalg.norm(np.asarray(b)))
    if ref_scale > 0.0:
        normalized, normalizer = raw / ref_scale, "norm(lambda_ref)"
    else:
        normalized, normalizer = raw, "unnormalized (zero reference spectrum)"
    status = "measured_zero" if normalized <= tolerance else "measured_nonzero"
    return _defect_base("spectral_drift", object_ids, defect_status=status,
                        raw_value=raw, norm_type="euclidean",
                        normalizer=normalizer, normalized_value=normalized,
                        tolerance=tolerance, pairing_rule=pairing)


def principal_angles(E_t, E_ref, aggregation="max", rank_rtol=1e-10):
    """Principal angles between subspaces spanned by the COLUMNS of E_t and
    E_ref (measurement helper -- embed in a spectral_signature or a defect
    payload; frame_drift_defect wraps this into a DefectRecord).

    Orthonormalizes via SVD and checks NUMERICAL rank: a rank-deficient
    input means the spanned subspace is smaller than its column count, and
    Gram-Schmidt noise directions would otherwise manufacture spurious
    small angles. rank_deficient=True in the output means treat this as a
    DimensionEvent / ill_conditioned, not as a drift measurement."""
    def orth(E):
        E = np.asarray(E, dtype=float)
        U, s, _ = np.linalg.svd(E, full_matrices=False)
        rank = int(np.sum(s > (s[0] * rank_rtol if s.size and s[0] > 0
                                else np.inf)))
        return U[:, :rank], rank, E.shape[1]
    Qt, rank_t, cols_t = orth(E_t)
    Qr, rank_r, cols_r = orth(E_ref)
    out = {"aggregation": aggregation,
           "column_counts": [cols_t, cols_r],
           "numerical_ranks": [rank_t, rank_r],
           "rank_deficient": bool(rank_t < cols_t or rank_r < cols_r)}
    if rank_t == 0 or rank_r == 0:
        out.update({"principal_angles": [], "principal_angle_agg": None})
        return out
    s = np.clip(np.linalg.svd(Qt.T @ Qr, compute_uv=False), -1.0, 1.0)
    angles = np.sort(np.arccos(s))
    agg = {"max": float(angles[-1]),
           "mean": float(np.mean(angles))}[aggregation]
    out.update({"principal_angles": [float(a) for a in angles],
                "principal_angle_agg": agg})
    return out


def frame_drift_defect(E_t, E_ref, object_ids=("E_t", "E_ref"),
                       aggregation="max", tolerance=1e-12, rank_rtol=1e-10):
    """Eigenframe/subspace drift as a DefectRecord (record producer):
    declared-aggregation principal angle between subspaces (paper 6.6).
    Rank deficiency or rank mismatch -> ill_conditioned / undefined_by_type
    instead of a fake angle."""
    pa = principal_angles(E_t, E_ref, aggregation=aggregation,
                          rank_rtol=rank_rtol)
    payload = {k: pa[k] for k in ("principal_angles", "column_counts",
                                  "numerical_ranks", "rank_deficient")}
    if pa["rank_deficient"]:
        return _defect_base("eigenframe_drift", object_ids,
                            defect_status="ill_conditioned",
                            aggregation=aggregation,
                            payload={**payload,
                                     "reason": "numerically rank-deficient "
                                     "input -- treat as DimensionEvent"})
    if pa["numerical_ranks"][0] != pa["numerical_ranks"][1]:
        return _defect_base("eigenframe_drift", object_ids,
                            defect_status="undefined_by_type",
                            aggregation=aggregation,
                            payload={**payload,
                                     "reason": "subspace dimensions differ "
                                     "-- emit a DimensionEvent"})
    raw = pa["principal_angle_agg"]
    status = "measured_zero" if raw <= tolerance else "measured_nonzero"
    return _defect_base("eigenframe_drift", object_ids, defect_status=status,
                        raw_value=raw, norm_type="principal_angle",
                        normalizer="radians (already scale-free)",
                        normalized_value=raw, tolerance=tolerance,
                        aggregation=aggregation, payload=payload)


# ------------------------------------------------------------- builders/IO

def _require(rec, fields, kind):
    missing = [f for f in fields if not rec.get(f)]
    if missing:
        raise ValueError(f"{kind} missing required fields: {missing}")


def capsule(**kw):
    kw.setdefault("capsule_id", f"cap-{uuid.uuid4().hex}")
    _require(kw, ("capsule_id", "value_ref", "object_kind", "carrier_space"),
             "capsule")
    gs = kw.get("grade_signature")
    if gs and "grade_norms" in gs and "norm_type" not in gs:
        raise ValueError("grade_signature.norm_type must be declared "
                         "(SKILL.md: load-bearing in indefinite signatures)")
    info = kw.get("information_signature")
    if info and not (info.get("random_variable") and info.get("estimator")):
        raise ValueError("no entropy without a declared random variable "
                         "AND estimator (paper 4.7)")
    return kw


def operation(**kw):
    kw.setdefault("operation_id", f"op-{uuid.uuid4().hex}")
    kw.setdefault("input_capsules", [])
    kw.setdefault("output_capsules", [])
    _require(kw, ("operation_id", "operation_name"), "operation")
    sev = kw.get("severity")
    if sev is not None and sev not in SEVERITIES:
        raise ValueError(f"severity must be one of {SEVERITIES}")
    return kw


def defect(**kw):
    kw.setdefault("defect_id", f"def-{uuid.uuid4().hex}")
    _require(kw, ("defect_id", "defect_type", "defect_status", "object_ids"),
             "defect")
    if kw["defect_status"] not in DEFECT_STATUSES:
        raise ValueError(f"defect_status must be one of {DEFECT_STATUSES}")
    return kw


def _json_default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not JSON serializable: {type(obj)}")


def emit(record, path):
    """Append one record as a JSONL line. Append-only: corrections are new
    records referencing the old id, never rewrites."""
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True, default=_json_default)
                + "\n")


# ---------------------------------------------------------------- selftest

def _selftest():
    ok = True

    def check(name, cond):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + name)
        ok = ok and cond

    # Paper Appendix B/C worked numbers (v1.1 corrected values). Field
    # names in this module follow paper sections 4/5/10, not the appendix
    # JSON as printed.
    norms = {0: 0.31, 1: 1.92, 2: 2.44, 3: 0.78, 4: 0.11}
    check("Appendix C entropy_before = 1.30 bits",
          abs(grade_entropy(norms, log_base=2) - 1.3039) < 5e-4)
    even = {k: v for k, v in norms.items() if k in (0, 2, 4)}
    check("Appendix C entropy_after = 0.14 bits",
          abs(grade_entropy(even, log_base=2) - 0.1382) < 5e-4)
    fr = grade_energy_fractions(norms)
    check("Appendix C loss_fraction = 0.41",
          abs(fr[1] + fr[3] - 0.4147) < 5e-4)

    # Commutator typing, status vocabulary, eps discipline.
    d = commutator_defect(np.diag([1.0, 2.0]), np.diag([3.0, 4.0]))
    check("commuting diagonals -> measured_zero",
          d["defect_status"] == "measured_zero")
    d = commutator_defect([[0, 1], [0, 0]], [[0, 0], [1, 0]])
    check("nilpotent pair -> measured_nonzero",
          d["defect_status"] == "measured_nonzero")
    d = commutator_defect(np.ones((2, 3)), np.ones((3, 2)))
    check("non-endomorphism -> undefined_by_type (never zero)",
          d["defect_status"] == "undefined_by_type"
          and "normalized_value" not in d)
    d = commutator_defect(np.zeros((2, 2)), np.eye(2), eps=1e-9)
    check("zero-norm operand with eps>0 -> still undefined_by_type",
          d["defect_status"] == "undefined_by_type")

    # Grade leakage on a Cl(3,1)-shaped coefficient vector (16 coeffs).
    coeffs = np.zeros(16)
    coeffs[0], coeffs[5] = 1.0, 2.0        # grade 0 and grade 2 content
    d = grade_leak(coeffs, allowed_grades={0, 2, 4})
    check("even-only multivector -> leak measured_zero",
          d["defect_status"] == "measured_zero")
    coeffs[1] = 1.0                        # add grade-1 content
    d = grade_leak(coeffs, allowed_grades=set(np.array([0, 2, 4])))
    check("odd content -> leak = 1/sqrt(6) (numpy grades accepted)",
          d["defect_status"] == "measured_nonzero"
          and abs(d["normalized_value"] - 1.0 / math.sqrt(6.0)) < 1e-12)
    json.dumps(d, default=_json_default)   # payload must serialize

    # Eigenvalue drift: tolerance discipline + count-mismatch handling.
    d = eigenvalue_drift_defect([1.0, 2.0 + 1e-15], [2.0, 1.0])
    check("permuted spectrum within tolerance -> measured_zero",
          d["defect_status"] == "measured_zero")
    d = eigenvalue_drift_defect([1.0, 2.0], [1.0, 2.0, 3.0])
    check("count mismatch -> undefined_by_type (DimensionEvent)",
          d["defect_status"] == "undefined_by_type")
    try:
        eigenvalue_drift_defect([1.0], [1.0], pairing="hungarian")
        check("unknown pairing raises ValueError", False)
    except ValueError:
        check("unknown pairing raises ValueError", True)

    # Principal angles: rotation recovered; rank deficiency flagged.
    theta = 0.3
    E1 = np.eye(4)[:, :2]
    R = np.eye(4)
    R[1, 1], R[1, 2], R[2, 1], R[2, 2] = (math.cos(theta), -math.sin(theta),
                                          math.sin(theta), math.cos(theta))
    d = frame_drift_defect(R @ E1, E1)
    check("rotated subspace -> max principal angle = theta",
          d["defect_status"] == "measured_nonzero"
          and abs(d["raw_value"] - theta) < 1e-9)
    bad = np.array([[1.0, 1.0], [1e-16, 0.0], [0.0, 0.0], [0.0, 0.0]])
    d = frame_drift_defect(bad, np.eye(4)[:, 1:3])
    check("rank-deficient input -> ill_conditioned, no fake angle",
          d["defect_status"] == "ill_conditioned"
          and "raw_value" not in d)

    # Builder guards.
    try:
        capsule(value_ref="x", object_kind="multivector",
                carrier_space="Cl(3,1)^1",
                grade_signature={"grade_norms": {}})
        check("capsule without norm_type rejected (presence-based)", False)
    except ValueError:
        check("capsule without norm_type rejected (presence-based)", True)
    try:
        defect(defect_type="noncommutativity", defect_status="broken",
               object_ids=["A"])
        check("bad defect_status rejected", False)
    except ValueError:
        check("bad defect_status rejected", True)

    # Schema validation of every record producer's output, including the
    # undefined/ill_conditioned branches (requires jsonschema).
    try:
        import jsonschema
        import os
        schema_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "..", "schemas",
                                   "defect_record.schema.json")
        with open(schema_path, encoding="utf-8") as f:
            schema = json.load(f)
        records = [
            commutator_defect(np.diag([1.0, 2.0]), np.diag([3.0, 4.0])),
            commutator_defect(np.ones((2, 3)), np.ones((3, 2))),
            commutator_defect(np.zeros((2, 2)), np.eye(2), eps=1e-9),
            grade_leak(coeffs, allowed_grades={0, 2, 4}),
            grade_leak(np.zeros(16), allowed_grades={0, 2, 4}),
            projection_loss([1.0, 2.0], [1.0, 0.0]),
            projection_loss([0.0, 0.0], [0.0, 0.0]),
            eigenvalue_drift_defect([1.0, 2.0], [2.0, 1.0]),
            eigenvalue_drift_defect([1.0, 2.0], [1.0, 2.0, 3.0]),
            frame_drift_defect(R @ E1, E1),
            frame_drift_defect(bad, np.eye(4)[:, 1:3]),
        ]
        errs = []
        for r in records:
            r2 = json.loads(json.dumps(r, default=_json_default))
            try:
                jsonschema.validate(r2, schema)
            except jsonschema.ValidationError as e:
                errs.append(f"{r.get('defect_type')}/"
                            f"{r.get('defect_status')}: {e.message}")
        check(f"all {len(records)} record-producer outputs validate "
              "against defect schema", not errs)
        for e in errs:
            print("  schema error:", e)
    except ImportError:
        print("SKIP schema validation (jsonschema not installed)")

    print("SELFTEST " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    print(__doc__)
