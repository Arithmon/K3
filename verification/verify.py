# -*- coding: utf-8 -*-
"""One-command verification of the certificates behind the paper.

    python3 verification/verify.py            # full verification
    python3 verification/verify.py --quick    # hashes only, no replay

TWO LEVELS, and the distinction is the point.

  REPLAY — the inexpensive certificates (under a second each) are RE-EXECUTED
    from the producers in `verification/producers/`, and we then check that
    every gate and every negative control is green and that the recorded
    outcome still carries its expected prefix. A regenerated artefact is NOT
    compared by hash: its provenance block records the commit it was built
    from, which changes with every commit, so hashing a regenerated file
    would test the version-control history rather than the mathematics.

  HASH — the expensive certificates (the bridge panel, the metric path, the
    face traversal) are not replayed here; their SHA-256 is checked against
    the recorded value. Reproducing them takes hours on a compute machine,
    which is not what a one-command check should ask of a reader.

The script is READ-ONLY on the tree. Replaying rewrites the artefacts, so the
original bytes are saved and RESTORED in a `finally` block: verifying must not
leave the repository dirty.

DEPENDENCIES. The producers use exact rational and interval arithmetic, and
two of them serialise the arithmetic backend into their provenance, so the
versions are pinned rather than merely recommended:

    python >= 3.10,  numpy,  sympy 1.14.0,  mpmath 1.3.0  (pinned)

The mpmath pin is enforced below: a different version would silently change
the last digits of every interval endpoint, and the certificates record the
version they were produced with.

Output is one line per certificate, a global verdict, and a nonzero exit code
if anything fails.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

ROOT = Path(__file__).resolve().parents[1]
CERTS = ROOT / "certificates"
PRODUCERS = Path(__file__).resolve().parent / "producers"

MPMATH_PIN = "1.3.0"

# (certificate, producer, expected outcome prefix) — replayed
REPLAY = [
    ("open_chart_theorem", "open_chart_theorem.py",
     "uniform_open_chart_theorem"),
    ("quantitative_atlas", "quantitative_atlas.py",
     "t2_fixed_k3_closed"),
    ("glue_obligations", "glue_obligations.py",
     "atlas_paper_glue_obligations_typed"),
    ("smoothness_and_transitions", "smoothness_and_transitions.py",
     "atlas_paper_three_nonzero_lemma"),
    # The sigma-floor correction. The paper states that writing it up
    # uncovered two defects in the results it was reporting, one of them the
    # sigma floor of the two charts (uniform radius 9.6e-10 -> 2.1e-12).
    # Without this certificate replayed, that claim of self-correction is
    # ASSERTED rather than CHECKABLE — which is precisely what a
    # one-command verifier exists to close.
    ("sigma_floor_correction", "sigma_floor_correction.py",
     "u1_sigma_floor_defect_confirmed"),
]

# (certificate, recorded SHA-256) — checked by hash, not replayed
HASHED = [
    ("bridge_atlas_panel",
     "8d01e64efb1bf52923f87de6b10d47b49e177e7285df294c0a678fbef0355150"),
    ("bridge_metric_path",
     "aa00dea0111577928869913ae4fa399655aac72577c65e6313e26c623baba316"),
    ("face_traversal_leaf",
     "35b5d0abd6f7a952f5970332ac9b013909485c7879fee3fed629b659278ecc54"),
]


def check_environment():
    """Refuse to run against an arithmetic backend the certificates disown."""
    try:
        import mpmath
    except ImportError:
        return "mpmath is not installed (required: %s)" % MPMATH_PIN
    if mpmath.__version__ != MPMATH_PIN:
        return ("mpmath %s is installed, but the certificates were produced "
                "with %s, which is serialised into their provenance; "
                "interval endpoints are not comparable across versions"
                % (mpmath.__version__, MPMATH_PIN))
    try:
        import numpy  # noqa: F401
        import sympy  # noqa: F401
    except ImportError as e:
        return "missing dependency: %s" % e.name
    return None


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check_json(path, prefix):
    d = json.loads(path.read_text(encoding="utf-8"))
    gp, gt = d.get("gates_passed"), d.get("gates_total")
    st = d.get("self_tests", {})
    ok_gates = gp is not None and gp == gt and gt > 0
    ok_neg = all(bool(v) for v in st.values()) if st else True
    ok_out = str(d.get("outcome", "")).startswith(prefix)
    detail = f"gates {gp}/{gt}"
    if st:
        detail += (f" · negative controls "
                   f"{sum(bool(v) for v in st.values())}/{len(st)}")
    return (ok_gates and ok_neg and ok_out), detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="hashes only, no replay")
    args = ap.parse_args()

    problem = check_environment()
    if problem:
        print("Certified K3 atlas — verification")
        print(f"\n  ENVIRONMENT: {problem}")
        print("\nVERDICT: FAIL — dependencies do not match the pin")
        return 2

    failures, lines = [], []

    if not args.quick:
        # Replaying rewrites the artefacts (provenance: timing, source
        # commit). Save the original bytes and restore them whatever happens.
        snapshot = {}
        for cert, _, _ in REPLAY:
            path = CERTS / f"{cert}.json"
            if path.exists():
                snapshot[path] = path.read_bytes()
        try:
            for cert, producer, prefix in REPLAY:
                sp = PRODUCERS / producer
                if not sp.exists():
                    failures.append(cert)
                    lines.append(f"  MISSING  {cert}  (producer absent)")
                    continue
                r = subprocess.run([sys.executable, str(sp)], cwd=str(ROOT),
                                   capture_output=True, text=True)
                if r.returncode != 0:
                    failures.append(cert)
                    lines.append(f"  FAIL     {cert}  "
                                 f"(replay: exit {r.returncode})")
                    continue
                ok, detail = check_json(CERTS / f"{cert}.json", prefix)
                lines.append(f"  {'PASS' if ok else 'FAIL'}     {cert}  "
                             f"[replay] {detail}")
                if not ok:
                    failures.append(cert)
        finally:
            for path, blob in snapshot.items():
                path.write_bytes(blob)
    else:
        lines.append("  (replay skipped: --quick)")

    for cert, expected in HASHED:
        path = CERTS / f"{cert}.json"
        if not path.exists():
            failures.append(cert)
            lines.append(f"  MISSING  {cert}")
            continue
        got = sha(path)
        ok = got == expected
        lines.append(f"  {'PASS' if ok else 'FAIL'}     {cert}  "
                     f"[hash] {got[:16]}…")
        if not ok:
            failures.append(cert)
            lines.append(f"           expected {expected[:16]}…")

    print("Certified K3 atlas — verification")
    print("\n".join(lines))
    if failures:
        print(f"\nVERDICT: FAIL — {len(failures)} certificate(s): "
              f"{', '.join(failures)}")
        return 1
    n = 0 if args.quick else len(REPLAY)
    print(f"\nVERDICT: PASS — {n} replayed, {len(HASHED)} hashed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
