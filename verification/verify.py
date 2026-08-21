# -*- coding: utf-8 -*-
"""One-command verification of the certificates behind the paper.

    python3 verification/verify.py            # full verification
    python3 verification/verify.py --quick    # hashes only, no replay

THREE LEVELS, and the distinction is the point.

  REPLAY — the inexpensive certificates (under a second each) are RE-EXECUTED
    from the producers in `verification/producers/`, and we then check that
    every gate and every negative control is green and that the recorded
    outcome still carries its expected prefix. A regenerated artefact is NOT
    compared by hash: its provenance block records the commit it was built
    from, which changes with every commit, so hashing a regenerated file
    would test the version-control history rather than the mathematics.

  RECOMPUTE — the coverage certificate (71.8 M boxes, about 70 s on four
  cores) is RECOMPUTED and its counters compared one by one with the shipped
  file. It predates the gates/self_tests convention, so reading it would only
  re-read a verdict; recomputing it would redden if a single counter moved.
  This is the computation carrying the pivot floor 4.8, so it is the one that
  most deserves to be reproduced rather than trusted.

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

Output is one line per checked certificate, a global verdict, and a nonzero
exit code if anything fails. The full run takes about two minutes, almost all
of it in the recomputed coverage certificate; `--quick` skips both replay and
recomputation and checks hashes only.
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
    # The four design certificates are deliberately NOT replayed. They do
    # have gates, self-tests and producers, and they pass; but shipping their
    # producers pulls twenty-three further artefacts into this repository,
    # through everything those producers read — contract amendments,
    # preregistrations, the internal chain those design documents rest on.
    # The repository would go from 14 certificates to 37, most of them
    # carrying no claim the paper makes. What the command covers is therefore
    # 9 of the 14 files shipped, and Appendix F says so without rounding.
]

# (certificate, producer, argv, compared fields) — RECOMPUTED, then compared
# field by field with the shipped certificate. The arguments are those of the
# published run; the box budget in particular is NOT optional: at the default
# (200 M) gauge 1 stops on `budget_abort` and the verdict falls to FAILED on 16
# unresolved boxes, while the other five gauges reproduce identically either
# way. `seconds` and `n_cores` are excluded from the comparison, being machine
# dependent; everything else is deterministic and must agree counter for counter.
RECOMPUTE = [
    ("atlas_coverage", "atlas_coverage.py",
     ["0.6", "4", "1e-3", "1000"],
     ("verdict", "tau", "w_min", "minor_floor_8tau", "total_boxes",
      "per_gauge", "statement", "domain", "arithmetic", "mu")),
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


def check_shipped(blob, prefix):
    """Is the SHIPPED certificate green? Replay writes and then RESTORES the
    original bytes, so a repository shipping a RED artefact still passed:
    only the regenerated file was examined, never the one a reader opens."""
    d = json.loads(blob.decode("utf-8"))
    gp, gt = d.get("gates_passed"), d.get("gates_total")
    st = d.get("self_tests", {})
    field = "outcome" if "outcome" in d else "issue"
    return (gp is not None and gp == gt and gt > 0
            and all(bool(v) for v in st.values())
            and bool(prefix) and str(d.get(field, "")).startswith(prefix))


def check_json(path, prefix):
    d = json.loads(path.read_text(encoding="utf-8"))
    gp, gt = d.get("gates_passed"), d.get("gates_total")
    st = d.get("self_tests", {})
    ok_gates = gp is not None and gp == gt and gt > 0
    ok_neg = all(bool(v) for v in st.values()) if st else True
    # Two naming conventions coexist: recent certificates serialise
    # `outcome`, the design ones serialise `issue`. Accept both EXPLICITLY,
    # rather than letting the second pass for an empty outcome.
    field = "outcome" if "outcome" in d else "issue"
    ok_out = bool(prefix) and str(d.get(field, "")).startswith(prefix)
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
        for cert in ([c for c, _, _ in REPLAY]
                     + [c for c, _, _, _ in RECOMPUTE]):
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
                shipped_ok = check_shipped(snapshot[CERTS / f"{cert}.json"],
                                           prefix)
                lines.append(f"  {'PASS' if ok and shipped_ok else 'FAIL'}"
                             f"     {cert}  [replay] {detail}")
                if not shipped_ok:
                    lines.append("           but the SHIPPED certificate is "
                                 "red — replay regenerates it and restores "
                                 "the original")
                if not ok or not shipped_ok:
                    failures.append(cert)
            for cert, producer, argv, fields in RECOMPUTE:
                sp = PRODUCERS / producer
                path = CERTS / f"{cert}.json"
                if not sp.exists() or not path.exists():
                    failures.append(cert)
                    lines.append(f"  MISSING  {cert}  (producer or certificate absent)")
                    continue
                shipped = json.loads(snapshot[path].decode("utf-8"))
                r = subprocess.run([sys.executable, str(sp), *argv],
                                   cwd=str(ROOT), capture_output=True, text=True)
                if r.returncode != 0:
                    failures.append(cert)
                    lines.append(f"  FAIL     {cert}  (recompute: exit {r.returncode})")
                    continue
                got = json.loads(path.read_text(encoding="utf-8"))
                # Without this check, a producer writing somewhere OTHER than
                # the file being compared makes the diff trivially empty: the
                # shipped certificate is compared with itself, and a mutated
                # counter passes. `seconds` changes on every run, so if it is
                # unchanged the recomputation did not land in this file.
                rewritten = got.get("seconds") != shipped.get("seconds")
                diff = [f for f in fields if got.get(f) != shipped.get(f)]
                ok = (rewritten and not diff
                      and got.get("verdict") == "CERTIFIED")
                if not rewritten:
                    lines.append(f"           recomputation did not rewrite "
                                 f"{path.name} — nothing was compared")
                detail = (f"{len(fields)} fields identical, "
                          f"{shipped.get('total_boxes', 0)} boxes recomputed"
                          if ok else f"fields differing: {diff}")
                lines.append(f"  {'PASS' if ok else 'FAIL'}     {cert}  "
                             f"[recompute] {detail}")
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
    m = 0 if args.quick else len(RECOMPUTE)
    print(f"\nVERDICT: PASS — {n} replayed, {m} recomputed, "
          f"{len(HASHED)} hashed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
