# -*- coding: utf-8 -*-
"""One-command verification of the certificates behind the paper.

    python3 verification/verify.py            # full verification
    python3 verification/verify.py --quick    # hashes only, no replay

THREE LEVELS, and the distinction is the point.

  REPLAY — the inexpensive certificates (under a second each) are re-executed
    from the producers in `verification/producers/`; every check and every
    negative control must be green and the recorded outcome must be exactly
    the expected one. A regenerated file is not compared by hash: its
    provenance records the commit it was built from.

  RECOMPUTE — the coverage certificate (71.8 M boxes, about 70 s on four
    cores) is recomputed and its counters compared one by one with the
    shipped file; it is the computation carrying the pivot floor 4.8.

  HASH — the expensive certificates (the bridge panel, the metric path, the
    face traversal) are checked by SHA-256 against the recorded value; their
    producers ship in `src/k3_atlas` and take hours to re-run.

Also checked: every `upstream` pin the producers wrote, against the shipped
bytes; the results index against what this command does; the hash table of
Appendix F against the files, in the markdown, the LaTeX and the PDF.

The script is read-only on the tree: replaying rewrites the certificates, so
the original bytes are saved and restored.

DEPENDENCIES. The reference environment is pinned in `requirements.lock`
(python 3.14.4, numpy 2.5.1, sympy 1.14.0, mpmath 1.3.0, among others). This
command needs numpy, sympy and mpmath; only the mpmath pin is
ENFORCED here, the others must merely import, so that the command also runs
outside the pinned environment (it has been run under python 3.12). The
producers use exact rational and interval arithmetic, and two of them
serialise the arithmetic backend into their provenance: a different mpmath
would silently change the last digits of every interval endpoint, and the
certificates record the version they were produced with.

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
import re
import subprocess
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

ROOT = Path(__file__).resolve().parents[1]
CERTS = ROOT / "certificates"
PRODUCERS = Path(__file__).resolve().parent / "producers"

MPMATH_PIN = "1.3.0"

# (certificate, producer, EXACT expected outcome — not a prefix) — replayed
REPLAY = [
    ("open_chart_theorem", "open_chart_theorem.py",
     "uniform_open_chart_theorem_certified"),
    ("quantitative_atlas", "quantitative_atlas.py",
     "t2_fixed_k3_closed_quantitative_existence_level"),
    ("glue_obligations", "glue_obligations.py",
     "atlas_paper_glue_obligations_typed_and_transition_generators_derived_word_length_le_4"),
    ("smoothness_and_transitions", "smoothness_and_transitions.py",
     "atlas_paper_three_nonzero_lemma_carries_smoothness_and_coverage_pivot_degree3_invariant_transitions_explicit"),
    # the correction reported in Appendix B, replayed so that it is checked
    # rather than asserted
    ("sigma_floor_correction", "sigma_floor_correction.py",
     "u1_sigma_floor_defect_confirmed_with_witness_radius_corrected_9p6e10_to_2p1e12_theorem_survives"),
    # The four design records are not replayed: their producers would pull
    # twenty-three further files into the repository (Appendix F).
]

# (certificate, producer, argv, compared fields) — recomputed, then compared
# field by field. The arguments are those of the published run (the box
# budget is not optional: at the default, gauge 1 stops on `budget_abort`).
# `seconds` and `n_cores` are machine dependent and excluded.
RECOMPUTE = [
    ("atlas_coverage", "atlas_coverage.py",
     ["0.6", "4", "1e-3", "1000"],
     ("verdict", "tau", "w_min", "minor_floor_8tau", "total_boxes",
      "per_gauge", "statement", "domain", "arithmetic", "mu")),
]

# (certificate, recorded SHA-256) — checked by hash, not replayed
HASHED = [
    ("bridge_atlas_panel",
     "58c06da48412c8b1dca175f7f7b44156b590d7f71be80b4831cdd4f4f1c5081c"),
    ("bridge_metric_path",
     "397c19b105819f79f93016fda4045e35b969375d8bfc2d880534c8096cc076be"),
    ("face_traversal_leaf",
     "9088e59844c7e0dd398cd5daa66ce11a1ddece7d4aca8fe8c07b47d08b9ae5c1"),
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
    """Is the SHIPPED certificate green (not only the regenerated one)?"""
    d = json.loads(blob.decode("utf-8"))
    gp, gt = d.get("checks_passed"), d.get("checks_total")
    st = d.get("perturbation_tests", {})
    field = "outcome" if "outcome" in d else "issue"
    return (gp is not None and gp == gt and gt > 0
            and all(bool(v) for v in st.values())
            and bool(prefix) and str(d.get(field, "")) == prefix)


def check_json(path, prefix):
    d = json.loads(path.read_text(encoding="utf-8"))
    gp, gt = d.get("checks_passed"), d.get("checks_total")
    st = d.get("perturbation_tests", {})
    ok_gates = gp is not None and gp == gt and gt > 0
    ok_neg = all(bool(v) for v in st.values()) if st else True
    # recent certificates serialise `outcome`, the design records `issue`
    field = "outcome" if "outcome" in d else "issue"
    ok_out = bool(prefix) and str(d.get(field, "")) == prefix
    detail = f"checks {gp}/{gt}"
    if st:
        detail += (f" · negative controls "
                   f"{sum(bool(v) for v in st.values())}/{len(st)}")
    return (ok_gates and ok_neg and ok_out), detail


# Appendix F publishes the SHA-256 prefix of each hash-verified certificate;
# a reader who runs `sha256sum` compares against that table. It is checked
# against the files in the markdown, the LaTeX and the deposited PDF.
PAPER_HASH_TABLE = {
    "bridge panel": "bridge_atlas_panel",
    "metric path": "bridge_metric_path",
    "face crossing": "face_traversal_leaf",
}
PAPER_FILES = ("paper/certified_k3_atlas.md",
               "paper/latex/certified_k3_atlas_body.tex")
PAPER_PDF = "paper/latex/certified_k3_atlas.pdf"
_PAPER_SEEN = []


def _pdf_text(path):
    import zlib
    raw = path.read_bytes()
    out = b""
    for m in re.finditer(rb"stream\r?\n(.*?)endstream", raw, re.S):
        try:
            out += zlib.decompress(m.group(1))
        except Exception:
            pass
    return out.decode("latin-1", "replace")


def check_paper_hashes():
    """Lines where the paper's published prefix differs from the shipped file.
    Every label must be found in every source, exactly once in the PDF."""
    bad, seen = [], []
    for rel in PAPER_FILES:
        f = ROOT / rel
        if not f.exists():
            bad.append((rel, 0, "(the whole table)", "file missing", ""))
            continue
        found = set()
        for n, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            for label, cert in PAPER_HASH_TABLE.items():
                if label not in line:
                    continue
                m = re.search(r"[0-9a-f]{16}", line)
                if not m:
                    continue
                found.add(label)
                path = CERTS / f"{cert}.json"
                want = sha(path)[:16] if path.exists() else "(missing)"
                if m.group(0) != want:
                    bad.append((rel, n, label, m.group(0), want))
        for label in PAPER_HASH_TABLE:
            if label not in found:
                bad.append((rel, 0, label, "no such row", "a row naming it"))
        seen.append(len(found))
        _PAPER_SEEN.append(len(found))
    pdf = ROOT / PAPER_PDF
    if pdf.exists():
        text = _pdf_text(pdf)
        for label, cert in PAPER_HASH_TABLE.items():
            path = CERTS / f"{cert}.json"
            want = sha(path)[:16] if path.exists() else "(missing)"
            if want in text:
                _PAPER_SEEN.append(1)
                if text.count(want) != 1:
                    bad.append((PAPER_PDF, 0, label,
                                f"{text.count(want)} occurrences", want))
            if want not in text:
                stale = [h for h in re.findall(r"[0-9a-f]{16}", text)
                         if h != want and h.count(h[0]) < 8]
                bad.append((PAPER_PDF, 0, label,
                            "absent" if not stale else "a different value",
                            want))
    else:
        bad.append((PAPER_PDF, 0, "(the deposited file)", "missing", ""))
    return bad


# `docs/RESULTS_INDEX.md` says, per certificate, how this command checks it.
# That column is a claim about this file, and is checked against it.
def check_results_index():
    """Rows where the index disagrees with what this command actually does."""
    import re
    f = ROOT / "docs" / "RESULTS_INDEX.md"
    if not f.exists():
        return [("docs/RESULTS_INDEX.md", 0, "(missing)", "", "")]
    actual = {c: "replayed" for c, _, _ in REPLAY}
    actual.update({c: "recomputed" for c, _, _, _ in RECOMPUTE})
    actual.update({c: "hashed" for c, _ in HASHED})
    bad, seen = [], set()
    for n, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
        cells = [c.strip() for c in line.split("|")]
        if len(cells) < 6 or not cells[3].startswith("`"):
            continue
        cert = cells[3].strip("`").removesuffix(".json")
        claimed = cells[4]
        seen.add(cert)
        want = actual.get(cert, "unverified")
        if claimed != want:
            bad.append(("docs/RESULTS_INDEX.md", n, cert, claimed, want))
        if not (CERTS / f"{cert}.json").exists():
            bad.append(("docs/RESULTS_INDEX.md", n, cert, "listed",
                        "no such certificate"))
    for cert in sorted(set(actual) - seen):
        bad.append(("docs/RESULTS_INDEX.md", 0, cert, "absent from the index",
                    actual[cert]))
    # every shipped certificate has a row, checked or not
    shipped = {p.stem for p in CERTS.glob("*.json")}
    for cert in sorted(shipped - seen - set(actual)):
        bad.append(("docs/RESULTS_INDEX.md", 0, cert, "absent from the index",
                    "a row (every shipped certificate has one)"))
    return bad


# Which shipped file an `upstream` key names (the keys are the producers' own).
PIN_TARGETS = {
    "d3_coverage_legacy": "certificates/atlas_coverage.json",
    "f1prime_u0": "certificates/uniform_chart_lemma.json",
    "amendment_b": "certificates/closure_skeleton.json",
    "u1": "certificates/open_chart_theorem.json",
    "u1_live": "certificates/open_chart_theorem.json",
    "u1_certificate": "certificates/open_chart_theorem.json",
    "t2": "certificates/quantitative_atlas.json",
    "bridge_panel": "certificates/bridge_atlas_panel.json",
    "generators": "certificates/glue_obligations.json",
    "f9_p0": "certificates/exact_transitions.json",
    "closeout_k_regional": "certificates/gluing_contract.json",
    "dyadic_cover": "src/k3_atlas/data/dyadic_cover.json",
}
# Pins that do not match the shipped bytes, with what was measured about the
# difference. A design record pins the bytes it was written from; its producer
# is not shipped, and some of its inputs were regenerated or translated after
# it. Each entry is printed on every run; any other mismatch fails.
DECLARED_STALE = {
    ("uniform_chart_lemma", "dyadic_cover"):
        "the input was regenerated after the record; compared leaf by leaf, "
        "0 numeric or boolean leaves differ, 3 prose leaves were translated "
        "and the counter block was renamed to `checks`",
    ("uniform_chart_lemma", "d3_coverage_legacy"):
        "the coverage certificate was recomputed after the record (this "
        "command recomputes it field by field); the pinned bytes are gone",
    ("uniform_chart_lemma", "amendment_b"):
        "closure_skeleton was translated after the record; the pin names the "
        "untranslated bytes",
    ("exact_transitions", "amendment_b"):
        "closure_skeleton was translated after the record; the pin names the "
        "untranslated bytes",
    ("exact_transitions", "f1prime_u0"):
        "uniform_chart_lemma was translated after the record; the pinned bytes "
        "are no longer available",
}
_GENERIC_KEY = re.compile(r"sha|hash|digest|path|file|npz|json", re.I)


def check_upstream_pins():
    """(cert, key, published, shipped, status) for every pin the table maps.

    status: 'ok', 'declared' (in DECLARED_STALE), or 'STALE' (a failure)."""
    rows = []
    for path in sorted(CERTS.glob("*.json")):
        d = json.loads(path.read_text(encoding="utf-8"))
        block = d.get("upstream", d.get("upstream_chain"))
        if not isinstance(block, (dict, list)):
            continue

        def walk(o, keys):
            if isinstance(o, dict):
                for k, v in o.items():
                    walk(v, keys + [str(k)])
            elif isinstance(o, list):
                for v in o:
                    walk(v, keys)
            elif isinstance(o, str):
                for h in re.findall(r"[0-9a-f]{64}", o):
                    named = [k for k in keys if not _GENERIC_KEY.search(k)]
                    name = named[-1] if named else (keys[-1] if keys else "")
                    if name not in PIN_TARGETS:
                        continue
                    target = ROOT / PIN_TARGETS[name]
                    got = sha(target) if target.exists() else "(missing)"
                    if got == h:
                        status = "ok"
                    elif (path.stem, name) in DECLARED_STALE:
                        status = "declared"
                    else:
                        status = "STALE"
                    rows.append((path.stem, name, h, got, status))
        walk(block, [])
    return rows


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
                # the recomputation must have rewritten THIS file
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

    pins = check_upstream_pins()
    stale = [r for r in pins if r[4] == "STALE"]
    for cert, name, h, got, _ in stale:
        lines.append(f"  FAIL     upstream pins  {cert}.{name}: pinned {h[:16]}…, "
                     f"shipped {got[:16]}…")
    for cert, name, h, got, st in pins:
        if st == "declared":
            lines.append(f"  declared {cert}.{name}: pinned {h[:16]}…, shipped "
                         f"{got[:16]}… — {DECLARED_STALE[(cert, name)]}")
    if stale:
        failures.append("upstream pins")
    else:
        lines.append(f"  PASS     upstream pins  [{len(pins)} pins checked against the "
                     f"shipped bytes, {sum(1 for r in pins if r[4] == 'declared')} declared stale]")

    index_bad = check_results_index()
    for rel, n, cert, claimed, want in index_bad:
        lines.append(f"  FAIL     results index  [{rel}:{n}] {cert}: "
                     f"says {claimed!r}, this command does {want!r}")
    if index_bad:
        failures.append("results index")
    else:
        lines.append("  PASS     results index  "
                     "[every row agrees with what this command does]")

    paper_bad = check_paper_hashes()
    for rel, n, label, got, want in paper_bad:
        lines.append(f"  FAIL     paper table  [{rel}:{n}] {label}: "
                     f"published {got}…, shipped {want}…")
    if paper_bad:
        failures.append("paper hash table")
    else:
        lines.append(f"  PASS     paper hash table  "
                     f"[{sum(_PAPER_SEEN)} rows read across "
                     f"{len(PAPER_FILES) + 1} files, all matching]")

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
