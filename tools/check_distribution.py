#!/usr/bin/env python3
"""Check that this repository is fit for distribution.

A published repository is read by people who were not in the room. Four
things must hold, and each is checked against the files as they are, never
against an intention:

  1. ENGLISH ONLY. Prose in .py, .md, .tex and .sh, command-line output and
     JSON keys are English. A French sentence here is a sentence the reader
     cannot use.
  2. NO DEVELOPMENT VOCABULARY. Milestone identifiers, internal phase names,
     machine names and review numbers say when a thing was built, not what it
     computes. They do not belong in a published result.

     Two things are NOT development vocabulary, and an early version of this
     rule wrongly flagged both: the AUTHOR'S NAME, and the DISCLOSURE OF
     ASSISTANCE FROM LANGUAGE MODELS. The second is required academic
     content, not jargon; a rule that removed it would push the paper
     towards concealing something it must state. They are exempt by name
     below, and the exemption is narrow: it covers the author line and the
     acknowledgements, not the code.
  3. NO UNRESOLVED REFERENCE. Every path a file names exists here; no file
     points into a workspace the reader does not have.
  4. NO OUTSIDE DEPENDENCY. Every import resolves to the standard library,
     to a pinned third-party package, or to this repository.

The check exits nonzero on any violation and prints one line per finding.
Its exemptions are few, named in this file, and PRINTED on every run with the
count they cover: an exemption nobody sees is a hole, one printed each time
is a declared limit. None of them carves our own vocabulary out of the rule
that names it.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TEXT_SUFFIXES = {".py", ".md", ".tex", ".sh", ".toml", ".cff", ".txt", ".lock"}
SKIP_DIRS = {".git", "__pycache__", "build", ".venv", "node_modules"}

# Words that mark development history rather than mathematics. The list is
# derived from what the workspace actually used; a repository that never
# used them is unaffected by it.
DEVELOPMENT_VOCABULARY = [
    r"\bk3_cap_\w+", r"\bb1e2iii\b", r"\bc1[0-9]{2}[a-z]?\b", r"\bd5\.\d\b",
    r"\bf1prime\b", r"\brface\b", r"\blot\b", r"\bledger\b", r"\bscout\b",
    r"\bfront\b", r"\bfreeze\b", r"\bgate[sd]?\b", r"\bcodex\b", r"\bgrok\b",
    r"\bkimi\b", r"\brevue\b",
    r"\b\d+(?:e|ᵉ|st|nd|rd|th)\s+revue\b", r"\bgift\b", r"\bk3-cap\b",
]

# Common French function words. Their presence in a line of prose is what
# distinguishes a French sentence from an English one carrying an accent.
#
# TWO of them make a sentence -- but only in a line long enough to hold two.
# A short French clause carries one, and it used to pass: `"phase B1 du
# contrat ... amendee ..."` shipped inside a certificate under that rule, and
# so did forty-three others. What separates them from English prose is not a
# second function word but an ACCENTED FRENCH WORD, so one marker plus one
# accented word now counts as French too. The accented-word rule needs its own
# exceptions, since English mathematics borrows accented names.
FRENCH_MARKERS = [
    r"\ble\b", r"\bla\b", r"\bles\b", r"\bdes\b", r"\bdu\b", r"\bune\b",
    r"\bun\b", r"\bqui\b", r"\bque\b", r"\bpour\b", r"\bavec\b", r"\bsur\b",
    r"\bdans\b", r"\bpas\b", r"\bplus\b", r"\best\b", r"\bsont\b",
    r"\bcette\b", r"\bce\b", r"\bnous\b", r"\bnon\b", r"\bmais\b",
    r"\bdonc\b", r"\bsans\b", r"\bleur\b", r"\bd'un\b", r"\bl'\w+",
]

# Accented words that belong to English mathematical prose: proper names and
# loanwords. The list is narrow on purpose -- anything not here that carries a
# French accent is French until shown otherwise.
ACCENTED_OK = {
    "kahler", "kähler", "poincare", "poincaré", "mobius", "möbius",
    "fourniere", "fournière", "hormander", "hörmander", "cech", "čech",
    "erdos", "erdős", "godel", "gödel", "resume", "résumé", "role", "rôle",
    "naive", "naïve", "hyperkahler", "hyperkähler", "plucker", "plücker",
    "ampere", "ampère", "kaehler", "bezout", "bézout", "frechet", "fréchet",
    "lebesgue", "levy", "lévy", "schrodinger", "schrödinger", "a", "e",
}
# The Latin-1 letter ranges, MINUS the two symbols they contain: U+00D7 is
# the multiplication sign and U+00F7 the division sign. Leaving them in made
# every `3x6` in a formula look like a French word.
ACCENTED = re.compile(
    r"\b\w*[\u00c0-\u00d6\u00d8-\u00dd\u00e0-\u00f6\u00f8-\u00ff"
    r"\u0152\u0153]\w*\b", re.I)


def french_score(s: str, fr) -> int:
    """How many INDEPENDENT signs of French this text carries.

    A distinct function word counts one. An accented word that is not a
    known English loan counts one. Two signs make French.
    """
    n = len(set(m.group(0).lower() for m in fr.finditer(s)))
    # An accented word that is not a known English loan is worth two on its
    # own. It had to be: a trailing French participle in an otherwise English
    # sentence scored one and shipped, because it carried no French function
    # word at all. The whitelist above is what keeps borrowed names from
    # tripping this.
    if any(w.lower() not in ACCENTED_OK for w in ACCENTED.findall(s)):
        n += 2
    return n


# Files that CANNOT be translated, and the reason, derived rather than
# asserted: `witness_manifest.json` is compared field for field against a
# manifest embedded inside the frozen witness archive, and the loader refuses
# the witness when the two differ. Translating the sidecar breaks that check --
# measured, not feared: the chain stopped on it. Rewriting the frozen archive
# instead would change its SHA-256, which the data files downstream record as
# `witness_sha256`.
# The prose stays as recorded, and this exemption says why rather than hiding
# it.
MIRRORED_IN_A_FROZEN_ARTEFACT = {
    "src/k3_atlas/data/witness_manifest.json":
        "mirrored inside k3_closedform_witness_kahler_v2.npz; the loader "
        "compares the two and refuses the witness if they differ",
}

# The inputs of the expensive chain, as shipped. Their bytes are pinned by the
# `upstream` blocks of the hash-verified certificates, so a word changed there
# changes a hash the paper publishes: what they record -- the run that
# produced each row, the paths of the retraction registry, a sentence left in
# French -- stays as recorded, and is COUNTED here rather than hidden.
FROZEN_INPUTS = "src/k3_atlas/data"

# Prose leaves of a hash-verified certificate that are known to be wrong
# (half French), by exact path. The certificate is pinned by two other
# hash-verified certificates, so a hand edit would break their pins; the
# producer's strings are corrected, and the next full replay of the chain
# removes these leaves. Until then they are declared, not silently passed.
HASH_FROZEN_PROSE = {
    ("certificates/bridge_atlas_panel.json", "/not_paid_here[2]"),
    ("certificates/bridge_atlas_panel.json", "/not_paid_here[3]"),
}

# Keys under which a certificate RECORDS where it came from. The name of a
# file that was read is provenance, and provenance is not rewritten to look
# tidier; development vocabulary under these keys is allowed and counted.
PROVENANCE_KEY = re.compile(
    r"(^upstream|provenance|^amends$|lineage|ledger|notes_read|_frozen$|"
    r"_ref$|built_from_head|self_sha256)", re.I)

# What the exemptions above let through on this run, printed at the end.
DECLARED: dict[str, int] = {
    "prose leaves of a hash-verified certificate": 0,
    "prose leaves of the frozen chain inputs": 0,
    "development vocabulary in the frozen chain inputs": 0,
    "development vocabulary under lineage keys of a certificate": 0,
}

STDLIB_OK = set(sys.stdlib_module_names)


def files():
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or p.suffix not in TEXT_SUFFIXES:
            continue
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        yield p


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT))


def check_english(findings):
    """Rule 1 — prose is English."""
    fr = re.compile("|".join(FRENCH_MARKERS), re.I)
    for p in files():
        for n, line in enumerate(p.read_text(encoding="utf-8",
                                             errors="replace").splitlines(), 1):
            s = line.strip()
            if not s or len(s) < 25:
                continue
            # Two independent signs of French in one line make prose, not a
            # stray token: one such sign can be a variable name, two of them
            # in a row make a sentence.
            if french_score(s, fr) >= 2:
                findings.append(("english", rel(p), n, s[:100]))


# Lines that legitimately carry a name a reader needs: authorship and the
# disclosure of language-model assistance. Matched on the LINE, so the
# exemption cannot spread into the code.
NAME_EXEMPT = re.compile(
    r"(author|thanks|acknowledg|Anthropic|OpenAI|Claude|correspondence|"
    r"^\*\*[A-Z])", re.I)
# On an exempt line only these NAMES pass; every other forbidden word on the
# same line is still a finding. The exemption used to cover the whole line,
# and a line opening with a name could carry six forbidden words behind it.
NAMES_ALLOWED_ON_EXEMPT_LINES = {"codex", "grok", "kimi"}


def check_vocabulary(findings):
    """Rule 2 — no development vocabulary."""
    voc = re.compile("|".join(DEVELOPMENT_VOCABULARY), re.I)
    for p in files():
        if rel(p) == "tools/check_distribution.py":
            continue        # the rule names the words it forbids
        for n, line in enumerate(p.read_text(encoding="utf-8",
                                             errors="replace").splitlines(), 1):
            exempt = bool(NAME_EXEMPT.search(line.strip()))
            for m in voc.finditer(line):
                if exempt and m.group(0).lower() in NAMES_ALLOWED_ON_EXEMPT_LINES:
                    continue
                findings.append(("vocabulary", rel(p), n, m.group(0)))


def check_json_keys(findings):
    """Rule 1 again, on the machine-readable surface.

    Keys AND string values: a certificate whose key is English and whose
    explanation is French is not readable by the audience it is published
    for.
    """
    fr = re.compile("|".join(FRENCH_MARKERS), re.I)

    def walk(o, path, p):
        if isinstance(o, dict):
            for k, v in o.items():
                if not str(k).isascii() or fr.search(str(k)):
                    findings.append(("json-key", rel(p), 0, f"{path}/{k}"))
                walk(v, f"{path}/{k}", p)
        elif isinstance(o, str) and len(o) > 20:
            if french_score(o, fr) >= 2:
                if (rel(p), path) in HASH_FROZEN_PROSE:
                    DECLARED["prose leaves of a hash-verified certificate"] += 1
                elif rel(p).startswith(FROZEN_INPUTS):
                    DECLARED["prose leaves of the frozen chain inputs"] += 1
                else:
                    findings.append(("json-text", rel(p), 0,
                                     f"{path} = {o[:70]}"))
        elif isinstance(o, list):
            # EVERY item. The first fifty used to be enough, and a French
            # sentence at index 100 of a 200-item list shipped under a green
            # verdict.
            for i, v in enumerate(o):
                walk(v, f"{path}[{i}]", p)

    for p in sorted(ROOT.rglob("*.json")):
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        if rel(p) in MIRRORED_IN_A_FROZEN_ARTEFACT:
            continue
        try:
            walk(json.loads(p.read_text(encoding="utf-8")), "", p)
        except (ValueError, OSError):
            findings.append(("json-key", rel(p), 0, "unreadable"))


def report_exemptions():
    """Say out loud what is exempt, and why. An exemption nobody sees is a
    hole; an exemption printed on every run is a declared limit."""
    for rel_path, why in sorted(MIRRORED_IN_A_FROZEN_ARTEFACT.items()):
        print(f"  exempt      {rel_path}: {why}")
    for rel_path, path in sorted(HASH_FROZEN_PROSE):
        print(f"  exempt      {rel_path}{path}: hash-verified and pinned twice; "
              "corrected in the producer, replaced at the next full replay")
    for what, n in DECLARED.items():
        print(f"  declared    {n:4d}  {what}")


def check_json_vocabulary(findings):
    """Rule 2, on the machine-readable surface.

    Development vocabulary in a certificate is a finding, except under the
    keys that RECORD provenance (counted), and in the frozen inputs of the
    expensive chain, whose bytes are pinned by published hashes (counted).
    """
    voc = re.compile("|".join(DEVELOPMENT_VOCABULARY), re.I)

    def record(p, path, word, keys):
        r = rel(p)
        if r.startswith(FROZEN_INPUTS):
            DECLARED["development vocabulary in the frozen chain inputs"] += 1
        elif any(PROVENANCE_KEY.search(k) for k in keys):
            DECLARED["development vocabulary under lineage keys of a certificate"] += 1
        else:
            findings.append(("json-vocabulary", r, 0, f"{path}: {word}"))

    def walk(o, path, p, keys):
        if isinstance(o, dict):
            for k, v in o.items():
                for m in voc.finditer(str(k)):
                    record(p, f"{path}/{k}", m.group(0), keys)
                walk(v, f"{path}/{k}", p, keys + [str(k)])
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}[{i}]", p, keys)
        elif isinstance(o, str):
            for m in voc.finditer(o):
                record(p, path, m.group(0), keys)

    for p in sorted(ROOT.rglob("*.json")):
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        if rel(p) in MIRRORED_IN_A_FROZEN_ARTEFACT:
            continue
        try:
            walk(json.loads(p.read_text(encoding="utf-8")), "", p, [])
        except (ValueError, OSError):
            pass        # rule 1 already reports an unreadable file


PATH_SUFFIXES = TEXT_SUFFIXES | {".log", ".aux", ".out", ".toc", ".json", ".cff"}


def path_files():
    """Rule 3 reads MORE than the prose rules: a LaTeX log is not prose, but
    it carried the build machine's home directory into a deposit once."""
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or p.suffix not in PATH_SUFFIXES:
            continue
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        yield p


def check_paths(findings):
    """Rule 3 — no path points outside this repository."""
    pat = re.compile(r"""["'`]((?:\.\./|/home/|/Users/|C:\\)[^"'`\n]{2,120})""")
    for p in path_files():
        for n, line in enumerate(p.read_text(encoding="utf-8",
                                             errors="replace").splitlines(), 1):
            for m in pat.finditer(line):
                target = m.group(1)
                # A relative path that RESOLVES inside the repository is not
                # an unresolved reference. Flagging it would teach the reader
                # to ignore this rule.
                if target.startswith("../"):
                    if (p.parent / target).resolve().exists():
                        continue
                findings.append(("path", rel(p), n, target))


# LaTeX names its paths in braces, not quotes, so rule 3 above never saw them.
# It cost a real defect: `\\graphicspath{{../../../figures/}}` resolved one
# directory ABOVE the repository, and the paper could not be rebuilt from a
# fresh clone. Braced paths are checked here, and they must resolve inside.
LATEX_PATH = re.compile(r"\\(?:graphicspath|includegraphics|input|include)"
                        r"(?:\[[^\]]*\])?\{+([^{}]+)\}+")


def check_latex_paths(findings):
    """Rule 3, on the surface LaTeX writes: every braced path resolves here.

    An included graphic is looked up through `\\graphicspath`, not relative to
    the file that names it, so the search roots are collected first. Flagging
    a figure that graphicspath resolves perfectly well would teach the reader
    to ignore this rule.
    """
    texs = [q for q in sorted(ROOT.rglob("*.tex"))
            if not any(part in SKIP_DIRS for part in q.parts)]
    roots = []
    for q in texs:
        for m in re.finditer(r"\\graphicspath\{((?:\{[^{}]*\})+)\}",
                             q.read_text(encoding="utf-8", errors="replace")):
            for g in re.findall(r"\{([^{}]*)\}", m.group(1)):
                roots.append(q.parent / g)
    for p in texs:
        for n, line in enumerate(p.read_text(encoding="utf-8",
                                             errors="replace").splitlines(), 1):
            for m in LATEX_PATH.finditer(line):
                target = m.group(1).strip()
                if not target or target.startswith("\\"):
                    continue
                cands = [p.parent / target] + [r / target for r in roots]
                if any(c.exists() or c.with_suffix(".tex").exists()
                       for c in cands):
                    continue
                findings.append(("path", rel(p), n,
                                 f"LaTeX path resolves nowhere in the "
                                 f"repository: {target}"))


# Prose that ADMITS ZERO where the code requires strict non-vanishing. The
# translation rendered a strictly-away-from-zero bound as a bound "by 0",
# which is not the same statement: division, branch choice and biholomorphy
# all need the value to stay AWAY from zero. The code was right; the certificate text was not,
# and it shipped. These spellings are forbidden so the weakening cannot come
# back silently.
WEAKENED_STRICTNESS = re.compile(
    r"bounded (?:below|above) (?:by|at|>=|≥|\\geq?) (?:0|zero|0\.0)\b"
    r"|bounded below > 0", re.I)


def check_strictness(findings):
    """A statement that admits zero where the check is strict.

    This is a REGRESSION GUARD: it knows the spellings that shipped and their
    plain variants (0, zero, 0.0; by, at, >=). It does not read every
    inequality of the text for its meaning, and should not be mistaken for a
    rule that could."""
    for p in files():
        for n, line in enumerate(p.read_text(encoding="utf-8",
                                             errors="replace").splitlines(), 1):
            if rel(p) == "tools/check_distribution.py":
                continue
            if WEAKENED_STRICTNESS.search(line):
                findings.append(("strictness", rel(p), n, line.strip()[:100]))


# Two adjacent prose lines that are IDENTICAL are the fingerprint of a
# translation that collapsed two DISTINCT statements into one copy of the
# second. It happened in the construction rule of the full-cell charts: the
# principal-branch row was deleted and the rotated one duplicated, so the
# docstring said both rows avoid the same half-axis -- not complementary, and
# contradicting the code. The structural fingerprint cannot see this by
# construction: it blanks string contents. A cheap textual rule can.
def check_duplicated_prose(findings):
    """Adjacent identical prose lines: a collapsed pair of statements."""
    for p in files():
        if p.suffix not in {".py", ".md", ".tex"}:
            continue
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
        for i in range(len(lines) - 1):
            a, b = lines[i].strip(), lines[i + 1].strip()
            # Markup is allowed to repeat: a LaTeX table preamble legitimately
            # names identical column specifications one after another. The rule
            # is about PROSE, so lines that begin as markup are skipped.
            if a.startswith(("\\", ">{", "|", "-", "=", "%")):
                continue
            if len(a) > 40 and a == b:
                findings.append(("duplicated-prose", rel(p), i + 1, a[:100]))


def check_imports(findings, allowed):
    """Rule 4 — every import resolves here, to stdlib, or to a pinned dep."""
    local = {q.stem for q in (ROOT / "src").rglob("*.py")} if (
        ROOT / "src").exists() else set()
    for p in sorted(ROOT.rglob("*.py")):
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        try:
            tree = ast.parse(p.read_text(encoding="utf-8", errors="replace"))
        except SyntaxError as e:
            findings.append(("import", rel(p), e.lineno or 0, "syntax error"))
            continue
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.level:      # relative import, inside the package
                    continue
                if node.module:
                    names = [node.module.split(".")[0]]
            for name in names:
                if (name in STDLIB_OK or name in allowed or name in local
                        or name == "k3_atlas"):
                    continue
                findings.append(("import", rel(p), node.lineno, name))


def pinned_dependencies():
    lock = ROOT / "requirements.lock"
    if not lock.exists():
        return set()
    return {l.split("==")[0].strip() for l in
            lock.read_text(encoding="utf-8").splitlines()
            if l.strip() and not l.startswith("#")}


# Appendix E states counts -- fields per certificate, the design records'
# counters, how many of the nine covered files carry `checks`. Every one of
# them was hand-counted once and one of them was wrong in a shipped revision
# (the design records were said to carry no counters; they all do). A count
# that a reader can falsify in one command must be RECOMPUTED by the lint
# from the certificates, and compared with what the paper prints.
APPENDIX_E_NINE = ("open_chart_theorem", "quantitative_atlas", "glue_obligations",
                   "smoothness_and_transitions", "sigma_floor_correction",
                   "atlas_coverage", "bridge_atlas_panel", "bridge_metric_path",
                   "face_traversal_leaf")
APPENDIX_E_DESIGN = ("closure_skeleton", "exact_transitions", "gluing_contract",
                     "uniform_chart_lemma")


def check_appendix_e(findings):
    """Appendix E counts, recomputed from certificates/*.json."""
    certs = {}
    for p in sorted((ROOT / "certificates").glob("*.json")):
        try:
            certs[p.stem] = json.loads(p.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            return
    md = ROOT / "paper" / "certified_k3_atlas.md"
    if not md.exists() or len(certs) != 14:
        findings.append(("appendix-e", "certificates", 0,
                         f"{len(certs)} certificates, appendix E counts fourteen"))
        return
    text = md.read_text(encoding="utf-8")
    start = text.find("## Appendix E")
    end = text.find("## Appendix F", start)
    sec = text[start:end] if start >= 0 else ""
    # the table: one row per field group, count = number of files carrying it
    for m in re.finditer(r"^\| ((?:`[a-z_]+`(?:, )?)+) \| (\d+) \|$", sec, re.M):
        fields = re.findall(r"`([a-z_]+)`", m.group(1))
        claimed = int(m.group(2))
        for f in fields:
            actual = sum(1 for c in certs.values() if f in c)
            if actual != claimed:
                findings.append(("appendix-e", "paper/certified_k3_atlas.md", 0,
                                 f"table says `{f}` = {claimed}, certificates say {actual}"))
    # the design records' counters
    counts = ", ".join(f"{certs[d].get('checks_passed')}/{certs[d].get('checks_total')}"
                       for d in APPENDIX_E_DESIGN)
    if f"({counts})" not in sec:
        findings.append(("appendix-e", "paper/certified_k3_atlas.md", 0,
                         f"design records carry checks {counts}; the text does not say so"))
    # E.2: the red-marker spellings, read from the replayed producers
    markers = set()
    for q in sorted((ROOT / "verification" / "producers").glob("*.py")):
        markers |= set(re.findall(r"_(checks_red|refused)\b", q.read_text(encoding="utf-8")))
    for mk in sorted(markers):
        if f"`*_{mk}`" not in sec:
            findings.append(("appendix-e", "paper/certified_k3_atlas.md", 0,
                             f"a producer collapses a red outcome to `*_{mk}`; E.2 does not name it"))
    # E.5: the field it names must exist, with the value it states
    m = re.search(r"`([a-z_]+)\.([a-z_]+): (true|false)`", sec)
    if not m:
        findings.append(("appendix-e", "paper/certified_k3_atlas.md", 0,
                         "E.5 names no typed field for the two-sided comparison"))
    else:
        val = certs.get("bridge_metric_path", {}).get(m.group(1), {}).get(m.group(2), "absent")
        if str(val).lower() != m.group(3):
            findings.append(("appendix-e", "paper/certified_k3_atlas.md", 0,
                             f"E.5 says {m.group(1)}.{m.group(2)}: {m.group(3)}; the certificate has {val}"))
    # the nine covered files
    with_checks = sum(1 for c in APPENDIX_E_NINE if "checks" in certs.get(c, {}))
    m = re.search(r"command covers, (\d+) carry `checks`", sec)
    if not m or int(m.group(1)) != with_checks:
        findings.append(("appendix-e", "paper/certified_k3_atlas.md", 0,
                         f"{with_checks} of the nine covered files carry `checks`; text says "
                         f"{m.group(1) if m else 'nothing'}"))


# Section 7.4 quotes the nerve of the bridge panel. Every number there is a
# counter of the certificate, and the edge total is a SUM of three of them:
# the 64 upper continuations join no node of the nerve and are counted
# separately, so the four counters of the block sum to more than the edge
# count the paper gives. Both readings are checked, so that neither can drift.
NERVE_COUNTERS = ("n_nodes", "n_lower_nodes", "n_bridge_nodes",
                  "n_lower_lower_edges_imported", "n_bridge_lower_edges_certified",
                  "n_bridge_bridge_edges_certified", "n_new_triples_involving_a_bridge")


def _spellings(n: int) -> set[str]:
    out = {str(n)}
    if n >= 1000:
        out.add(f"{n // 1000}{{,}}{n % 1000:03d}")
    return out


def check_nerve_counts(findings):
    """Section 7.4 against the `nerve` block of the bridge panel."""
    p = ROOT / "certificates" / "bridge_atlas_panel.json"
    md = ROOT / "paper" / "certified_k3_atlas.md"
    if not p.exists() or not md.exists():
        return
    try:
        nerve = json.loads(p.read_text(encoding="utf-8"))["nerve"]
    except (ValueError, OSError, KeyError):
        findings.append(("nerve-counts", rel(p), 0, "no nerve block"))
        return
    text = md.read_text(encoding="utf-8")
    start = text.find("### 7.4")
    end = text.find("## 8.", start)
    sec = text[start:end] if start >= 0 else ""
    for k in NERVE_COUNTERS:
        if not any(s in sec for s in _spellings(int(nerve[k]))):
            findings.append(("nerve-counts", "paper/certified_k3_atlas.md", 0,
                             f"section 7.4 does not quote {k} = {nerve[k]}"))
    edges = (int(nerve["n_lower_lower_edges_imported"])
             + int(nerve["n_bridge_lower_edges_certified"])
             + int(nerve["n_bridge_bridge_edges_certified"]))
    four = edges + int(nerve["n_bridge_upper_cont_certified"])
    for n, what in ((edges, "certified edges of the nerve"),
                    (four, "sum of the four counters")):
        if not any(s in sec for s in _spellings(n)):
            findings.append(("nerve-counts", "paper/certified_k3_atlas.md", 0,
                             f"section 7.4 does not quote {n} ({what})"))


# Three shipped inputs record only a 16-character prefix as their own source
# fingerprint; they predate the full-hash rule. Their bytes are pinned in
# full by the certificates that read them. The README names them; this check
# keeps the name list equal to the measured one, in both directions.
SHORT_SELF_SHA = ("chart_selection_criterion.json", "dyadic_cover.json",
                  "residual_closure.json")


def check_self_sha(findings):
    """Every `self_sha256` is 64 characters, except the three declared ones."""
    readme = (ROOT / "README.md").read_text(encoding="utf-8") if (
        ROOT / "README.md").exists() else ""
    short = []
    for p in sorted(list((ROOT / "certificates").glob("*.json"))
                    + list((ROOT / FROZEN_INPUTS).glob("*.json"))):
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            continue
        s = d.get("self_sha256") or (d.get("provenance") or {}).get("self_sha256")
        if isinstance(s, str) and len(s) != 64:
            short.append(p.name)
            if len(s) != 16:
                findings.append(("self-sha", rel(p), 0, f"self_sha256 of length {len(s)}"))
    if tuple(short) != SHORT_SELF_SHA:
        findings.append(("self-sha", "tools/check_distribution.py", 0,
                         f"declared {list(SHORT_SELF_SHA)}, measured {short}"))
    for name in SHORT_SELF_SHA:
        if f"`{name}`" not in readme:
            findings.append(("self-sha", "README.md", 0,
                             f"{name} records a 16-character prefix and the README does not say so"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", action="store_true",
                    help="counts per rule, no per-finding lines")
    ap.add_argument("--rule", choices=("english", "vocabulary", "json-key",
                                       "json-text", "json-vocabulary", "path",
                                       "strictness", "duplicated-prose", "import",
                                       "appendix-e", "nerve-counts", "self-sha"),
                    help="report a single rule")
    args = ap.parse_args()

    findings: list[tuple[str, str, int, str]] = []
    check_english(findings)
    check_vocabulary(findings)
    check_json_keys(findings)
    check_paths(findings)
    check_latex_paths(findings)
    check_strictness(findings)
    check_duplicated_prose(findings)
    check_imports(findings, pinned_dependencies())
    check_appendix_e(findings)
    check_json_vocabulary(findings)
    check_nerve_counts(findings)
    check_self_sha(findings)

    if args.rule:
        findings = [f for f in findings if f[0] == args.rule]

    by_rule: dict[str, int] = {}
    by_file: dict[str, int] = {}
    for rule, path, _n, _what in findings:
        by_rule[rule] = by_rule.get(rule, 0) + 1
        by_file[path] = by_file.get(path, 0) + 1

    print("Distribution check")
    report_exemptions()
    for rule in ("english", "vocabulary", "json-key", "json-text", "path",
                 "strictness", "duplicated-prose", "import", "appendix-e"):
        if args.rule and rule != args.rule:
            continue
        print(f"  {rule:12s} {by_rule.get(rule, 0):6d}")
    if not args.summary:
        for rule, path, n, what in findings[:200]:
            print(f"  {rule:11s} {path}:{n}  {what}")
        if len(findings) > 200:
            print(f"  … and {len(findings) - 200} more")
    print()
    print(f"  worst files: " + ", ".join(
        f"{f} ({c})" for f, c in sorted(by_file.items(), key=lambda x: -x[1])[:5]))
    print()
    print(f"VERDICT: {'PASS' if not findings else f'FAIL — {len(findings)} finding(s)'}")
    return 0 if not findings else 1


if __name__ == "__main__":
    sys.exit(main())
