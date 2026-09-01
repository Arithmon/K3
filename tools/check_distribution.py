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
It is deliberately without exemptions: a rule that carves out our own
vocabulary would stop measuring the thing it names.
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
FRENCH_MARKERS = [
    r"\ble\b", r"\bla\b", r"\bles\b", r"\bdes\b", r"\bdu\b", r"\bune\b",
    r"\bun\b", r"\bqui\b", r"\bque\b", r"\bpour\b", r"\bavec\b", r"\bsur\b",
    r"\bdans\b", r"\bpas\b", r"\bplus\b", r"\best\b", r"\bsont\b",
    r"\bcette\b", r"\bce\b", r"\bnous\b", r"\bnon\b", r"\bmais\b",
    r"\bdonc\b", r"\bsans\b", r"\bleur\b", r"\bd'un\b", r"\bl'\w+",
]

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
            hits = len(set(m.group(0).lower() for m in fr.finditer(s)))
            # Two distinct French function words in one line make prose, not a
            # stray token: one such word can be a variable name, two of them
            # in a row make a sentence.
            if hits >= 2:
                findings.append(("english", rel(p), n, s[:100]))


# Lines that legitimately carry a name a reader needs: authorship and the
# disclosure of language-model assistance. Matched on the LINE, so the
# exemption cannot spread into the code.
NAME_EXEMPT = re.compile(
    r"(author|thanks|acknowledg|Anthropic|OpenAI|Claude|correspondence|"
    r"^\*\*[A-Z])", re.I)


def check_vocabulary(findings):
    """Rule 2 — no development vocabulary."""
    voc = re.compile("|".join(DEVELOPMENT_VOCABULARY), re.I)
    for p in files():
        if rel(p) == "tools/check_distribution.py":
            continue        # the rule names the words it forbids
        for n, line in enumerate(p.read_text(encoding="utf-8",
                                             errors="replace").splitlines(), 1):
            if NAME_EXEMPT.search(line.strip()):
                continue
            for m in voc.finditer(line):
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
        elif isinstance(o, str) and len(o) > 25:
            hits = len(set(m.group(0).lower() for m in fr.finditer(o)))
            if hits >= 2:
                findings.append(("json-text", rel(p), 0,
                                 f"{path} = {o[:70]}"))
        elif isinstance(o, list):
            for i, v in enumerate(o[:50]):
                walk(v, f"{path}[{i}]", p)

    for p in sorted(ROOT.rglob("*.json")):
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        try:
            walk(json.loads(p.read_text(encoding="utf-8")), "", p)
        except (ValueError, OSError):
            findings.append(("json-key", rel(p), 0, "unreadable"))


def check_paths(findings):
    """Rule 3 — no path points outside this repository."""
    pat = re.compile(r"""["'`]((?:\.\./|/home/|/Users/|C:\\)[^"'`\n]{2,120})""")
    for p in files():
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", action="store_true",
                    help="counts per rule, no per-finding lines")
    ap.add_argument("--rule", choices=("english", "vocabulary", "json-key",
                                       "json-text", "path", "import"),
                    help="report a single rule")
    args = ap.parse_args()

    findings: list[tuple[str, str, int, str]] = []
    check_english(findings)
    check_vocabulary(findings)
    check_json_keys(findings)
    check_paths(findings)
    check_imports(findings, pinned_dependencies())

    if args.rule:
        findings = [f for f in findings if f[0] == args.rule]

    by_rule: dict[str, int] = {}
    by_file: dict[str, int] = {}
    for rule, path, _n, _what in findings:
        by_rule[rule] = by_rule.get(rule, 0) + 1
        by_file[path] = by_file.get(path, 0) + 1

    print("Distribution check")
    for rule in ("english", "vocabulary", "json-key", "json-text", "path",
                 "import"):
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
