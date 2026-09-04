# Certified K3 atlas

Machine-checkable certificates for the results of the accompanying paper,
`paper/latex/certified_k3_atlas.pdf`.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22047469.svg)](https://doi.org/10.5281/zenodo.22047469)

The surface is the complete intersection of three diagonal quadrics in
`P^5`, with the six parameters `(1, 2, 3, 5, 7, 11)`. The paper builds an
explicit atlas on it, with certified chart radii, certified transition maps,
and a coverage argument; this repository lets a reader re-run the arguments
that are cheap to re-run and check the fingerprints of the ones that are not.

## Verifying

```
python3 verification/verify.py            # full verification
python3 verification/verify.py --quick    # hashes only, no replay
```

Reference environment (`requirements.lock`, `pyproject.toml`, `Dockerfile`;
details in `docs/ENVIRONMENT.md`):

| package | version |
|---|---|
| python | 3.14.4 |
| numpy | 2.5.1 |
| scipy | 1.18.0 |
| sympy | 1.14.0 |
| matplotlib | 3.11.1 |
| mpmath | **1.3.0 — pinned, enforced** |

The environment is pinned, not recommended: `pip install -e .` requires
Python 3.14.4 exactly. The verification command itself imports only numpy,
sympy and mpmath (matplotlib serves the figure script), and has also been run outside the pinned
environment, under Python 3.12 with mpmath 1.3.0.

The mpmath pin is checked before anything runs, and is not a suggestion.
The certificates use directed rounding, and they serialise the arithmetic
backend into their own provenance; a different mpmath would move the last
digits of every interval endpoint, so a "verification" against another
version would be comparing two different computations.

## The three levels, and why they differ

**Replayed** (five certificates, under a second each). The producer in
`verification/producers/` is re-executed and the resulting certificate is
checked: every check green, every negative control green, and the recorded
outcome equal, exactly, to the expected one.

A replayed certificate is deliberately **not** compared by hash. Each one
records the source revision it was built from, so its bytes change with
every commit to the repository that produced it; hashing a regenerated file
would test version-control history rather than mathematics. What is checked
is the content: the checks and the outcome.

**Recomputed** (one certificate: the coverage enumeration). Its
branch-and-bound is re-run — 71,807,792 boxes, about seventy seconds on four
cores — and ten fields are compared one by one with the shipped file, gauge by
gauge. This is the computation carrying the uniform pivot floor, so reading
its verdict back would check nothing; the recomputation reddens if a single
counter moves. The verifier also checks that the recomputation actually
rewrote the file it compares, so that a producer writing elsewhere cannot make
the comparison vacuous.

**Hashed** (three certificates: the bridge panel, the metric path, the face
traversal). These take hours on a compute machine and are not re-run here;
their SHA-256 is checked against the recorded value. Reproducing them is
possible — the producers are ordinary Python — but it is not something a
one-command check should ask of a reader.

Verification is read-only on the tree: replaying rewrites the certificates,
so their original bytes are saved and restored.

## Negative controls

Each producer carries, besides its checks, a set of *negative controls*: the
computation is deliberately perturbed — a coefficient moved, a sign flipped,
a bound loosened — and the check that should catch the perturbation is
required to fail. A check that passes on correct input tells you little; a
check that also fails on wrong input is the one worth reading. Both counts
appear in the verification output.

## Layout

```
paper/           the paper, as markdown, LaTeX source and PDF
certificates/    the 14 certificates, as JSON
figures/         the paper's figures
verification/
  verify.py      the one-command entry point
  producers/     the programs that emit the replayed and recomputed certificates
    model.py     the surface itself, defined once and imported
src/k3_atlas/    the producers of the three hash-verified certificates and
                 their shared numerical kernels, as an installable package
  data/          the 25 input artifacts that chain reads (41 MB)
docs/
  RESULTS_INDEX.md   every result of the paper: generator, artifact, mode, cost
  ENVIRONMENT.md     the pinned environment and what each certificate recorded
tools/           the distribution check (language, vocabulary, paths, counts)
Dockerfile       the pinned environment, materialised
pyproject.toml   the package, with the same pins
requirements.lock
```

`model.py` holds the defining data of the surface. Every certificate
derives the quadrics and the Vandermonde blocks from it rather than
restating them, so that changing the surface would make the checks fail
rather than let them quietly disagree.

## Provenance, and what was removed

Each certificate carries a provenance block: the SHA-256 of its own source,
the versions of the arithmetic backends, timings, and the revision of the
repository it was built from. That revision field points at the private
working repository where these results were produced, so the identifiers
will not resolve publicly. It is kept because it is what the producers
actually recorded, and rewriting a provenance field to look tidier would
make it worthless.

The same rule governs three things a reader will notice. Three shipped inputs
of the expensive chain (`chart_selection_criterion.json`, `dyadic_cover.json`,
`residual_closure.json`) record only a 16-character prefix as their own
source fingerprint: they predate the full-hash rule, and their bytes are
pinned in full by the `upstream` blocks of the certificates that read them.
The inputs and the lineage fields keep their development vocabulary as
recorded: the run that produced each row of `exact_gluing.json`, the
workspace paths of the retraction registry, the names of the files a design
record read. And three prose leaves of the bridge panel (`not_paid_here`) are
half French: the panel is hash-verified and pinned by two other hash-verified
certificates, so a hand edit would break those pins; the producer already
carries the corrected text, and the next full replay of the chain will ship
it. The distribution check (`tools/check_distribution.py`) counts what these
allowances let through and prints the counts on every run; anywhere else the
same words are a failure.

Four design records pin inputs that were regenerated or translated after the
records were written, and their producers are not shipped. `verify.py`
checks every `upstream` pin it can map against the shipped bytes, prints the
stale ones with what was measured about the difference, and fails on any
other mismatch.

One field was removed. The open-chart certificate recorded the fingerprint
of an internal review note, purely to timestamp it; the note plays no
mathematical role in the certificate. Publishing an internal review note
would expose the workshop rather than the proof, so both the note and the
field that hashed it are absent here. This is the only content difference
between these producers and the ones that were run.

## Citing

See `paper/latex/certified_k3_atlas.pdf` for the full statement of results and the
bibliography. The archival deposit is the concept DOI, which always resolves to the latest
version:

> de La Fournière, B. *Certified Analytic Geometry on an Explicit K3 Surface.*
> Zenodo. https://doi.org/10.5281/zenodo.22047469

Cite the concept DOI rather than a version DOI unless you need to refer to one
specific deposited version.
