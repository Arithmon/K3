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
python3 verification/verify.py            # full verification, about two minutes
python3 verification/verify.py --quick    # hashes only, no replay
```

The command needs numpy, sympy and mpmath **1.3.0** (the pin is enforced:
the certificates use directed rounding and record the arithmetic backend
they were produced with). It has been run under Python 3.12 and 3.14. The
reference environment for reproducing the expensive certificates is pinned
in `requirements.lock`, `pyproject.toml` and the `Dockerfile`, and described
in `docs/ENVIRONMENT.md`.

## What the command checks

**Replayed** (five certificates, under a second each). The producer in
`verification/producers/` is re-executed and the result is checked: every
check green, every negative control green, the recorded outcome equal to
the expected one. A regenerated file is deliberately not compared by hash:
its provenance records the commit it was built from.

**Recomputed** (the coverage enumeration: 71,807,792 boxes, about seventy
seconds on four cores). Ten fields are compared one by one with the shipped
file, gauge by gauge. This is the computation carrying the uniform pivot
floor; reading its verdict back would check nothing.

**Hashed** (the bridge panel, the metric path, the face traversal). These
take hours; their SHA-256 is checked against the value recorded here and
published in Appendix F of the paper, in the markdown, the LaTeX and the
PDF. Their producers ship in `src/k3_atlas` with every input they read
(Appendix F.4 of the paper).

The command also checks every `upstream` pin a certificate wrote against the
shipped bytes (five pins of the design records are stale, printed with what
was measured about the difference), and `docs/RESULTS_INDEX.md` against what
it actually does. It is read-only on the tree.

**Negative controls.** Each producer carries, besides its checks, deliberate
perturbations of its own computation, and the check that should catch each
one is required to fail. Both counts appear in the output.

## Layout

```
paper/               the paper: markdown source, LaTeX, PDF
certificates/        the 14 certificates, as JSON
figures/             the paper's figures
docs/
  RESULTS_INDEX.md   every result of the paper: generator, artifact, how verified
  ENVIRONMENT.md     the pinned environment and what each certificate recorded
verification/
  verify.py          the one-command entry point
  producers/         the programs behind the replayed and recomputed certificates
    model.py         the surface itself, defined once and imported
src/k3_atlas/        the producers behind the hash-verified certificates and
                     their shared numerical kernels, as an installable package
  data/              the 25 inputs that chain reads (41 MB)
Dockerfile, pyproject.toml, requirements.lock    the pinned environment
```

## Provenance

Each certificate carries a provenance block: the SHA-256 of its own source,
the versions of the arithmetic backends, timings, and the revision of the
repository it was built from. That revision points at the private working
repository where the results were produced and will not resolve publicly; it
is kept because it is what the producers recorded. Three inputs of the
expensive chain (`chart_selection_criterion.json`, `dyadic_cover.json`,
`residual_closure.json`) record a 16-character prefix as their own source
fingerprint; their bytes are pinned in full by the certificates that read
them.

One field was removed: the open-chart certificate recorded the fingerprint
of an internal review note, purely to timestamp it. Both the note and the
field are absent here. This is the only content difference between these
producers and the ones that were run.

## Citing

See `paper/latex/certified_k3_atlas.pdf` for the full statement of results
and the bibliography. The archival deposit is the concept DOI, which always
resolves to the latest version:

> de La Fournière, B. *Certified Analytic Geometry on an Explicit K3 Surface.*
> Zenodo. https://doi.org/10.5281/zenodo.22047469
