# Contributing to K3

This repository exists so that the quantitative claims of the paper can be
checked rather than believed. The most valuable contribution is therefore the
one that makes a certificate fail: a gate that passes when it should not, a
bound that does not hold, a replay that does not reproduce. Adversarial input
is welcome by design.

Org-wide rules (what helps, what does not, house style) are in the
[organization CONTRIBUTING](https://github.com/arithmon/.github/blob/main/CONTRIBUTING.md).
This file covers the four procedures this repository owns.

## Report a certificate that does not replay

Five certificates are re-executed by `python3 verification/verify.py`. If one
of them fails on your machine, that is worth an issue even when the cause
turns out to be local.

1. Open an issue titled `replay: <certificate name>`.
2. Give the full output of `python3 verification/verify.py`, your Python
   version, and the versions of numpy, sympy and mpmath.
3. Check the mpmath pin first. The certificates record the arithmetic backend
   they were produced with, and version 1.3.0 is enforced before anything
   runs, because directed rounding is not comparable across versions.
4. A replay that fails for a reason other than the environment is a defect in
   the certificate, and it is fixed in the open with the old value left
   visible.

## Challenge a hashed certificate

Three certificates (the bridge panel, the metric path, the face traversal)
are checked by SHA-256 rather than replayed, because reproducing them takes
hours on a compute machine. That is a convenience, not a privilege.

1. Their producers are ordinary Python and can be re-run in full.
2. If your run disagrees with the recorded hash, report the disagreement with
   your output, whatever it says.
3. A hashed certificate that cannot be reproduced independently is treated as
   unsupported until it is, and the paper is corrected rather than defended.

## Report a defect in the mathematics

The paper states its own limits: no Ricci-flat or hyperkähler metric is
claimed, the witnesses are found on deterministic rational grids rather than
by optimisation, and several constants are deliberately loose. A defect
report is most useful when it separates these three cases.

1. Say which the claim is: an error in a proof, a constant that is not sharp,
   or a hypothesis that is asserted rather than certified.
2. Point to the specific certificate and gate, or to the section of the paper.
3. An error in a proof is a correction. A loose constant is an improvement. A
   hypothesis found to be asserted rather than certified is the most serious
   of the three, and is treated as such.

One such defect was found while the paper was being written, and it is
recorded rather than quietly repaired: an upstream argument bounded a product
where the minimum was required, which moved the uniform certified radius from
9.60e-10 to 2.11e-12. The statement survived, the constant did not. See
`certificates/sigma_floor_correction.json`, which carries an explicit witness
on the surface and a negative control showing that the witness detects the
substitution itself and not merely a floor set too high.

## Improve the negative controls

Each producer carries gates and, separately, negative controls: the
computation is deliberately perturbed and the gate that should catch the
perturbation is required to fail. A gate that passes on correct input tells
you little. A gate that also fails on wrong input is the one worth reading.

1. A proposed negative control should name the perturbation and the gate it
   targets.
2. It is accepted when it fails on the perturbed input and passes on the
   real one, and when it is not already covered by an existing control.
3. A control that no existing gate can catch is the most useful of all: it
   means a gate is missing, not that the control is wrong.

---

Siblings: [K7](https://github.com/arithmon/K7) ·
[Program](https://github.com/arithmon/program) ·
[Sieve](https://github.com/arithmon/sieve) ·
[Atlas](https://github.com/arithmon/atlas) ·
[Lean](https://github.com/arithmon/lean)

<sub>K₇ (formerly GIFT) is the founding framework of the Arithmon program.
Program: [arithmon.com](https://arithmon.com) ·
[github.com/arithmon](https://github.com/arithmon)</sub>
