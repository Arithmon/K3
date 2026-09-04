# Results index

Every computational claim of the paper, the program that produces it, the artifact it produces, and the reproduction mode that rebuilds it.

This file is generated, not written: it is derived from an inventory of the producing workspace (each certificate's producer and the transitive closure of what that producer reads) and regenerated at every release. Do not edit it by hand: the next generation overwrites it.

| Paper result | Generator | Reference artifact | Mode | Current check | Closure (programs) | Closure (lines) | Inputs to ship |
| --- | --- | --- | --- | --- | ---: | ---: | ---: |
| Sec 7.4 (Exact overlap closure and the certified nerve) -- certified bridge domains, nerve and overlap closure | `bridge_continuation.py` | `bridge_atlas_panel.json` | full | hashed | 19 | 13815 | 9 |
| Sec 9 (Validated Metric Compatibility on Selected Continuation Domains) -- validated metric compatibility on continuation domains | `metric_transport.py` | `bridge_metric_path.json` | full | hashed | 20 | 14456 | 9 |
| Sec 4.3 (The atlas theorem, and the glued space as a corollary) -- design record behind the closure theorem | `closure_skeleton.py` | `closure_skeleton.json` | full | unverified | 2 | 1782 | 6 |
| Sec 5.2 (Exact cocycle law) -- design record behind the exact cocycle law | `exact_transitions.py` | `exact_transitions.json` | full | unverified | 2 | 1751 | 2 |
| Sec 7.4 (Exact overlap closure and the certified nerve) -- face crossing, control experiment | `face_continuation.py` | `face_traversal_leaf.json` | full | hashed | 20 | 14506 | 8 |
| Figures (Figures) -- manifest of every plotted number | `figures.py` | `figures_manifest.json` | full | unverified | 1 | 315 | 4 |
| Sec 4.2 (Uniformity and finiteness) -- design record behind the regional gluing contract | `gluing_contract.py` | `gluing_contract.json` | full | unverified | 1 | 299 | 16 |
| Sec 3.4 (A quantitative inverse-function lemma) -- design record behind the chart lemma | `uniform_chart_lemma.py` | `uniform_chart_lemma.json` | full | unverified | 2 | 1783 | 4 |
| Sec 4.1 (Global admissibility) -- exhaustive box enumeration behind the pivot floor | `atlas_coverage.py` | `atlas_coverage.json` | standard | recomputed | 1 | 367 | 0 |
| Sec 5.3 (Finite transition generators) -- transition generators, words of length at most 4 | `glue_obligations.py` | `glue_obligations.json` | standard | replayed | 1 | 343 | 3 |
| Sec 3.6 (The certified radius) -- certified radius, uniform open-chart theorem | `open_chart_theorem.py` | `open_chart_theorem.json` | standard | replayed | 2 | 1807 | 2 |
| Sec 4.3 (The atlas theorem, and the glued space as a corollary) -- atlas theorem, X_atlas isomorphic to X | `quantitative_atlas.py` | `quantitative_atlas.json` | standard | replayed | 1 | 131 | 5 |
| Appendix B (Appendix B -- Proof of Proposition 3.2) -- sigma floor defect, 9.6e-10 to 2.1e-12 | `sigma_floor_correction.py` | `sigma_floor_correction.json` | standard | replayed | 1 | 242 | 1 |
| Sec 2.2 (Smoothness and the K3 property) -- three-nonzero lemma, smoothness and transitions | `smoothness_and_transitions.py` | `smoothness_and_transitions.json` | standard | replayed | 1 | 257 | 2 |

## Reproduction modes

`compare` checks integrity and schema only, and claims no reproduction. `standard` rebuilds every result that is reasonably short. `full` additionally rebuilds the expensive chain; it takes hours, which is the point of separating it from the daily command.

## Measured cost of `full`

- programs in the transitive closure: **32** (18224 lines)
- shared numerical kernels: `kahler_metric`, `invariant_quotient_ring`, `interval_arithmetic`, `spectral_basis`, `taylor_models`
- artifacts that must be present before the run: **42** (45.63 MB, binary data included)
- of those, rebuildable by also running their own generator: **42**
- artifacts the run rewrites: **0** (the certificates themselves)

A generator imported for its functions does not rebuild its own artifact: importing is not running. Every input above must either ship, or be rebuilt by an explicit earlier step of `full`.
- largest single input: 29.59 MB

## The four design certificates

Four of the shipped certificates are design records rather than claims of the theorems. Making them reproducible pulls **23** further artifacts into the repository (transitive closure of the four, minus the 14 shipped), or **22** counting only what the verified chain does not already read: contract amendments, preregistrations, and an internal chain none of the theorems rest on.

Generated for release 2.0.0 of this repository.
