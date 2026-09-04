# Results index

Every computational claim of the paper, the program that produces it, the artifact it produces, and how `verification/verify.py` checks it. Generated at release from an inventory of the producing workspace; not edited by hand.

| Paper result | Generator | Artifact | Verified by |
| --- | --- | --- | --- |
| Sec 7.4 (Exact overlap closure and the certified nerve) -- certified bridge domains, nerve and overlap closure | `bridge_continuation.py` | `bridge_atlas_panel.json` | hashed |
| Sec 9 (Validated Metric Compatibility on Selected Continuation Domains) -- validated metric compatibility on continuation domains | `metric_transport.py` | `bridge_metric_path.json` | hashed |
| Sec 4.3 (The atlas theorem, and the glued space as a corollary) -- design record behind the closure theorem | `closure_skeleton.py` | `closure_skeleton.json` | unverified |
| Sec 5.2 (Exact cocycle law) -- design record behind the exact cocycle law | `exact_transitions.py` | `exact_transitions.json` | unverified |
| Sec 7.4 (Exact overlap closure and the certified nerve) -- face crossing, control experiment | `face_continuation.py` | `face_traversal_leaf.json` | hashed |
| Figures (Figures) -- manifest of every plotted number | `figures.py` | `figures_manifest.json` | unverified |
| Sec 4.2 (Uniformity and finiteness) -- design record behind the regional gluing contract | `gluing_contract.py` | `gluing_contract.json` | unverified |
| Sec 3.4 (A quantitative inverse-function lemma) -- design record behind the chart lemma | `uniform_chart_lemma.py` | `uniform_chart_lemma.json` | unverified |
| Sec 4.1 (Global admissibility) -- exhaustive box enumeration behind the pivot floor | `atlas_coverage.py` | `atlas_coverage.json` | recomputed |
| Sec 5.3 (Finite transition generators) -- transition generators, words of length at most 4 | `glue_obligations.py` | `glue_obligations.json` | replayed |
| Sec 3.6 (The certified radius) -- certified radius, uniform open-chart theorem | `open_chart_theorem.py` | `open_chart_theorem.json` | replayed |
| Sec 4.3 (The atlas theorem, and the glued space as a corollary) -- atlas theorem, X_atlas isomorphic to X | `quantitative_atlas.py` | `quantitative_atlas.json` | replayed |
| Appendix B (Appendix B -- Proof of Proposition 3.2) -- sigma floor defect, 9.6e-10 to 2.1e-12 | `sigma_floor_correction.py` | `sigma_floor_correction.json` | replayed |
| Sec 2.2 (Smoothness and the K3 property) -- three-nonzero lemma, smoothness and transitions | `smoothness_and_transitions.py` | `smoothness_and_transitions.json` | replayed |

`replayed`: the producer in `verification/producers/` is re-executed and the result checked. `recomputed`: re-run and compared field by field with the shipped file. `hashed`: the SHA-256 of the shipped file is checked; the producer ships in `src/k3_atlas` with its inputs and takes hours. `unverified`: a design record, shipped for reading.
