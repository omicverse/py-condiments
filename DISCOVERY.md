# Discovery — py-condiments (scaffold; v0.0 not started)

## 1. Is this package already ported?

`python -m engine.discover_omicverse_deps --check condiments` → **No existing omicverse port found.**

## 2. Dependency audit + scope assessment

condiments (Roux de Bezieux et al. *Nat. Commun.* 2024, 38 citations) provides differential-trajectory analysis between conditions on a Slingshot+tradeSeq output:
- `topologyTest` — does trajectory topology differ between conditions?
- `progressionTest` — do cells advance at different rates?
- `differentiationTest` — do cells make different lineage choices?
- `imbalance_score` — per-cell imbalance metric
- `weights_from_pst` — convert pseudotime to lineage weights
- `merge_sds` — utility
- `create_differential_topology` — synthesize a test fixture

**Size**: 7 exports, **2123 R LOC**. Most tractable of the remaining 3.

**Critical dependency**: condiments builds ON top of `tradeSeq::fitGAM` (uses GAM Wald statistics for differential-progression tests). So **py-condiments depends on py-tradeSeq being functional** — and py-tradeSeq v0.1-alpha is WIP (`associationTest` passes; `startVsEndTest` 0.45 Spearman, needs joint-fit fix).

**Reusable from omicverse**:
- Eventually `py-tradeSeq` once it ships v0.1
- `omicverse.single._pyslingshot` for the upstream Slingshot output
- `scipy.stats` for the Kolmogorov–Smirnov tests in `progressionTest`

## 3. Decision

**Wait on py-tradeSeq v0.1**. py-condiments is mostly downstream statistical-test wrappers (KS tests, classification tests, permutation tests) — the algorithmic body is small (~500-700 LOC of "real" code excluding S4 boilerplate). Once py-tradeSeq's `associationTest` + `startVsEndTest` are validated, py-condiments can be ported in ~3-5 days.

## 4. v0.1 roadmap (for a future session)

1. Wait for py-tradeSeq v0.1 to ship.
2. Port `weights_from_pst`, `merge_sds`, `imbalance_score` (utilities, ~200 LOC, no GAM).
3. Port `topologyTest` — KS-style test, fixture-light.
4. Port `progressionTest` — KS / classification on pseudotime per condition.
5. Port `differentiationTest` — classification on lineage weights per condition.
6. Three notebooks; release as `pycondiments 0.1.0`.

Expected effort: ~3-5 days after py-tradeSeq dependency is satisfied.
