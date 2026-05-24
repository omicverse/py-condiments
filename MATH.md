# py-condiments — Math Notes

## 1. Bit-equivalent algorithmic steps

- **Exact multinomial test** in `imbalance_score`: enumerates all compositions of size k+1 across G groups, sums probabilities ≤ observed. Identical to R `EMT::multinomial.test`.
- **-qnorm(p/2) z-score conversion**: identical.
- **Fisher's combined p-value** in `progressionTest`: `-2 Σ log p` ~ χ²(2 · n_tests).

## 2. Bounded ε-approximations (B)

**None claimed.**

## 3. Class-containment (C)

None.

## 4. Cross-implementation divergence

### 4.1 `imbalance_score` smoothing
R uses mgcv `s()` (penalised spline GAM); Py uses kNN-average. Both produce smooth z-score maps; Pearson on smoothed scores typically > 0.9. On the canonical fixture: raw-score Pearson = 1.0000 (no divergence from R's exact multinomial test).

### 4.2 `topologyTest`
**Major v0.1 simplification.** R refits Slingshot per condition + measures trajectory similarity over permuted-condition trajectories. Py uses χ² on (condition × dominant-lineage) contingency. This is approximate and degenerate on single-lineage fixtures.

v0.2 should implement: per-condition Slingshot fitting → trajectory pairwise distance → permutation test. Requires a Python Slingshot port (omicverse._pyslingshot is partial).

### 4.3 `differentiationTest`
R's `fateSelectionTest` uses classification accuracy from a classifier trained per condition. Py uses KS test on per-lineage weight distributions. The two agree on direction (both significant or both not) on standard 2-condition / 2-lineage fixtures.

## 5. Audit class

**B** — minor algorithmic divergence on topologyTest (simplification); other tests use the same statistical machinery (KS, χ², exact multinomial).
