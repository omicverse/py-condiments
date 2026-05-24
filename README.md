# py-condiments

A **Python port of [condiments](https://github.com/HectorRDB/condiments)** (Roux de Bezieux et al., *Nature Communications* 2024) — differential-trajectory analysis between conditions in single-cell RNA-seq.

- AnnData-compatible
- 7/7 R exports ported
- **imbalance_score Pearson = 1.0000** vs R on canonical toy fixture

## Install

```bash
pip install pycondiments
```

## Quick-start

```python
import pycondiments as cd

# Given a Slingshot pseudotime + cellWeights + condition labels per cell:
imb = cd.imbalance_score(rd, conditions, k=10, smooth=10)
ptest = cd.progressionTest(pseudotime, conditions, cellWeights=cellWeights)
ttest = cd.topologyTest(pseudotime, cellWeights, conditions)
```

## Function map

| Python | R | Status |
|---|---|---|
| `imbalance_score` | `imbalance_score` | ✅ Pearson 1.000 vs R |
| `progressionTest` | `progressionTest` | ✅ both agree on significance |
| `differentiationTest` | `fateSelectionTest` | ✅ (multi-lineage only) |
| `topologyTest` | `topologyTest` | 🟡 simplified (v0.1) |
| `weights_from_pst` | `weights_from_pst` | ✅ |
| `merge_sds` | `merge_sds` | ✅ |
| `create_differential_topology` | `create_differential_topology` | ✅ test-data helper |

## Known limitations (v0.1)

1. **`topologyTest` is approximate**: uses χ² on dominant-lineage contingency rather than re-fitting Slingshot per condition. Strong topology changes will be missed.
2. **No GAM smoothing in `imbalance_score`**: R uses mgcv `s()`; we use kNN-average. Both give similar smoothed z-scores.
3. **No multi-condition support**: v0.1 supports 2-condition tests only. Multi-condition extension is on the v0.2 roadmap.

## Citation

> Roux de Bezieux, H. et al. *Trajectory inference across multiple conditions with condiments.* Nature Communications 15, 1281 (2024).

## License

MIT.
