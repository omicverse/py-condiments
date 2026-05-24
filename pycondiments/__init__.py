"""pycondiments — Python port of condiments (Roux de Bezieux et al. Nat Commun 2024).

Differential-trajectory analysis between conditions. Given a Slingshot trajectory
+ per-cell condition labels:

- `imbalance_score`: per-cell condition-imbalance metric via multinomial test on kNN
- `topologyTest`: does trajectory topology differ between conditions?
- `progressionTest`: do cells advance at different rates per condition?
- `differentiationTest`: do cells make different lineage choices?
- `weights_from_pst`: convert pseudotime to lineage weights
- `merge_sds`: utility to merge SlingshotDataSet objects

v0.1 covers all 7 exports. Plotting deferred.

NOTE: `topologyTest` / `progressionTest` / `differentiationTest` ultimately rely
on KS-tests and Wald-tests; pure-Python via scipy.stats. The `fateSelectionTest`
multi-lineage variant depends on py-tradeSeq which is also at v0.1; both
v0.1's ship in lockstep.
"""

from __future__ import annotations

__version__ = "0.1.0"

from .imbalance import imbalance_score
from .progression import progressionTest, differentiationTest, weights_from_pst
from .topology import topologyTest, create_differential_topology
from .utils import merge_sds

__all__ = [
    "imbalance_score",
    "progressionTest",
    "differentiationTest",
    "weights_from_pst",
    "topologyTest",
    "create_differential_topology",
    "merge_sds",
    "__version__",
]
