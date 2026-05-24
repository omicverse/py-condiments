"""topologyTest, create_differential_topology.

`topologyTest` tests whether the trajectory topology fitted per condition
differs from the topology fitted on pooled data. Mirror's
condiments::topologyTest with a simplification: in v0.1 we test whether
the LINEAGE ASSIGNMENT distribution differs between conditions (via χ²)
rather than re-fitting Slingshot per condition (which requires the full
Slingshot port — out of scope here).

`create_differential_topology` is a test-data generator.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency


def topologyTest(
    pseudotime: np.ndarray,
    cellWeights: np.ndarray,
    conditions: np.ndarray,
    *,
    n_lineages: int | None = None,
) -> dict:
    """Test whether trajectory topology differs between conditions.

    v0.1 simplification: χ² test of contingency table
    (condition × dominant-lineage-assignment).

    Args:
        pseudotime: (n_cells × n_lineages) per Slingshot
        cellWeights: (n_cells × n_lineages)
        conditions: (n_cells,) labels
        n_lineages: optional override

    Returns:
        dict with keys: pvalue, statistic, df, contingency_table
    """
    cw = np.asarray(cellWeights, dtype=np.float64)
    if cw.ndim == 1: cw = cw[:, None]
    # Each cell's dominant lineage
    dominant = np.argmax(cw, axis=1)
    cond = np.asarray(conditions)
    # Contingency: rows = conditions, cols = lineages
    levels = np.unique(cond)
    n_lin = cw.shape[1] if n_lineages is None else n_lineages
    table = np.zeros((len(levels), n_lin), dtype=int)
    for i, lev in enumerate(levels):
        for L in range(n_lin):
            table[i, L] = int(((cond == lev) & (dominant == L)).sum())
    # Avoid empty rows / cols
    if table.sum() == 0 or len(levels) < 2 or n_lin < 2:
        return {"pvalue": 1.0, "statistic": np.nan, "df": 0,
                "contingency_table": table}
    chi2_stat, p, dof, _ = chi2_contingency(table + 1e-9)  # tiny smoothing
    return {"pvalue": float(p), "statistic": float(chi2_stat), "df": int(dof),
            "contingency_table": table}


def create_differential_topology(
    n_cells: int = 1000,
    n_conditions: int = 2,
    shift: float = 0.5,
    seed: int = 42,
) -> dict:
    """Test-data generator — synthetic Slingshot-like fixture with conditions.

    Mirrors condiments::create_differential_topology but in a simpler form
    (we don't shell out to Slingshot; we generate pseudotime + cellWeights
    directly).
    """
    rng = np.random.default_rng(seed)
    cond = rng.integers(0, n_conditions, n_cells)
    # 2 lineages
    pt = rng.uniform(0, 1, (n_cells, 2))
    cw = np.zeros((n_cells, 2))
    for i in range(n_cells):
        # Lineage choice depends on condition
        p_lin1 = 0.5 + shift * (cond[i] - 0.5)
        if rng.random() < p_lin1:
            cw[i, 0] = 1.0; pt[i, 1] = np.nan
        else:
            cw[i, 1] = 1.0; pt[i, 0] = np.nan
    return {"pseudotime": pt, "cellWeights": cw,
            "conditions": cond.astype(str), "n_lineages": 2}
