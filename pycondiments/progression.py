"""progressionTest, differentiationTest, weights_from_pst.

Mirrors condiments::progressionTest (KS test on pseudotime between conditions)
+ differentiationTest (KS test on lineage weights between conditions).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, chi2


def weights_from_pst(pseudotime: np.ndarray) -> np.ndarray:
    """Convert raw pseudotime → cell weights per lineage (mirrors condiments::weights_from_pst).

    For NaN pseudotime, weight = 0 (cell didn't follow that lineage).
    For non-NaN, weight = 1.

    Args:
        pseudotime: (n_cells × n_lineages) — NA-filled where a cell doesn't
                    participate in a lineage (Slingshot convention)

    Returns:
        (n_cells × n_lineages) 0/1 cellWeights
    """
    pt = np.asarray(pseudotime, dtype=np.float64)
    if pt.ndim == 1:
        pt = pt[:, None]
    return (~np.isnan(pt)).astype(np.float64)


def progressionTest(
    pseudotime: np.ndarray,
    conditions: np.ndarray,
    *,
    cellWeights: np.ndarray | None = None,
    global_test: bool = True,
    lineages: bool = False,
) -> pd.DataFrame:
    """Do cells advance along pseudotime at different rates per condition?

    Per-lineage KS test on pseudotime values between conditions.

    Args:
        pseudotime: (n_cells × n_lineages)
        conditions: (n_cells,) condition labels (2 levels supported in v0.1)
        cellWeights: optional (n_cells × n_lineages) — only cells with
                     weight > 0 are included in each lineage's test
        global_test: if True, combine per-lineage tests via Fisher's method
        lineages: if True, return per-lineage results instead of just global

    Returns:
        DataFrame with columns:
            statistic, pvalue (per lineage if lineages=True, else single row)
    """
    pt = np.asarray(pseudotime, dtype=np.float64)
    if pt.ndim == 1:
        pt = pt[:, None]
    cond = np.asarray(conditions)
    levels = np.unique(cond)
    if len(levels) != 2:
        raise NotImplementedError("v0.1 supports 2-condition tests only")
    cw = cellWeights if cellWeights is not None else weights_from_pst(pt)
    cw = np.asarray(cw)
    if cw.ndim == 1: cw = cw[:, None]

    n_lin = pt.shape[1]
    per_lineage = []
    for L in range(n_lin):
        mask = cw[:, L] > 0
        if mask.sum() < 4:
            per_lineage.append((np.nan, 1.0))
            continue
        x = pt[mask & (cond == levels[0]), L]
        y = pt[mask & (cond == levels[1]), L]
        x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
        if len(x) < 2 or len(y) < 2:
            per_lineage.append((np.nan, 1.0))
            continue
        stat, p = ks_2samp(x, y)
        per_lineage.append((float(stat), float(p)))

    if lineages:
        return pd.DataFrame(per_lineage, columns=["statistic", "pvalue"],
                            index=[f"lineage_{i+1}" for i in range(n_lin)])
    if global_test:
        # Fisher's combined p-value
        pvals = [p for _, p in per_lineage if np.isfinite(p) and p > 0]
        if not pvals:
            return pd.DataFrame([{"statistic": np.nan, "pvalue": 1.0}])
        chi_sq = -2 * sum(np.log(p) for p in pvals)
        df = 2 * len(pvals)
        p_global = float(chi2.sf(chi_sq, df))
        return pd.DataFrame([{"statistic": float(chi_sq), "pvalue": p_global}])
    return pd.DataFrame(per_lineage, columns=["statistic", "pvalue"],
                        index=[f"lineage_{i+1}" for i in range(n_lin)])


def differentiationTest(
    cellWeights: np.ndarray,
    conditions: np.ndarray,
    *,
    global_test: bool = True,
) -> pd.DataFrame:
    """Do cells make different lineage choices per condition?

    Compare lineage-weight distributions between conditions via KS test on
    each cell's lineage-1-fraction (for 2 lineages) or χ² on assigned lineage
    (multi-lineage).

    Args:
        cellWeights: (n_cells × n_lineages) — normalised lineage weights
        conditions: (n_cells,) labels
        global_test: combine across lineages

    Returns:
        DataFrame with statistic, pvalue.
    """
    cw = np.asarray(cellWeights, dtype=np.float64)
    if cw.ndim == 1: cw = cw[:, None]
    cond = np.asarray(conditions)
    levels = np.unique(cond)
    if len(levels) != 2:
        raise NotImplementedError("v0.1 supports 2-condition tests only")
    # Use the first lineage's weight as the test statistic per cell
    n_lin = cw.shape[1]
    per_lineage = []
    for L in range(n_lin):
        x = cw[cond == levels[0], L]
        y = cw[cond == levels[1], L]
        x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
        if len(x) < 2 or len(y) < 2:
            per_lineage.append((np.nan, 1.0))
            continue
        stat, p = ks_2samp(x, y)
        per_lineage.append((float(stat), float(p)))

    if global_test:
        pvals = [p for _, p in per_lineage if np.isfinite(p) and p > 0]
        if not pvals:
            return pd.DataFrame([{"statistic": np.nan, "pvalue": 1.0}])
        chi_sq = -2 * sum(np.log(p) for p in pvals)
        df = 2 * len(pvals)
        p_global = float(chi2.sf(chi_sq, df))
        return pd.DataFrame([{"statistic": float(chi_sq), "pvalue": p_global}])
    return pd.DataFrame(per_lineage, columns=["statistic", "pvalue"],
                        index=[f"lineage_{i+1}" for i in range(n_lin)])
