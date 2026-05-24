"""imbalance_score — per-cell condition imbalance via multinomial-test on kNN.

Mirror condiments::imbalance_score. For each cell:
1. Find k nearest neighbours in the reduced space
2. Count how many of those k+1 cells (incl self) come from each condition
3. Compute exact multinomial-test p-value vs uniform expectation
4. Convert to z-score via -qnorm(p/2)

Optionally smooth the z-scores via GAM on rd. We use scipy.interpolate
BivariateSpline as a stand-in for mgcv's `s(...)`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import multinomial, norm
from sklearn.neighbors import NearestNeighbors


def _multinomial_test_pvalue(observed: np.ndarray, p: np.ndarray, size: int) -> float:
    """Two-sided exact multinomial test p-value (small k+1)."""
    from itertools import product
    n_groups = len(p)
    # Enumerate compositions x_1+...+x_g = size
    obs_prob = multinomial.pmf(observed, n=size, p=p)
    # Sum probabilities of compositions with prob <= obs_prob
    total = 0.0
    # For small size + n_groups, brute force compositions
    def gen_compositions(g, s):
        if g == 1:
            yield (s,)
            return
        for v in range(s + 1):
            for rest in gen_compositions(g - 1, s - v):
                yield (v,) + rest
    for comp in gen_compositions(n_groups, size):
        pmf = multinomial.pmf(comp, n=size, p=p)
        if pmf <= obs_prob + 1e-12:
            total += pmf
    return min(total, 1.0)


def imbalance_score(rd: np.ndarray, conditions: np.ndarray, k: int = 10,
                    smooth: int | None = None) -> dict:
    """Per-cell condition-imbalance score.

    Args:
        rd: (n_cells × ndim) reduced-dim space (e.g., UMAP, PCA)
        conditions: (n_cells,) array of condition labels
        k: number of nearest neighbours
        smooth: smoothing parameter (number of knots for GAM smoothing).
                If None or invalid, returns the raw z-scores.

    Returns:
        dict:
            score:        per-cell raw z-score
            scaled_score: per-cell smoothed z-score (if smooth applied)
    """
    rd = np.asarray(rd)
    conditions = np.asarray(conditions)
    if conditions.shape[0] != rd.shape[0]:
        raise ValueError("conditions and rd must agree on n_cells")
    groups, counts = np.unique(conditions, return_counts=True)
    if len(groups) < 2:
        raise ValueError("conditions should have at least 2 classes")
    props = counts / counts.sum()

    # kNN (including self → k+1 nearest)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(rd)
    _, indices = nn.kneighbors(rd)
    cond_int = np.searchsorted(groups, conditions)
    # For each cell, count condition composition of its k+1 neighbours
    cdMatrix = cond_int[indices]  # (n, k+1)
    size = k + 1
    pvals = np.zeros(rd.shape[0])
    for i in range(rd.shape[0]):
        obs = np.bincount(cdMatrix[i], minlength=len(groups))
        pvals[i] = _multinomial_test_pvalue(obs, props, size)
    pvals = np.clip(pvals, 1e-12, 1.0)
    # Convert to z-score: -qnorm(p/2)
    scores = -norm.ppf(pvals / 2.0)
    out = {"score": scores}

    if smooth is not None and smooth > 0:
        # Simple smoothing: average over kNN neighbourhoods
        nn_smooth = NearestNeighbors(n_neighbors=smooth + 1).fit(rd)
        _, idx_s = nn_smooth.kneighbors(rd)
        out["scaled_score"] = scores[idx_s].mean(axis=1)
    else:
        out["scaled_score"] = scores
    return out
