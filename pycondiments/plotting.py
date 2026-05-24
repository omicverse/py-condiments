"""Visualization helpers — 1:1 port of condiments vignette plots using ggplot2-python.

condiments doesn't export named ``plot.*`` functions; the package's vignettes
build plots inline with ggplot2 calls. We package the most-common patterns
here so users can reproduce the vignette plots with one call.

All functions return a ``ggplot2_py.GGPlot``. Render with::

    from ggplot2_py import ggsave
    ggsave("out.png", plot=p, width=6, height=4, dpi=120)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ggplot2_py import (
    aes,
    geom_density,
    geom_path,
    geom_point,
    ggplot,
    labs,
    scale_color_viridis_c,
    theme_classic,
)


def plotConditionsOnEmbedding(
    coords: np.ndarray | pd.DataFrame,
    conditions,
    curves: pd.DataFrame | None = None,
    alpha: float = 0.5,
    point_size: float = 1.5,
):
    """Scatter of 2-D embedding coloured by condition, with optional curve overlay.

    Args:
        coords: (n × 2) array — rows are cells.
        conditions: per-cell categorical labels.
        curves: optional DataFrame with columns ``Dim1``, ``Dim2``, ``lineages``,
            and ``Order`` (the output of ``slingCurves(sds, as.df=TRUE)``).
        alpha, point_size: point aesthetics.
    """
    arr = np.asarray(coords, dtype=np.float64)[:, :2]
    df = pd.DataFrame(arr, columns=["Dim1", "Dim2"])
    df["conditions"] = list(conditions)

    p = (
        ggplot(df, aes(x="Dim1", y="Dim2", colour="conditions"))
        + geom_point(alpha=alpha, size=point_size)
        + theme_classic()
        + labs(x="Dim1", y="Dim2", colour="Condition")
    )

    if curves is not None:
        curves = curves.sort_values(["lineages", "Order"]) if "Order" in curves.columns else curves
        p = p + geom_path(
            aes(x="Dim1", y="Dim2", group="lineages"),
            data=curves,
            size=1.2,
            colour="black",
        )

    return p


def plotImbalanceScore(
    coords: np.ndarray | pd.DataFrame,
    scores: np.ndarray,
    title: str | None = None,
    point_size: float = 1.5,
):
    """Scatter of 2-D embedding coloured by condition-imbalance score."""
    arr = np.asarray(coords, dtype=np.float64)[:, :2]
    df = pd.DataFrame(arr, columns=["Dim1", "Dim2"])
    df["score"] = np.asarray(scores, dtype=np.float64)
    p = (
        ggplot(df, aes(x="Dim1", y="Dim2", colour="score"))
        + geom_point(size=point_size)
        + theme_classic()
        + scale_color_viridis_c()
        + labs(x="Dim1", y="Dim2", colour="Imbalance score")
    )
    if title is not None:
        p = p + labs(title=title)
    return p


def plotPseudotimeByCondition(
    pseudotime,
    conditions,
    alpha: float = 0.5,
):
    """Density plot of pseudotime stratified by condition."""
    df = pd.DataFrame(
        {
            "pseudotime": np.asarray(pseudotime, dtype=np.float64),
            "conditions": list(conditions),
        }
    )
    return (
        ggplot(df, aes(x="pseudotime", fill="conditions"))
        + geom_density(alpha=alpha)
        + theme_classic()
        + labs(x="Pseudotime", y="Density", fill="Condition")
    )


def plotWeightByCondition(
    weights: np.ndarray,
    conditions,
    lineage_index: int = 0,
    alpha: float = 0.5,
):
    """Density plot of cell-weights for a chosen lineage, stratified by condition."""
    w = np.asarray(weights, dtype=np.float64)
    w_col = w[:, lineage_index] if w.ndim == 2 else w
    df = pd.DataFrame(
        {
            "weight": w_col,
            "conditions": list(conditions),
        }
    )
    return (
        ggplot(df, aes(x="weight", fill="conditions"))
        + geom_density(alpha=alpha)
        + theme_classic()
        + labs(
            x=f"Lineage {lineage_index + 1} weight",
            y="Density",
            fill="Condition",
        )
    )
