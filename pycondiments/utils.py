"""Utility functions: merge_sds."""

from __future__ import annotations

import numpy as np


def merge_sds(*sds_list) -> dict:
    """Merge multiple SlingshotDataSet-like objects.

    Each `sds` is a dict with keys `pseudotime` and `cellWeights`. Concatenates
    them row-wise. Useful when fitting Slingshot per condition then needing a
    unified object for downstream tests.
    """
    if not sds_list:
        raise ValueError("merge_sds: no inputs")
    pts = []
    cws = []
    for s in sds_list:
        pt = np.asarray(s.get("pseudotime"))
        cw = np.asarray(s.get("cellWeights"))
        if pt.ndim == 1: pt = pt[:, None]
        if cw.ndim == 1: cw = cw[:, None]
        pts.append(pt); cws.append(cw)
    # Pad to common n_lineages
    n_lin = max(p.shape[1] for p in pts)
    def pad(a, n_lin):
        if a.shape[1] < n_lin:
            pad = np.full((a.shape[0], n_lin - a.shape[1]), np.nan)
            return np.hstack([a, pad])
        return a
    pts_padded = [pad(p, n_lin) for p in pts]
    cws_padded = [pad(p, n_lin) for p in cws]
    return {"pseudotime": np.vstack(pts_padded),
            "cellWeights": np.vstack(cws_padded)}
