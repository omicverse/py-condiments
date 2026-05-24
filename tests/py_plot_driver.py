"""Render py-condiments plots on the same R-dumped fixture."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))

import pycondiments
from pycondiments.plotting import plotConditionsOnEmbedding, plotImbalanceScore
from ggplot2_py import ggsave


def main():
    out_dir = Path(sys.argv[1])
    df = pd.read_csv(out_dir / "df.csv")
    coords = df[["Dim1", "Dim2"]].to_numpy(dtype=np.float64)
    conditions = df["conditions"].tolist()

    p = plotConditionsOnEmbedding(coords, conditions)
    ggsave(str(out_dir / "Py_embedding.png"), plot=p, width=6, height=4, dpi=100)

    # Recompute scaled imbalance scores in Py
    sc = pycondiments.imbalance_score(coords, conditions=conditions)
    p2 = plotImbalanceScore(coords, scores=sc["scaled_score"])
    ggsave(str(out_dir / "Py_imbalance.png"), plot=p2, width=6, height=4, dpi=100)
    print("done")


if __name__ == "__main__":
    main()
