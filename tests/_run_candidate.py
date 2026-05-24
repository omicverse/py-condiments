"""Candidate runner — loads R sidecars + runs pycondiments + dumps JSON."""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))
import pycondiments


def main():
    fixture_path, output_path = sys.argv[1], sys.argv[2]
    pp = Path(fixture_path)
    rd = pd.read_csv(pp.with_name(pp.stem + "_reducedDim.csv")).to_numpy(dtype=np.float64)
    pt = pd.read_csv(pp.with_name(pp.stem + "_pseudotime.csv")).to_numpy(dtype=np.float64)
    cw = pd.read_csv(pp.with_name(pp.stem + "_cellWeights.csv")).to_numpy(dtype=np.float64)
    cond = pd.read_csv(pp.with_name(pp.stem + "_conditions.csv"))["condition"].values
    print(f"[cand] rd {rd.shape}, pt {pt.shape}, cw {cw.shape}, cond {cond.shape}")

    np.random.seed(42)
    imb = pycondiments.imbalance_score(rd, cond, k=10, smooth=10)
    np.random.seed(42)
    ptest = pycondiments.progressionTest(pt, cond, cellWeights=cw)
    np.random.seed(42)
    ttest = pycondiments.topologyTest(pt, cw, cond)

    out = {
        "imbalance": {
            "score": imb["score"].tolist(),
            "scaled_score": imb["scaled_score"].tolist(),
        },
        "topology": {"pvalue": float(ttest["pvalue"])},
        "progression": {
            "pvalue": float(ptest["pvalue"].iloc[0]),
            "statistic": float(ptest["statistic"].iloc[0]),
        },
        "differentiation": {"pvalue": None},
    }
    with open(output_path, "w") as f:
        json.dump(out, f)
    print(f"[cand] wrote {output_path}")


if __name__ == "__main__":
    main()
