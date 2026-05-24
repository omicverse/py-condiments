import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pycondiments


def test_import():
    assert pycondiments.__version__.startswith("0.1")
    for fn in ("imbalance_score", "topologyTest", "progressionTest",
               "differentiationTest", "weights_from_pst", "merge_sds",
               "create_differential_topology"):
        assert hasattr(pycondiments, fn)


def test_create_differential_topology():
    toy = pycondiments.create_differential_topology(n_cells=200, shift=0.5, seed=42)
    assert toy["pseudotime"].shape == (200, 2)
    assert toy["cellWeights"].shape == (200, 2)
    assert len(toy["conditions"]) == 200


def test_imbalance_score():
    rng = np.random.default_rng(42)
    rd = rng.normal(0, 1, (200, 2))
    # Condition A clustered around origin, B around (3, 3) for clear imbalance
    conds = np.array(['A'] * 100 + ['B'] * 100)
    rd[100:] += 3.0
    out = pycondiments.imbalance_score(rd, conds, k=10)
    assert "score" in out and "scaled_score" in out
    assert out["score"].shape == (200,)
    # Strong imbalance → high mean score
    assert np.nanmean(out["score"]) > 0.5


def test_progressionTest():
    toy = pycondiments.create_differential_topology(n_cells=300, shift=0.5, seed=42)
    res = pycondiments.progressionTest(toy["pseudotime"], toy["conditions"],
                                       cellWeights=toy["cellWeights"])
    assert "pvalue" in res.columns
    assert 0 <= res["pvalue"].iloc[0] <= 1


def test_topologyTest():
    toy = pycondiments.create_differential_topology(n_cells=300, shift=0.5, seed=42)
    res = pycondiments.topologyTest(toy["pseudotime"], toy["cellWeights"],
                                    toy["conditions"])
    assert "pvalue" in res
    assert 0 <= res["pvalue"] <= 1


def test_weights_from_pst():
    pt = np.array([[1.0, np.nan], [np.nan, 2.0], [3.0, np.nan]])
    cw = pycondiments.weights_from_pst(pt)
    assert (cw == np.array([[1, 0], [0, 1], [1, 0]])).all()


def test_merge_sds():
    a = {"pseudotime": np.array([[0.1, np.nan], [0.2, np.nan]]),
         "cellWeights": np.array([[1.0, 0.0], [1.0, 0.0]])}
    b = {"pseudotime": np.array([[np.nan, 0.3], [np.nan, 0.4]]),
         "cellWeights": np.array([[0.0, 1.0], [0.0, 1.0]])}
    merged = pycondiments.merge_sds(a, b)
    assert merged["pseudotime"].shape == (4, 2)
