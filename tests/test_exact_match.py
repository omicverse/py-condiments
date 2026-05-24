"""Parity gate against R condiments on canonical fixture."""
import json, sys
from pathlib import Path
import numpy as np
import pytest
import yaml

PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PORT))
sys.path.insert(0, str(PORT.parent / "omicverse-rebuildr" / "engine"))
from parity_metrics import compute_parity


@pytest.fixture(scope="session")
def manifest():
    return yaml.safe_load((PORT / "data" / "manifest.yaml").read_text())


@pytest.fixture(scope="session")
def outputs():
    ref = PORT / "data" / "reference_output.json"
    cand = PORT / "data" / "candidate_output.json"
    if not (ref.exists() and cand.exists()):
        pytest.skip("Run r_reference_driver.R + _run_candidate.py first")
    return json.loads(ref.read_text()), json.loads(cand.read_text())


def test_imbalance_score_parity(manifest, outputs):
    r, p = outputs
    r_imb = np.array(r["imbalance"]["score"]); p_imb = np.array(p["imbalance"]["score"])
    mask = np.isfinite(r_imb) & np.isfinite(p_imb)
    m = compute_parity(r_imb[mask], p_imb[mask], "ordinal")
    spec = next(o for o in manifest["outputs"] if o["name"] == "imbalance_score")
    assert m >= spec["threshold"], f"imbalance Pearson {m:.4f} < {spec['threshold']}"
