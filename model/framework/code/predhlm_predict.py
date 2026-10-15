"""Thin wrapper around the PredHLM ensemble (Jang et al., bioRxiv 2026).

Mirrors predict.py of joshua-416/PredHLM (RDKit2D + XGBoost, 5 folds) without
changing the prediction path. Differences are confined to input handling:
molecules that cannot be parsed, or that give a non-finite NaN in a descriptor
the model uses, return NaN instead of aborting the whole batch.
"""
import os

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.logger().setLevel(RDLogger.CRITICAL)

ROOT = os.path.dirname(os.path.abspath(__file__))
CKPT = os.path.join(ROOT, "..", "..", "checkpoints")
N_FOLDS = 5
COLUMNS = ["pred_log10_half_life", "pred_std", "pred_half_life_min", "sdc"]

_state = {}


def _descriptor_functions():
    """RDKit 2D descriptors in the order used by training (Ipc excluded, SPS de-duplicated)."""
    seen, out = set(), {}
    for name, fn in Descriptors._descList:
        if name == "Ipc" or name in seen:
            continue
        seen.add(name)
        out[name] = fn
    return out


def _ecfp4(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def _load():
    if _state:
        return _state
    names = pd.read_csv(os.path.join(CKPT, "rdkit2d_names.csv")).iloc[:, 0].tolist()
    funcs = _descriptor_functions()
    _state["fns"] = [funcs[n] for n in names]
    _state["models"] = [
        (
            joblib.load(os.path.join(CKPT, "fold_%d" % i, "scaler.pkl")),
            joblib.load(os.path.join(CKPT, "fold_%d" % i, "model.pkl")),
        )
        for i in range(N_FOLDS)
    ]
    train = pd.read_csv(os.path.join(CKPT, "training_reference.csv"))["smiles"].astype(str)
    fps = []
    for s in train:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            fps.append(_ecfp4(mol))
    _state["train_fps"] = fps
    return _state


def _sdc(mol, train_fps):
    sims = np.asarray(DataStructs.BulkTanimotoSimilarity(_ecfp4(mol), train_fps), dtype=np.float64)
    td = 1.0 - sims
    with np.errstate(divide="ignore", invalid="ignore"):
        expo = np.where(td < 1.0, -3.0 * td / (1.0 - td), -np.inf)
    return float(np.sum(np.exp(expo)))


def predict(smiles_list):
    """Return an (n, 4) float array: log10 t1/2, ensemble SD, t1/2 (min), SDC (NaN rows for failures)."""
    st = _load()
    out = np.full((len(smiles_list), len(COLUMNS)), np.nan, dtype=np.float64)
    rows, idx, mols = [], [], []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(str(smi)) if isinstance(smi, str) and smi.strip() else None
        if mol is None:
            continue
        try:
            vec = [float(fn(mol)) for fn in st["fns"]]
        except Exception:
            continue
        vec = np.asarray(vec, dtype=np.float64)
        if np.isnan(vec).any():
            continue
        # predict.py maps +-inf -> NaN -> 0 after column selection
        vec[np.isinf(vec)] = 0.0
        rows.append(vec)
        idx.append(i)
        mols.append(mol)
    if not rows:
        return out
    X = np.stack(rows)
    preds = np.stack([m.predict(sc.transform(X)) for sc, m in st["models"]])
    mean, std = preds.mean(axis=0), preds.std(axis=0)
    for j, i in enumerate(idx):
        out[i, 0] = mean[j]
        out[i, 1] = std[j]
        out[i, 2] = 10.0 ** mean[j]
        out[i, 3] = _sdc(mols[j], st["train_fps"])
    return out
