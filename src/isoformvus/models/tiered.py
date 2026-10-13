"""Tiered interpretable models.

Five model tiers share the *same* learner, splits, and variants; they differ
only in which feature tiers are visible (see features/spec.md):

    M1  variant only
    M2  + phenotype
    M3  + cell type
    M4  + isoform            (IsoformVUS)
    M5  + isoform effect     (IsoformVUS+)

The M4-vs-M3 (and M5-vs-M4) comparison is the critical test: any delta is
attributable to isoform information, not to data or modeling choices.

Learners (in order of interpretability):
    1. L2 logistic regression (coefficients directly interpretable)
    2. LightGBM (gain/SHAP importance; still auditable)

A neural network is deliberately NOT included: with ~10^3-10^4 variants and
gene-level splits, an MLP adds opacity without a justified capacity need.
It may be revisited only if linear+GBM both underfit the validation set.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

# Feature tiers in column-name terms. Concrete column lists are built by
# features/build.py from spec.md; these are the tier *prefixes* each model sees.
TIER_COLUMNS: dict[str, list[str]] = {
    # Exact column names are registered in features/build.py FEATURE_TIERS.
    "M1": ["variant"],
    "M2": ["variant", "phenotype"],
    "M3": ["variant", "phenotype", "celltype"],
    "M4": ["variant", "phenotype", "celltype", "isoform"],
    "M5": ["variant", "phenotype", "celltype", "isoform", "isoform_effect"],
}

MODEL_IDS = ["M1", "M2", "M3", "M4", "M5"]


@dataclass
class TieredResult:
    model_id: str
    learner: str
    feature_columns: list[str] = field(default_factory=list)
    val_auroc: float = float("nan")
    test_scores: np.ndarray = field(default_factory=lambda: np.array([]))


def make_learner(name: str, seed: int = 20261006):
    """Instantiate an interpretable learner by name."""
    if name == "logreg":
        return LogisticRegression(C=1.0, max_iter=5000,
                                  class_weight="balanced",
                                  random_state=seed)
    if name == "lightgbm":
        from lightgbm import LGBMClassifier
        return LGBMClassifier(n_estimators=1000, learning_rate=0.03,
                              num_leaves=31, min_child_samples=50,
                              subsample=0.8, subsample_freq=1,
                              colsample_bytree=0.8, reg_lambda=1.0,
                              class_weight="balanced", random_state=seed,
                              verbose=-1)
    raise ValueError(f"Unknown learner: {name!r}")


def fit_tiered_resolved(X_train: pd.DataFrame, y_train: pd.Series,
                        X_val: pd.DataFrame, y_val: pd.Series,
                        resolved: dict[str, list[str]],
                        learner_name: str = "logreg",
                        seed: int = 20261006) -> dict[str, TieredResult]:
    """Fit tiers given an explicit model_id -> column list mapping."""
    results: dict[str, TieredResult] = {}
    for model_id in MODEL_IDS:
        cols = resolved[model_id]
        missing = [c for c in cols if c not in X_train.columns]
        if missing:
            raise KeyError(f"{model_id}: missing columns {missing}")
        # Strict nesting check: each tier must extend the previous one.
        idx = MODEL_IDS.index(model_id)
        if idx > 0:
            prev = set(resolved[MODEL_IDS[idx - 1]])
            if not prev.issubset(set(cols)):
                raise ValueError(f"{model_id} is not a strict superset of "
                                 f"{MODEL_IDS[idx - 1]}; ablation invalid")
        clf = make_learner(learner_name, seed)
        clf.fit(X_train[cols], y_train)
        proba = clf.predict_proba(X_val[cols])[:, 1]
        results[model_id] = TieredResult(
            model_id=model_id, learner=learner_name,
            feature_columns=list(cols),
            val_auroc=float(roc_auc_score(y_val, proba)))
    return results
