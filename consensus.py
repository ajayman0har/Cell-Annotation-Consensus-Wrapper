"""
consensus.py
Generic weighted-vote consensus scorer.

MATH (documented explicitly - this is what to show a PI/reviewer):

  Every method's "confidence" on a cluster is measured the SAME way -
  CLUSTER PURITY - regardless of what that method's own native output
  looks like:

      purity(m, c) = (# cells in cluster c where method m's harmonized
                       label equals m's own MAJORITY label for c)
                     / (# cells in cluster c where m cast a real,
                        non-masked, non-placeholder label)

  This replaces an earlier approach (percentile-ranking each method's
  native confidence score within its own distribution) that was needed
  because the six methods' native outputs are NOT the same kind of
  number: CellTypist/scANVI/SCINA report real classifier probabilities
  (bounded 0-1), SingleR reports a correlation/similarity score (also
  0-1-ish but answering a different question), and scType reports an
  unbounded raw marker-enrichment sum (can exceed 20). Percentile-ranking
  made these comparable in RELATIVE terms (this cluster vs. this same
  method's other clusters) but never made them comparable to EACH OTHER
  in absolute terms - purity(m, c) = 0.9 always means the same thing
  ("90% of this cluster's cells got the same raw label from method m"),
  regardless of which method m is. No per-method normalization needed.

  Each method also has a fixed TRUST WEIGHT w_m (set in config.py,
  reflecting how much we trust that method/reference combination overall
  - e.g. the DC-subtype refinement is weighted slightly below a full
  vote, since it's a narrow specialist, not a general-purpose classifier).

  For cell i (belonging to cluster c), each eligible method m (not masked
  out, not a placeholder/no-call value) casts its harmonized label with
  weight:

      vote_weight(i, m) = w_m * purity(m, c)

  votes(i, label) = sum of vote_weight(i, m) over methods m voting for
                    that label

  consensus_label(i)  = argmax_label votes(i, label)              (plurality)
  consensus_score(i)  = votes(i, consensus_label) / sum(vote_weight(i, m)
                         over all eligible methods for cell i)
  consensus_n_methods(i) = count of methods that cast a real vote for cell i
"""

from collections import defaultdict
import numpy as np
import pandas as pd

from harmonize import harmonize_series, check_for_placeholders


def _compute_cluster_purity(harmonized_col: pd.Series, cluster_values: np.ndarray) -> dict:
    """
    Computes, for each cluster, this method's own cluster purity: the
    fraction of the cluster's cells (that this method actually labeled -
    excluding masked-out/placeholder cells) which agree with the
    method's own majority label for that cluster.

    Uses the HARMONIZED label (not the raw string), so purity reflects
    agreement on the shared category the vote actually happens in, not
    incidental raw-vocabulary differences.

    Returns a dict {cluster_value: purity}. A cluster where this method
    never cast a real label at all (fully masked out) gets NaN.
    """
    df = pd.DataFrame({"cluster": cluster_values, "label": harmonized_col.values})
    purity = {}
    for cl, group in df.groupby("cluster"):
        labels = group["label"].dropna()
        if len(labels) == 0:
            purity[cl] = np.nan
            continue
        top_count = labels.value_counts().iloc[0]
        purity[cl] = top_count / len(labels)
    return purity


def compute_consensus(adata, methods, cluster_col="leiden_r05", verbose=True):
    """
    Parameters
    ----------
    adata : AnnData
        Must have all methods[i].name columns already present in adata.obs,
        plus a cluster column (default 'leiden_r05') used to compute each
        method's per-cluster purity.
    methods : list[MethodConfig]
        The method registry (see config.build_registry()).
    cluster_col : str
        Column in adata.obs identifying cluster membership.

    Returns
    -------
    adata (modified in place) with new obs columns:
        consensus_label, consensus_score, consensus_n_methods
    """
    obs = adata.obs
    label_cols = [m.name for m in methods if m.name in obs.columns]
    missing = [m.name for m in methods if m.name not in obs.columns]
    if missing and verbose:
        print(f"NOTE: {len(missing)} configured methods have no data yet "
              f"and will be skipped: {missing}")

    if cluster_col not in obs.columns:
        raise ValueError(f"cluster_col '{cluster_col}' not found in adata.obs - "
                          f"required to compute per-cluster purity.")

    if verbose:
        placeholder_report = check_for_placeholders(obs, label_cols)
        for col, n in placeholder_report.items():
            print(f"  WARNING: {col} has {n} placeholder/no-call values "
                  f"(e.g. 'pending') - these will be excluded from voting.")

    # ── Apply masks and harmonize each method's labels ────────────────────
    harmonized = pd.DataFrame(index=obs.index)
    base_weights = {}
    purity_by_method = {}
    cluster_values = obs[cluster_col].astype(str).values

    for m in methods:
        if m.name not in obs.columns:
            continue

        col_values = obs[m.name].copy()
        if m.mask_fn is not None:
            eligible = m.mask_fn(adata)
            col_values = col_values.where(eligible, np.nan)
        harmonized[m.name] = harmonize_series(col_values)
        base_weights[m.name] = m.weight

        purity_by_method[m.name] = _compute_cluster_purity(harmonized[m.name], cluster_values)

    vote_cols = list(harmonized.columns)

    # ── Build a per-cell, per-method WEIGHT matrix: w_m * purity(m, cluster(i)) ──
    weight_matrix = np.zeros((len(obs), len(vote_cols)))
    for j, col in enumerate(vote_cols):
        purity_map = purity_by_method[col]
        purity_series = pd.Series(cluster_values).map(purity_map).astype(float).values
        weight_matrix[:, j] = base_weights[col] * purity_series

    values = harmonized.values

    n = len(obs)
    consensus_labels = np.empty(n, dtype=object)
    consensus_scores = np.empty(n, dtype=float)
    n_methods_voted = np.empty(n, dtype=int)

    for i in range(n):
        row = values[i]
        weights_row = weight_matrix[i]
        votes = defaultdict(float)
        n_voted = 0
        total_possible = 0.0
        for label, w in zip(row, weights_row):
            if pd.notna(label) and not np.isnan(w):
                votes[label] += w
                n_voted += 1
                total_possible += w
        if votes:
            best_label = max(votes, key=votes.get)
            score = votes[best_label] / total_possible if total_possible > 0 else 0.0
        else:
            best_label, score = "Unresolved", 0.0
        consensus_labels[i] = best_label
        consensus_scores[i] = score
        n_methods_voted[i] = n_voted

    obs["consensus_label"] = consensus_labels
    obs["consensus_score"] = consensus_scores
    obs["consensus_n_methods"] = n_methods_voted
    adata.obs = obs

    if verbose:
        print("\n-- Consensus label distribution --")
        print(obs["consensus_label"].value_counts())
        print("\n-- Consensus score summary --")
        print(obs["consensus_score"].describe())

    return adata