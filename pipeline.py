"""
pipeline.py
Orchestrator: reads the method registry (config.py), runs each configured
method that hasn't already been run, merges results into adata.obs, then
computes the weighted consensus vote.

This replaces the earlier workflow of manually running separate scripts
(singler_immgen.R, xiong_add_immgen_labels.py, xiong_celltypist_dc_subtypes.py,
consensus_scoring.py) by hand in sequence with hardcoded paths. Adding a
new tool now means adding one MethodConfig entry (config.py) and, if it's
a genuinely new kind of tool, one runner function - nothing else in this
file changes.
"""

import scanpy as sc

from config import build_registry
from consensus import compute_consensus
from runners.r_runner import run_r_method


def ensure_normalized(adata):
    """CellTypist (and typically SingleR) expect log1p-normalised-to-10k
    expression in .X. Ensure that's the case before running any method.
    """
    if "log1p_norm" in adata.layers:
        adata.X = adata.layers["log1p_norm"]
    else:
        print("NOTE: assuming adata.X is already log1p-normalised. "
              "If a method raises a normalisation error below, check this.")
    return adata


def run_pipeline(h5ad_path: str, output_path: str, tissue: str = "liver",
                  skip_existing: bool = True, selected_methods: list = None):
    """
    Parameters
    ----------
    selected_methods : list[str], optional
        Method NAMES (matching MethodConfig.name in config.py) to actually
        run and include in the consensus vote - e.g. ["ct_label",
        "scanvi_label", "singler_label"]. If None (the default), uses the
        full tissue registry, same as before.

        This exists because the combination_search.py analysis showed the
        single BEST-scoring method combination differs depending on the
        dataset (e.g. the Xiong vs. GSE166504 rankings didn't match) - so
        rather than hardcoding one "best" set, the full registry is kept
        as a sensible DEFAULT, while a user who has done their own
        combination analysis on their specific data can override it here.
    """
    print(f"Tissue: {tissue}")
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    full_registry = build_registry(tissue)

    if selected_methods is not None:
        by_name = {m.name: m for m in full_registry}
        unknown = [n for n in selected_methods if n not in by_name]
        if unknown:
            raise ValueError(
                f"Unknown method name(s): {unknown}. "
                f"Valid method names for tissue '{tissue}': {list(by_name.keys())}"
            )
        methods = [by_name[n] for n in selected_methods]
        print(f"Using USER-SELECTED methods ({len(methods)}/{len(full_registry)}): "
              f"{[m.name for m in methods]}")
    else:
        methods = full_registry
        print(f"Using DEFAULT full registry ({len(methods)} methods): "
              f"{[m.name for m in methods]}")

    adata = ensure_normalized(adata)

    for m in methods:
        if skip_existing and m.name in adata.obs.columns:
            print(f"\n[{m.name}] already present in obs - skipping (skip_existing=True)")
            continue

        print(f"\n[{m.name}] running ({m.kind})...")

        if m.runner is None:
            print(f"  WARNING: no runner configured for '{m.name}' yet - skipping. "
                  f"(This is expected for methods still under development, "
                  f"e.g. scANVI / the future mouse-native immune reference.)")
            continue

        if m.kind == "python":
            labels, confidence = m.runner(adata, **m.kwargs)
        elif m.kind == "r":
            labels, confidence = run_r_method(adata, script_path=m.runner, **m.kwargs)
        else:
            raise ValueError(f"Unknown method kind: {m.kind}")

        adata.obs[m.name] = labels.reindex(adata.obs_names)
        adata.obs[f"{m.name}_conf"] = confidence.reindex(adata.obs_names)

    print("\nComputing consensus...")
    adata = compute_consensus(adata, methods)

    print(f"\nSaving to {output_path}...")
    adata.write_h5ad(output_path)
    print("Done!")
    return adata


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tissue", default="liver",
                         help="Which tissue's method registry to use (liver, adipose, brain, ...). "
                              "Defaults to 'liver'.")
    parser.add_argument("--methods", default=None,
                         help="Comma-separated list of method names to use instead of the "
                              "full default registry, e.g. "
                              "--methods ct_label,scanvi_label,singler_label. "
                              "If omitted, uses the full tissue registry (the default). "
                              "The 'best' combination has been shown to differ by dataset "
                              "(see combination_search.py) - use that script on your own "
                              "data first if you want an informed choice here, rather than "
                              "guessing.")
    args = parser.parse_args()
    selected = args.methods.split(",") if args.methods else None
    run_pipeline(args.input, args.output, tissue=args.tissue, selected_methods=selected)