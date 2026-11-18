# Archive

Superseded exploratory scripts, kept for development-history reference
rather than deleted - these informed decisions documented in
`METHODS.md` (e.g. why rank-based confidence normalization was adopted
over absolute thresholds or min-max scaling) but are not part of the
current, documented analysis pipeline.

- `normalized_confidence_comparison.py`, `rank_based_confidence_comparison.py`,
  `three_tool_confidence_comparison.py` - the three normalization
  approaches tested before adopting cluster purity (see METHODS.md
  Section 7.3).
- `combination_search_figure.py` - an earlier, simpler version of
  `combination_justification_figure.py` (mean score vs. number of
  methods only, without the ambiguous-cluster-count dimension that made
  the "what did the winner miss" argument concrete).
- `fig1_two_step_annotation.py`, `two_step_figure_generic.py` - the
  original two-step annotation approach (liver atlas, then immune-atlas
  refinement for low-confidence clusters), superseded once the immune
  atlas was found to decrease accuracy (see `immune_atlas_sanity_check.py`)
  and the project moved to the current weighted-consensus design instead.
