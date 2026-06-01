# IJMS Computational Robustness Summary

## Threshold sensitivity

Canonical thresholds recovered 7/7 selected genes; across all tested thresholds the minimum recovery was 0/7.

Interpretation: Stable at the canonical and less stringent thresholds; stricter log2FC cutoffs reduce recovery.

## Normalization sensitivity

All normalization methods retained positive severity/primary directions for at least 7/7 selected genes.

Interpretation: Directionally stable across log2CPM, median-ratio, and upper-quartile transforms.

## Covariate-design sensitivity

Selected genes passing the primary DE rule by design: age_sex=5/7; full_available=7/7; unadjusted=6/7.

Interpretation: Most candidates remain directionally positive, but adjusted significance is design-sensitive.

## Leave-one-stage stability

Positive severity direction by omitted stage: NPDR=7/7; NPDR/PDR + DME=7/7; diabetic=7/7; healthy control=7/7. Positive and FDR<0.1: NPDR=7/7; NPDR/PDR + DME=0/7; diabetic=7/7; healthy control=3/7.

Interpretation: Direction is stable; significance weakens when endpoint stages are omitted.
