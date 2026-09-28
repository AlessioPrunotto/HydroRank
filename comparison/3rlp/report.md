# HydraRank–SSTMap comparison: HSP90–3RLP

## Executive result

Under harmonized preprocessing, HydraRank and SSTMap identify nearly the same core hydration-site geometry: all **19 HydraRank sites** have unique SSTMap matches within 1.0 Å, while SSTMap identifies **21 sites**, including two additional lower-occupancy sites. Matched-center separation is 0.103 Å on average and 0.304 Å at maximum.

This is strong evidence of implementation-level spatial concordance for 3RLP. It is not yet scientific validation: it covers one system and the two methods share related inhomogeneous-solvation/nearest-neighbor ideas.

## Quantitative comparison

| Quantity (19 matched sites) | Pearson r | Spearman ρ | HydraRank − SSTMap mean | MAE | RMSE |
|---|---:|---:|---:|---:|---:|
| Occupancy | 0.948 | 0.973 | -0.031 | 0.033 | 0.093 |
| −TΔS translational (kcal/mol) | 0.956 | 0.989 | -0.181 | 0.242 | 0.288 |
| −TΔS orientational (kcal/mol) | 0.978 | 1.000 | 0.097 | 0.187 | 0.246 |
| −TΔS total (kcal/mol) | 0.994 | 0.991 | -0.084 | 0.267 | 0.356 |
| Solute H-bonds (native definitions) | 0.943 | 0.923 | 0.279 | 0.279 | 0.507 |

Positive entropy values above are the unfavorable cost **−TΔS**. HydraRank's dimensionless entropies were converted with −RTΔS/R at 303.15 K. SSTMap's reported 300 K TΔS values were sign-inverted and multiplied by 303.15/300; this is an exact temperature normalization for SSTMap's linear T·S conversion and does not require rerunning its trajectory analysis. The bootstrap intervals in `metrics.json` describe across-site variation only and are **not trajectory-sampling confidence intervals**.

## What can and cannot be compared

- Occupancy, site positions, and translational/orientational/total entropy are direct matched-site comparisons after unit/sign normalization.
- Solute hydrogen bonds are informative but not estimator-equivalent: HydraRank uses a >130° donor–H–acceptor threshold, while SSTMap uses the stricter equivalent of >150°.
- SSTMap provides water–solute, water–water, and total interaction energies; HydraRank has no energy analogue.
- HydraRank provides residence time, enclosure, category, and a displacement ranking; SSTMap has no direct analogues. Correlations between the HydraRank score and SSTMap energies are therefore exploratory, not a validation of the rank.
- SSTMap's `f_enc` is undefined in this HSA implementation, so it cannot validate HydraRank enclosure.

## Ranking comparison

SSTMap does not provide a displacement ranking, so a direct native-rank comparison is impossible. As a diagnostic, applying HydraRank's score form to SSTMap outputs—SSTMap −TΔS minus 1 kcal/mol per SSTMap solute H-bond—gives Pearson r=0.660 and Spearman ρ=0.682 against the HydraRank score. The top-3, top-5, and top-10 overlaps are 2/3, 3/5, and 8/10. This is meaningful moderate agreement, but the proxy is not an independent SSTMap ranking and its H-bond definition differs.

## Alignment sensitivity

The locally aligned definitive SSTMap calculation gives 21 sites; the earlier whole-protein-aligned control gives 19. There are 17 matched sites within 1 Å. Their center shift is 0.476 Å on average (maximum 0.842 Å); total entropy has Pearson r=0.999 between alignments. This confirms that local-pocket alignment is not a cosmetic choice and that only the local-fit run should be used for the primary controlled comparison.

## Interpretation and next validation steps

For 3RLP, the tools agree exceptionally well on where the principal waters are, but quantitative thermodynamic agreement must be judged from the statistics above rather than inferred from spatial overlap. Any systematic entropy offset is scientifically plausible because the tools do not have identical estimators and implementation details.

The next defensible validation step is to repeat the controlled analysis for 3RLQ and 3RLR, then measure convergence by trajectory blocks or independent replicas. Experimental displacement data or ligand-series affinity changes would be needed to test whether HydraRank's product-level ranking is predictive; SSTMap alone is a computational comparator, not ground truth.

## Reproducibility

Both methods used 25,001 frames at 2 ps spacing, the same 61-atom local-pocket fit, a 5 Å hydration region, 1 Å sites, 0.0333 Å⁻³ bulk density, and a 2× density threshold. Entropy costs are reported at the simulation temperature of 303.15 K. Exact input checksums are recorded in `metrics.json`. Recreate every table, metric, and figure with:

```bash
uv run python comparison/3rlp/compare.py
```
