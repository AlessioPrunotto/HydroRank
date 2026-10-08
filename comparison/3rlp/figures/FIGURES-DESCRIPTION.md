### Description of each figure present in this folder:

 - [site_overlay.png](/Users/aprunotto/Documents/alessio-code/water-entropy/comparison/3rlp/figures/site_overlay.png)
Shows the 3D positions of hydration-site centres: blue circles for HydraRank and orange triangles for SSTMap, with lines connecting matched sites. All 19 HydraRank sites have SSTMap counterparts; SSTMap detects two additional sites. The matched centres are only 0.103 Å apart on average, showing close spatial agreement.

 - [occupancy_comparison.png](/Users/aprunotto/Documents/alessio-code/water-entropy/comparison/3rlp/figures/occupancy_comparison.png)
Compares how frequently each matched site contains water, with SSTMap on the horizontal axis and HydraRank on the vertical axis. Most sites lie near the diagonal, although HydraRank reports substantially lower occupancy for two sites. Overall correlation is strong (Pearson r = 0.948). Here, only the matched sites are reported (the extra 2 sites that sstmap identified are not on this plot)

 - [entropy_comparison.png](/Users/aprunotto/Documents/alessio-code/water-entropy/comparison/3rlp/figures/entropy_comparison.png)
Compares translational, orientational, and total entropy costs in three panels. Values are expressed as −TΔS in kcal/mol, normalized to 303.15 K: larger positive values mean a greater cost of ordering water relative to bulk. Agreement is strong, particularly for the total (r = 0.994), though individual values differ. Again, data are for the 19 matched sites

 - [hbond_comparison.png](/Users/aprunotto/Documents/alessio-code/water-entropy/comparison/3rlp/figures/hbond_comparison.png)
Compares the average number of hydrogen bonds between site waters and the solute. HydraRank generally counts more. The methods use different angular thresholds—HydraRank >130°, SSTMap >150°—so this figure compares their native definitions rather than identical measurements.

 - [ranking_context.png](/Users/aprunotto/Documents/alessio-code/water-entropy/comparison/3rlp/figures/ranking_context.png)
Places the HydraRank displacement score against three SSTMap quantities: a derived score proxy, solute–water interaction energy, and total interaction energy. The proxy uses the same score form, entropy cost minus a hydrogen-bond penalty, and shows moderate agreement (r = 0.660). The energy panels explore whether the score tracks interaction energies; SSTMap does not supply a native displacement ranking. For the central and right panel, it makes sense to have values that refer to sstmap only, because it depends on how the hydration sites were computed and identified.

 - [alignment_sensitivity.png](/Users/aprunotto/Documents/alessio-code/water-entropy/comparison/3rlp/figures/alignment_sensitivity.png)
Shows how SSTMap results change when frames are aligned to the local pocket versus the whole protein. The first two panels compare occupancy and total entropy cost; the third shows the distribution of site-centre shifts. Across 17 matched sites, local alignment generally increases occupancy and entropy cost, while centres shift by 0.476 Å on average, demonstrating that alignment affects the results.