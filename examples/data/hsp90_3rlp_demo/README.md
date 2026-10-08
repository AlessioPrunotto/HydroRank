# HSP90 3RLP software fixture

This directory contains the original topology and the first 101 contiguous frames
(0–200 ps at 2 ps spacing) from the 3RLP production trajectory. It exists only to exercise
HydraRank's real-file workflow from a fresh checkout. The slice is far too short for
scientific interpretation of residence times, entropy, or displacement rankings.

The trajectory slice was created without changing atom order or coordinates:

```bash
printf 'System\n' | gmx trjconv \
    -s sample_traj/3rlp/step5.tpr \
    -f sample_traj/3rlp/step5.xtc \
    -o examples/data/hsp90_3rlp_demo/step5_0_200ps.xtc \
    -b 0 -e 200
```

Scientific validation uses the complete 50 ns HSP90 trajectories described under
`validation/hsp90/` and in `paper/manuscript.md`.
