# HSP90 apo–holo validation

This directory contains the reproducible comparison of the apo HSP90 hydration map with
the 3RLP, 3RLQ, and 3RLR holo trajectories and crystallographic structures. It addresses
the water-displacement sequence discussed by Kung et al. (2011) without treating three
complexes as a statistically meaningful affinity model.

Run from the repository root after the four HydraRank analyses have been generated:

```bash
.venv/bin/python validation/hsp90/run_validation.py
```

The script reads the compressed observation caches, not the complete trajectories. It
builds a common 3RLP pocket frame, maps crystallographic waters, measures fixed-site
occupancy, checks 10 ns blocks, evaluates one-factor-at-a-time sensitivity, and writes a
hash manifest for the inputs that determine the report.

The main interpretation is in `summary.md`. CSV files retain the site-level data used to
make the figures. Important limitations are part of the result: only one trajectory is
available per state, the apo system was initialized with crystallographic waters, and a
crystal-water position is structural evidence rather than a thermodynamic label.
