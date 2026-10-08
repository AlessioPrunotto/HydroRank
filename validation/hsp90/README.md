# HSP90 apo–holo validation

This directory contains the reproducible comparison of the apo HSP90 hydration map with
the 3RLP, 3RLQ, and 3RLR holo trajectories and crystallographic structures. It addresses
the water-displacement sequence discussed by Kung et al. (2011) without treating three
complexes as a statistically meaningful affinity model.

The four `hsp90_*.yaml` configurations in this directory define the full 50 ns
analyses. They require the local production topology and trajectory files under
`sample_traj/3rlp/`, `sample_traj/3rlq/`, `sample_traj/3rlr/`, and
`sample_traj/apo_form/`. These inputs are ignored by Git and must be supplied
separately; they are not included in a fresh checkout.

Generate the four HydraRank analyses from the repository root:

```bash
uv run hydrarank analyse -c validation/hsp90/hsp90_apo.yaml
uv run hydrarank analyse -c validation/hsp90/hsp90_3rlp.yaml
uv run hydrarank analyse -c validation/hsp90/hsp90_3rlq.yaml
uv run hydrarank analyse -c validation/hsp90/hsp90_3rlr.yaml
```

The configurations write to `output/hsp90/<system>/`. Then run the validation:

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

For a small software demonstration with included inputs, use
[`examples/hsp90_demo.yaml`](../../examples/hsp90_demo.yaml). Its 200 ps trajectory
is too short for scientific validation.
