# Reproducible performance benchmark

For the controlled HydraRank–SSTMap HSA/GIST pilot, resource monitoring, raw
measurements, and implementation limitations, see
[`comparison/performance/`](../comparison/performance/README.md). The commands
below are lightweight HydraRank-only smoke timings, not a cross-tool comparison.

The bundled 200 ps 3RLP slice is a software smoke-test and performance fixture, not a
scientifically converged dataset. Run benchmarks from the repository root on an otherwise
idle machine and record the package version, Python version, hardware, and command output.

```bash
uv sync --all-extras
uv run hydrarank --version
uv run python --version

/usr/bin/time -p uv run hydrarank preprocess \
  -c examples/hsp90_demo.yaml \
  -o output/benchmark-observations.npz \
  2> output/preprocess-time.txt

/usr/bin/time -p uv run hydrarank sites \
  -c examples/hsp90_demo.yaml \
  --observations output/benchmark-observations.npz \
  --json > output/benchmark-sites.json \
  2> output/sites-time.txt

/usr/bin/time -p uv run hydrarank rank \
  -c examples/hsp90_demo.yaml \
  --observations output/benchmark-observations.npz \
  --json > output/benchmark-ranking.json \
  2> output/ranking-time.txt
```

The compact fixture is useful for reproducible smoke timings, not meaningful scaling curves.
Use a full trajectory and repeat with increasing `--stop` values for scaling measurements,
holding `--start` and `--step` fixed. Do not compare scientific values across differently
sampled trajectories; these commands measure runtime and memory behavior only.
