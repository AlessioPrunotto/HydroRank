# Python API

This doc shows how to use hydrarank in python (rather than the more common command line version that is described in the README).
All coordinates use ångström, trajectory times use picoseconds, and energies use kcal/mol.

## Loading and configuration

- `hydrarank.PreprocessConfig`: validated configuration with YAML round-trip support.
- `hydrarank.io.load_universe(config)`: load topology and trajectories.
- `hydrarank.io.describe_system(universe, config)`: inspect selections and sampling.
- `hydrarank.qc.run_preprocess_qc(universe, config)`: run PBC/alignment diagnostics.

## Analysis pipeline

The supported high-level API runs the complete workflow and writes a reproducible bundle:

```python
from hydrarank import PreprocessConfig, run_analysis

config = PreprocessConfig(
    topology="system.tpr",
    trajectory=["trajectory.xtc"],
    ligand_selection="resname LIG",
)
result = run_analysis(config, output_dir="hydrarank-results")
print(result.format(top=10))
```

For custom pipelines, the individual stages remain public:

```python
from hydrarank.analysis import analyse_sites
from hydrarank.clustering import cluster_hydration_sites
from hydrarank.config import PreprocessConfig
from hydrarank.io import load_universe
from hydrarank.preprocess import run_preprocess
from hydrarank.ranking import rank_sites

config = PreprocessConfig.from_yaml("analysis.yaml")
observations = run_preprocess(load_universe(config), config)
sites = cluster_hydration_sites(
    observations,
    radius=config.site_radius,
    density_factor=config.density_factor,
)
analysis = analyse_sites(
    observations,
    sites,
    temperature=config.temperature,
    max_gap=config.max_gap,
)
ranking = rank_sites(analysis, hbond_penalty=config.hbond_penalty)
```

`WaterObservations` is the persistent boundary between trajectory processing and analysis.
Use `save(path)` and `WaterObservations.load(path)` to avoid rereading a trajectory.
The high-level workflow validates its cache against the input paths, file sizes and
modification times, plus every preprocessing parameter.

For an apo trajectory, provide a bound reference whose protein atom names and residue
numbering match the apo topology:

```python
config = PreprocessConfig(
    topology="apo.tpr",
    trajectory=["apo.xtc"],
    reference_structure="bound.tpr",
    reference_coordinates="bound.xtc",
    reference_ligand_selection="resname LIG",
)
result = run_analysis(config, output_dir="hydrarank-apo-results")
```

The reference ligand defines the spatial region only and is excluded from apo interaction
counts. `reference_coordinates` is optional for coordinate-bearing formats such as PDB or
GRO, but required when `reference_structure` is a TPR.

## Results and exports

- `SiteAnalysis.rows()` and `SiteRanking.rows()` return JSON/CSV-friendly dictionaries.
- `format_analysis()` and `format_ranking()` create terminal tables.
- `write_csv(rows, path)` writes either table to CSV.
- `write_site_coordinates(rows, path)` writes site centres to PDB or mmCIF.
- `write_analysis_plots(analysis, directory)` writes occupancy, residence, and convergence plots.
- `write_ranking_plot(ranking, directory)` writes a spatial score map.

The plotting API requires the `plots` optional dependency. The ranking score is a heuristic
prioritisation signal, not a computed displacement free energy.

## Progress callbacks

`run_preprocess` accepts an optional `progress(current_frame, total_frames)` callback. Library
calls remain silent by default; the CLI supplies a throttled stderr reporter unless
`--no-progress` is used.
