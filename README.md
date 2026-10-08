# HydraRank

HydraRank is a light-weight hydration-site analysis for supporting ligand design: given an MD trajectory of a
protein–ligand complex, it identifies the water molecules around the ligand that would be
**entropically favourable to displace**.

Existing tools for this are either commercial (WaterMap), unmaintained and pinned to
old Python (SSTMap, WaterKit), or computationally expensive (GIST). This package aims to be a small,
installable, testable alternative built on MDAnalysis, NumPy and SciPy.

> **Status: research software with retrospective validation.** The complete analysis
> workflow is implemented and covered by automated tests. Retrospective HSP90 validation
> reproduces the W3 and W249 displacement patterns, while exposing limitations for W286
> and apo sampling; a matched-trajectory SSTMap comparison supports agreement in site
> positions and entropy ordering. Ranking remains a prioritisation heuristic, not a
> displacement free energy or binding-affinity prediction. Independent-replica
> reproducibility and prospective predictive performance remain to be established.
> See the [HSP90 validation](validation/hsp90/summary.md) and
> [SSTMap comparison](comparison/3rlp/report.md) for evidence and limitations.

## Installation

From a checkout, install HydraRank and its dependencies into the project environment with:

```bash
uv sync
```

`uv sync` creates `.venv` when it does not exist and safely reuses it when it does.
HydraRank supports Python 3.11 or newer. To explicitly select a Python version for a new
environment, use `uv sync --python 3.14`. Alternatively, use `pip install .`. Install the
optional plotting dependency with `uv sync --extra plots` or `pip install '.[plots]'`.

Run the command through the project environment with `uv run hydrarank`. If you prefer to
invoke `hydrarank` directly, either activate the environment or install it as a uv tool:

```bash
source .venv/bin/activate
# or, without activating an environment:
uv tool install .
```

## One-command analysis

For a prepared protein–ligand MD trajectory, the complete workflow is one command:

```bash
uv run hydrarank analyse \
    --topology system.tpr \
    --trajectory trajectory.xtc \
    --ligand-selection "resname LIG" \
    --output-dir hydrarank-results
```

`analyze` is accepted as an alias for `analyse`. HydraRank inspects the system, applies
PBC treatment and binding-site alignment, collects QC while preprocessing, clusters and
analyses hydration sites, ranks them, prints a human-readable report, and writes:

- `report.txt`: the same readable system/QC/ranking report;
- `results.json`: complete machine-readable results and provenance;
- `ranking.csv`: the ranked hydration-site table;
- `sites.pdb`: site centres for PyMOL, ChimeraX, or VMD;
- `observations.npz`: the reusable preprocessed observations;
- `config.yaml`: the effective configuration used for the run.

On a repeat run, preprocessing is skipped only when the cached input file signatures and
all preprocessing parameters still match. Use `--force` to rebuild it. Add `--plots` to
write diagnostic figures when the optional plotting dependency is installed.

The topology must provide atom identities and masses, molecular bonds, water hydrogens,
and periodic box information; the trajectory supplies coordinates and time. A TPR/XTC
pair is a convenient GROMACS choice, while PSF/DCD and other combinations supported by
MDAnalysis also work when they carry the same required information.

### Apo trajectories

An apo trajectory has no ligand from which to define a reproducible hydration region.
Supply a protein-ligand reference structure and select its ligand instead:

```bash
uv run hydrarank analyse \
    --topology apo.tpr \
    --trajectory apo.xtc \
    --reference-structure bound.tpr \
    --reference-coordinates bound.xtc \
    --reference-ligand-selection "resname LIG" \
    --output-dir hydrarank-apo-results
```

HydraRank maps the reference pocket atoms onto the apo protein, places the reference
ligand pose in the apo frame, and uses it only to define the pocket and hydration volume.
The ligand is not added to the simulation and does not contribute hydrogen bonds or
enclosure. Protein residue numbering and atom names must match uniquely. A standalone
PDB, GRO, or mmCIF can be used without `--reference-coordinates`; a TPR reference requires
explicit coordinates because its embedded coordinates are not a reliable spatial reference.

## Demonstration and individual stages

> **Demonstration data only.** The bundled 3RLP HSP90 trajectory slice contains 101
> contiguous frames spaced 2 ps apart (200 ps total). It exercises the complete software
> workflow, but is far too short for reliable residence times, entropy estimates, or
> displacement rankings. The numerical results below demonstrate the output format and
> must not be interpreted as scientific conclusions about this ligand or binding site.

The bundled system can run through the new workflow directly:

```bash
uv run hydrarank analyse -c examples/hsp90_demo.yaml --output-dir output/demo
```

The lower-level commands remain available for debugging and custom pipelines. Inspect a
system and check that the selections do what you expect:

```bash
uv run hydrarank info \
    -s examples/data/hsp90_3rlp_demo/step5.tpr \
    -f examples/data/hsp90_3rlp_demo/step5_0_200ps.xtc \
    -l "resname 3RP"
```

```
System
  atoms              : 62822
  frames             : 101 (dt = 2 ps)
  box                : 80.00 x 80.00 x 80.00 A, angles 90.0/90.0/90.0 (orthorhombic)
  bonds in topology  : yes
Water
  molecules          : 14799 (resname OPC)
  model              : 4-site (TIP4P/OPC-like)
Ligand
  selection          : resname 3RP
  atoms              : 29 (18 heavy)
Warnings
  - only 101 frames selected; per-site statistics will be noisy (a few thousand frames is a reasonable target)
```

Then check that the periodic-boundary treatment and the binding-site fit are sound:

```bash
uv run hydrarank check -c examples/hsp90_demo.yaml
```

```
Preprocessing QC
  frames analysed    : 101
  fit atoms          : 61 (reference frame 0)
  fit RMSD           : mean 0.53 A, max 0.66 A, last 0.52 A
  longest bond       : 1.93 A (solute whole)
  waters within 5 A: mean 18.4 (min 14, max 21)
```

Then cluster the first-shell waters into hydration sites:

```bash
uv run hydrarank sites -c examples/hsp90_demo.yaml
```

```
Hydration sites (radius 1 A, 101 frames at 2 ps, T = 303.15 K)
  21 sites, 1342 of 1860 water observations assigned

 site       x       y       z   wat  occup  res/ps  encl    hb    S_tr    S_or    -TdS
    1   38.05   43.32   38.51   101   1.00     202  27.0  1.82   -3.52   -4.43    4.79
    2   39.23   45.81   39.32   101   1.00     202  21.6  2.02   -3.53   -4.19    4.65
    3   45.11   43.29   42.02   101   1.00     202  31.7  4.00   -4.40   -5.35    5.87
  ...

S_tr, S_or: excess translational / orientational entropy per water relative to
bulk, in units of the gas constant (negative = more ordered than bulk).
-TdS: free-energy cost of that ordering at the given temperature, in kcal/mol.
encl: mean nearby solute heavy atoms; hb: mean solute hydrogen bonds per water.
```

Rank the sites by entropy released minus a configurable hydrogen-bond penalty:

```bash
uv run hydrarank rank -c examples/hsp90_demo.yaml \
    --observations output/observations.npz --top 5
```

```
rank site  occup  encl    hb   -TdS  score  category
   1   19   0.36  19.3  1.00   4.32   3.32  replace-hbonds
   2    6   0.91  24.3  0.99   4.09   3.10  displaceable
   3    1   1.00  27.0  1.82   4.79   2.97  replace-hbonds
```

`displaceable` sites are ordered and weakly hydrogen bonded to the solute.
`replace-hbonds` sites are ordered but should only be targeted if the ligand can replace
their interactions. These labels explain how to interpret a sufficiently sampled
analysis; on this demonstration trajectory they are illustrative only. The numeric score
is a heuristic tie-breaker, not a free energy.

For scientific use, analyse a trajectory saved frequently enough to resolve water motion
and containing enough effectively independent observations for each site. Inspect the
warnings from `info`, the alignment diagnostics from `check`, and convergence across
trajectory blocks before interpreting residence, entropy, or ranking values. Required
sampling depends on the system, so the package deliberately does not declare a universal
minimum trajectory length.

Extracting the water observations is the expensive part, so it can be cached:

```bash
uv run hydrarank preprocess -c examples/hsp90_demo.yaml -o output/observations.npz
uv run hydrarank sites -c examples/hsp90_demo.yaml --observations output/observations.npz
uv run hydrarank rank -c examples/hsp90_demo.yaml --observations output/observations.npz
```

Preprocessing reports progress on stderr. Pass `--no-progress` for quiet batch jobs.

## Exporting results

Both `sites` and `rank` can write CSV tables, molecular-viewer coordinates, and diagnostic
plots alongside their terminal or JSON output:

```bash
uv run hydrarank rank -c examples/hsp90_demo.yaml \
    --observations output/observations.npz \
    --csv output/ranking.csv \
    --site-coordinates output/sites.pdb \
    --plot-dir output/plots
```

Use a `.pdb`, `.cif`, or `.mmcif` suffix for `--site-coordinates`. In coordinate files,
occupancy is stored as occupancy and `-TΔS` as the B-factor, making the sites directly
viewable with the protein in PyMOL, ChimeraX, or VMD. Plot output includes site occupancy,
residence-time distributions, entropy convergence, and—for `rank`—a spatial score map.

The same in Python:

```python
from hydrarank import PreprocessConfig
from hydrarank.analysis import analyse_sites, format_analysis
from hydrarank.clustering import cluster_hydration_sites
from hydrarank.io import describe_system, load_universe
from hydrarank.preprocess import run_preprocess
from hydrarank.ranking import format_ranking, rank_sites

config = PreprocessConfig(
    topology="examples/data/hsp90_3rlp_demo/step5.tpr",
    trajectory=["examples/data/hsp90_3rlp_demo/step5_0_200ps.xtc"],
    ligand_selection="resname 3RP",
)
universe = load_universe(config)
report, water_topology = describe_system(universe, config)
print(report)

observations = run_preprocess(universe, config)
sites = cluster_hydration_sites(observations)
analysis = analyse_sites(observations, sites)
print(format_analysis(analysis))
print(format_ranking(rank_sites(analysis)))
```

Configuration can also live in YAML — see [examples/hsp90_demo.yaml](examples/hsp90_demo.yaml).

## Package layout

| module | role |
| --- | --- |
| [config.py](src/hydrarank/config.py) | `PreprocessConfig`, YAML round-trip, validation |
| [io.py](src/hydrarank/io.py) | universe loading, frame ranges, `SystemReport` |
| [selections.py](src/hydrarank/selections.py) | ligand / water / pocket selection, water-model detection |
| [pbc.py](src/hydrarank/pbc.py) | unwrap → centre → wrap-by-residue transformation stack |
| [alignment.py](src/hydrarank/alignment.py) | binding-site reference frame and rigid-body fitting |
| [preprocess.py](src/hydrarank/preprocess.py) | orchestration: raw trajectory → aligned water observations |
| [data.py](src/hydrarank/data.py) | `WaterObservations`, the contract between stages, with `.npz` I/O |
| [clustering.py](src/hydrarank/clustering.py) | density-peak clustering into hydration sites |
| [entropy.py](src/hydrarank/entropy.py) | nearest-neighbour translational and orientational entropy estimators |
| [analysis.py](src/hydrarank/analysis.py) | occupancy, residence times and the final site table |
| [hbonds.py](src/hydrarank/hbonds.py) | water-solute hydrogen bonds and enclosure proxy |
| [ranking.py](src/hydrarank/ranking.py) | heuristic displacement scoring and categories |
| [export.py](src/hydrarank/export.py) | CSV and PDB/mmCIF artifact export |
| [plotting.py](src/hydrarank/plotting.py) | optional occupancy, residence, convergence and ranking plots |
| [qc.py](src/hydrarank/qc.py) | per-frame diagnostics of the above |
| [workflow.py](src/hydrarank/workflow.py) | unified analysis, cache validation, provenance, and result bundle |
| [cli.py](src/hydrarank/cli.py) | `hydrarank analyse` and the lower-level stage commands |

The stages communicate through one array table, `WaterObservations`: one row per
(frame, water) pair, holding the oxygen and hydrogen positions in the aligned frame
plus the frame index, display residue id, and globally unique topology residue index,
so later stages never touch the trajectory again. Version-2 caches remain readable;
new caches use format 3 and preserve unique water identity across repeated residue
numbers in different segments.

Options passed explicitly on the command line override values loaded from a YAML
configuration file.

## Method and limitations

- **Inputs and selections:** complete analysis requires masses, bonds, explicit water
  hydrogens, and periodic-box information. Water residue names are configurable, and
  ambiguous ligand selections stop the analysis.
- **PBC and local alignment:** HydraRank makes the solute whole, centres the ligand or
  apo pocket anchor, and wraps other molecules by residue before fitting the local
  protein pocket. Longest-bond and fit-RMSD diagnostics help identify preprocessing
  problems; passing them does not establish scientific reliability.
- **Fixed analysis region:** water observations are collected around a fixed reference
  ligand pose in the aligned frame. For apo analysis, that pose defines the region
  without adding ligand interactions to the simulation.
- **Site definition:** density-peak clustering uses a configurable radius (default
  1.0 Å) and observation-count threshold (default twice bulk density). The resulting
  number and identity of sites depend on these choices.
- **Entropy and ranking:** first-order excess-entropy proxies describe water ordering
  relative to bulk. The score subtracts a heuristic water–solute hydrogen-bond penalty;
  enclosure is reported separately. Interaction energies and higher-order solvent
  correlations are omitted, so scores and categories do not establish displacement
  free energies or affinity changes.
- **Sampling:** frame-count and spacing warnings flag potential problems, but do not
  establish adequate sampling. Inspect block convergence and, where available,
  independent replicas before interpreting occupancy, residence, entropy, or ranking.
  Ten valid observations permit an entropy estimate; they do not guarantee convergence.

See [Methodology and interpretation](docs/methodology.md) for the design rationale,
defaults, and diagnostic limits, and the [HSP90 validation](validation/hsp90/summary.md)
for the current scientific evidence.

## Development

```bash
uv run pytest
uv run ruff check src tests
uv run ruff format src tests
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution requirements,
[docs/api.md](docs/api.md) for the Python API, and
[docs/benchmarking.md](docs/benchmarking.md) for reproducible performance commands.
Release changes are recorded in [CHANGELOG.md](CHANGELOG.md); citation metadata is provided
in [CITATION.cff](CITATION.cff). A journal-neutral draft software manuscript and its BibTeX
library are available in [paper/manuscript.md](paper/manuscript.md) and
[paper/references.bib](paper/references.bib). The reproducible apo–holo HSP90 validation,
including crystallographic-water recovery, block convergence, sensitivity tests, and the
paper-target analysis for W3, W249, and W286, is in
[validation/hsp90/summary.md](validation/hsp90/summary.md).

## Roadmap

1. ~~Scaffold, IO and selection layer~~
2. ~~PBC treatment (unwrap / centre / wrap by residue) and alignment to the binding site~~
3. ~~Clustering of water oxygen positions into hydration sites~~
4. ~~Occupancy, persistence and entropy proxies~~
5. ~~Heuristic ranking using water-solute hydrogen bonds and enclosure~~
6. ~~Retrospective HSP90 validation and matched-trajectory SSTMap comparison~~
7. Independent-replica and prospective validation on a larger ligand series

## License

[MIT](LICENSE)
