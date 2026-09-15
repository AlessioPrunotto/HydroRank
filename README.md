# water-entropy

Light-weight hydration-site analysis for ligand design: given an MD trajectory of a
protein–ligand complex, identify the water molecules around the ligand that would be
**entropically favourable to displace**.

Existing tools for this are either commercial (WaterMap), unmaintained and pinned to
old Python (SSTMap, WaterKit), or expensive (GIST). This package aims to be a small,
installable, testable alternative built on MDAnalysis, NumPy and SciPy.

> **Status: early development.** The pipeline runs end to end: preprocessing, hydration
> sites, occupancy, residence times, entropy proxies, solute hydrogen bonds and a
> heuristic ranking of displaceable waters. The ranking is intended for prioritisation,
> not as a rigorously validated displacement free energy.

## Installation

```bash
uv venv --python 3.14
uv sync
```

## Quick start

> **Demonstration data only.** The bundled kinase trajectory contains 301 frames spaced
> 100 ps apart. This is useful for exercising the complete software workflow, but it is
> too short and too coarsely sampled for reliable residence times, entropy estimates, or
> displacement rankings. The numerical results below demonstrate the output format and
> must not be interpreted as scientific conclusions about this ligand or binding site.

Inspect a system and check that the selections do what you expect:

```bash
uv run water-entropy info \
    -s sample_traj/gromacs_traj/step3_input.psf \
    -f sample_traj/gromacs_traj/step5.xtc \
    -l "resname 547"
```

```
System
  atoms              : 77907
  frames             : 301 (dt = 100 ps)
  box                : 86.00 x 86.00 x 86.00 A, angles 90.0/90.0/90.0 (orthorhombic)
  bonds in topology  : yes
Water
  molecules          : 18337 (resname OPC)
  model              : 4-site (TIP4P/OPC-like)
Ligand
  selection          : resname 547
  atoms              : 56 (30 heavy)
Warnings
  - effective frame spacing is 100 ps, which is longer than typical water residence
    times (1-100 ps); persistence and entropy estimates will be unreliable
```

Then check that the periodic-boundary treatment and the binding-site fit are sound:

```bash
uv run water-entropy check -c examples/kinase.yaml
```

```
Preprocessing QC
  frames analysed    : 301
  fit atoms          : 55 (reference frame 0)
  fit RMSD           : mean 1.04 A, max 1.35 A, last 1.35 A
  longest bond       : 1.93 A (solute whole)
  waters within 5 A: mean 32.1 (min 18, max 40)
```

Then cluster the first-shell waters into hydration sites:

```bash
uv run water-entropy sites -c examples/kinase.yaml
```

```
Hydration sites (radius 1 A, 301 frames at 100 ps, T = 300 K)
  21 sites, 3124 of 7585 water observations assigned

 site       x       y       z   wat  occup  res/ps  encl    hb    S_tr    S_or    -TdS
    1   46.64   39.64   45.43   248   0.82     420  21.9  1.87   -2.63   -2.67    3.16
    2   41.15   46.12   38.44   231   0.77     608  14.4  0.03   -2.81   -3.33    3.66
    3   43.39   46.12   40.21   215   0.71    7167  19.6  1.13   -3.10   -4.53    4.55
  ...

S_tr, S_or: excess translational / orientational entropy per water relative to
bulk, in units of the gas constant (negative = more ordered than bulk).
-TdS: free-energy cost of that ordering at the given temperature, in kcal/mol.
encl: mean nearby solute heavy atoms; hb: mean solute hydrogen bonds per water.
```

Rank the sites by entropy released minus a configurable hydrogen-bond penalty:

```bash
uv run water-entropy rank -c examples/kinase.yaml \
    --observations output/observations.npz --top 5
```

```
rank site  occup  encl    hb   -TdS  score  category
   1    2   0.77  14.4  0.03   3.66   3.63  displaceable
   2    3   0.71  19.6  1.13   4.55   3.41  replace-hbonds
   3    9   0.50  17.1  0.68   3.34   2.66  displaceable
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
uv run water-entropy preprocess -c examples/kinase.yaml -o output/observations.npz
uv run water-entropy sites -c examples/kinase.yaml --observations output/observations.npz
uv run water-entropy rank -c examples/kinase.yaml --observations output/observations.npz
```

The same in Python:

```python
from water_entropy import PreprocessConfig
from water_entropy.analysis import analyse_sites, format_analysis
from water_entropy.clustering import cluster_hydration_sites
from water_entropy.io import describe_system, load_universe
from water_entropy.preprocess import run_preprocess
from water_entropy.ranking import format_ranking, rank_sites

config = PreprocessConfig(
    topology="step3_input.psf",
    trajectory=["step5.xtc"],
    ligand_selection="resname 547",
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

Configuration can also live in YAML — see [examples/kinase.yaml](examples/kinase.yaml).

## Package layout

| module | role |
| --- | --- |
| [config.py](src/water_entropy/config.py) | `PreprocessConfig`, YAML round-trip, validation |
| [io.py](src/water_entropy/io.py) | universe loading, frame ranges, `SystemReport` |
| [selections.py](src/water_entropy/selections.py) | ligand / water / pocket selection, water-model detection |
| [pbc.py](src/water_entropy/pbc.py) | unwrap → centre → wrap-by-residue transformation stack |
| [alignment.py](src/water_entropy/alignment.py) | binding-site reference frame and rigid-body fitting |
| [preprocess.py](src/water_entropy/preprocess.py) | orchestration: raw trajectory → aligned water observations |
| [data.py](src/water_entropy/data.py) | `WaterObservations`, the contract between stages, with `.npz` I/O |
| [clustering.py](src/water_entropy/clustering.py) | density-peak clustering into hydration sites |
| [entropy.py](src/water_entropy/entropy.py) | nearest-neighbour translational and orientational entropy estimators |
| [analysis.py](src/water_entropy/analysis.py) | occupancy, residence times and the final site table |
| [hbonds.py](src/water_entropy/hbonds.py) | water-solute hydrogen bonds and enclosure proxy |
| [ranking.py](src/water_entropy/ranking.py) | heuristic displacement scoring and categories |
| [qc.py](src/water_entropy/qc.py) | per-frame diagnostics of the above |
| [cli.py](src/water_entropy/cli.py) | `water-entropy info` / `check` / `preprocess` / `sites` / `rank` |

The stages communicate through one array table, `WaterObservations`: one row per
(frame, water) pair, holding the oxygen and hydrogen positions in the aligned frame
plus the frame index, display residue id, and globally unique topology residue index,
so later stages never touch the trajectory again. Version-2 caches remain readable;
new caches use format 3 and preserve unique water identity across repeated residue
numbers in different segments.

Options passed explicitly on the command line override values loaded from a YAML
configuration file.

## Design notes

- **Water is matched by an explicit resname list**, not by the built-in `water`
  keyword: MDAnalysis and MDTraj do not know about `OPC`, which silently yields an
  empty selection (the same bug that makes SSTMap fail on OPC systems).
- **Atom roles are assigned by mass, not by name.** Water oxygens are called `OW`,
  `OH2` or `O` depending on the force field, and 4-/5-site models carry a massless
  virtual site (`MW`, `EPW`) that must be excluded from all geometry.
- **The tool refuses to guess.** Ligand auto-detection raises rather than picking one
  of several candidates, and empty selections raise with the list of residue names
  actually present.
- **Sampling adequacy is reported, not assumed.** Frame spacing and frame count are
  checked against water residence timescales up front.
- **The PBC stack has one correct order**: unwrap the solute, centre the ligand in the
  box, then wrap everything else *by residue*. Wrapping first splits molecules; wrapping
  by atom splits waters. After this, plain Euclidean distances are valid in the pocket,
  and the rest of the pipeline can ignore periodicity. Cutoffs are checked against half
  the shortest box vector so the minimum-image convention is never violated.
- **Frames are fitted on the binding-site backbone**, not on the whole protein (domain
  motions smear the pocket) and not on the ligand alone (unstable for small or symmetric
  ligands, and it would make the protein move instead). The per-frame fit RMSD is
  reported because a poor fit inflates the apparent positional spread of the waters,
  which downstream looks exactly like disorder and produces false "easy to displace"
  hits.
- **The hydration shell is defined against the reference ligand pose**, in the aligned
  frame, not against the moving ligand. A region that follows the ligand's wobble would
  make site occupancies depend on ligand motion rather than on water behaviour.
- **Sites come from density-peak clustering**, the scheme WaterMap and SSTMap use:
  repeatedly take the position with the most neighbours within 1 Å, call it a site,
  remove the waters it claims, and stop when no remaining peak is denser than bulk water
  (0.0329 molecules Å⁻³). It needs only a KD-tree, is deterministic, and yields sites of
  a fixed physical radius — unlike k-means, which needs the number of sites up front, or
  DBSCAN, whose clusters can grow into elongated blobs spanning several real sites.
- **Entropies are measured against bulk water, not against nothing.** Both estimators are
  nearest-neighbour (Kozachenko-Leonenko) estimates of the excess entropy per water: the
  translational term against a uniform fluid at 0.0329 molecules Å⁻³, the orientational
  term against uniformly random orientations. Zero therefore means "already bulk-like,
  nothing to gain", which is exactly the question being asked. This is the first-order
  inhomogeneous-solvation-theory approximation of WaterMap and SSTMap: water-water
  correlations are ignored, which is what makes it cheap.
- **The two hydrogens are treated as indistinguishable.** A half turn about the dipole
  axis maps a water onto itself, so the orientational estimator searches the symmetry
  orbit and folds the reference measure by the same factor. Skipping this makes every
  site look ~0.69 gas-constant units more ordered than it is.
- **Estimators return `nan` rather than a plausible-looking number** when a site holds
  fewer than ten observations.
- **Ranking combines ordering and solute interactions.** The displacement score subtracts
  a configurable penalty per mean water-solute hydrogen bond from `-TdS`. Categories carry
  the main interpretation; the heuristic score is only a tie-breaker.

## Development

```bash
uv run pytest
uv run ruff check src tests
uv run ruff format src tests
```

## Roadmap

1. ~~Scaffold, IO and selection layer~~
2. ~~PBC treatment (unwrap / centre / wrap by residue) and alignment to the binding site~~
3. ~~Clustering of water oxygen positions into hydration sites~~
4. ~~Occupancy, persistence and entropy proxies~~
5. ~~Heuristic ranking using water-solute hydrogen bonds and enclosure~~
6. Scientific validation, uncertainty estimates and benchmark comparisons

## License

MIT
