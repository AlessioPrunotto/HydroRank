# HSP90 validation figures: captions and reading guide

These four figures compare hydration in one apo simulation (HSP90 without a bound
ligand) and three holo simulations (HSP90 bound to the ligands from 3RLP, 3RLQ,
and 3RLR). Each trajectory is 50 ns long, with 25,001 frames at 2 ps spacing.
Coordinates are aligned into a common 3RLP pocket frame before comparison.

**Occupancy** is the fraction of analysed frames in which a region contains at
least one water oxygen. An occupancy of 0.8 means that region contains water in
80% of frames. It does not mean that the same water molecule stays there for 80%
of the trajectory, and it does not directly measure entropy or binding affinity.

W3, W249, and W286 label experimental water positions: W3 is taken from the 3RLP
crystal structure, and W249 and W286 from 3RLQ. These labels identify locations,
rather than water molecule IDs tracked through the simulations. The apo analysis
uses the 3RLP ligand pose as a spatial reference; that ligand is absent from the
apo simulation.

## 1. `target_water_occupancy.png`

![Occupancy around the three experimental water positions](target_water_occupancy.png)

**Caption.** Fraction of frames containing a water oxygen within 1.4 Å of each
fixed crystallographic water position, across four separate 50 ns simulations.
Blue, orange, green, and red bars represent APO, 3RLP, 3RLQ, and 3RLR,
respectively. Zero-height bars are invisible; the very small green W286 bar is
close to zero, rather than exactly zero.

**What we are seeing.** Each group asks whether water occupies the same
experimental location when the bound ligand changes. It does not require
HydraRank to detect a hydration-site center at that location.

| Experimental position | APO | 3RLP | 3RLQ | 3RLR |
|---|---:|---:|---:|---:|
| W3 | 0.469 | 1.000 | 0.000 | 0.000 |
| W249 | 0.984 | 0.994 | 0.991 | 0.000 |
| W286 | 0.377 | 0.546 | 0.003 | 0.000 |

W3 is occupied in every analysed 3RLP frame and becomes empty with the larger
3RLQ and 3RLR scaffolds. W249 stays almost continuously occupied through 3RLQ,
then becomes empty in 3RLR. These patterns support the displacement sequence
described in the validation summary.

W286 is less clear: it is already almost empty in the 3RLQ simulation, even
though its experimental reference position comes from that structure. These
trajectories therefore cannot validate a further W286 displacement from 3RLQ
to 3RLR. Near-zero simulated occupancy does not establish that the experimental
water is absent or thermodynamically unfavourable.

**How to interpret it.** The bars compare separate simulations, not successive
stages of one trajectory. They describe sampled occupancy without uncertainty
bars; there is only one trajectory per system. Figure 2 shows how much these
averages hide variation over time.

Source: [`paper_target_waters.csv`](../paper_target_waters.csv).

## 2. `target_water_block_convergence.png`

![Occupancy in five consecutive trajectory blocks](target_water_block_convergence.png)

**Caption.** Occupancy around W3, W249, and W286 in five consecutive blocks of
approximately 10 ns, using the same fixed 1.4 Å spheres as Figure 1. Each panel
shows one experimental location; each colored line shows one system. Block 1
covers approximately 0–10 ns, block 2 covers 10–20 ns, and so on to 50 ns.

**What we are seeing.** Unlike Figures 1 and 3, this horizontal axis represents
time. Each point is an average over its block, rather than an instantaneous
measurement or a cumulative average from the start of the simulation.

- **W3:** APO occupancy changes substantially: 0.754, 0.787, 0.448, 0.067,
  and 0.288. The overall APO value of 0.469 therefore conceals a strong change
  during the trajectory. In contrast, 3RLP stays at 1.0 and both 3RLQ and 3RLR
  stay at zero.
- **W249:** APO, 3RLP, and 3RLQ remain close to 1.0 in every block; 3RLR stays
  at zero. The three upper lines largely overlap.
- **W286:** APO fluctuates around 0.33–0.43 and 3RLP around 0.48–0.60.
  3RLQ stays nearly empty, and 3RLR stays completely empty.

Some lines are hidden by others. For example, the red 3RLR line covers the green
3RLQ line at zero in the W3 panel. This is overlap, not missing data.

**How to interpret it.** Large block differences, especially for apo W3, warn
that the estimated occupancy depends on which part of the trajectory is used.
Flat curves show stability over the sampled interval, but do not prove equilibrium
or reproducibility. A region that remains empty could also be poorly sampled;
independent simulations are needed to test that possibility.

Source: [`block_target_occupancy.csv`](../block_target_occupancy.csv), selecting
the rows with `radius_A = 1.4`.

## 3. `apo_site_occupancy_series.png`

![Occupancy at the same 19 apo-derived centers across the four systems](apo_site_occupancy_series.png)

**Caption.** Occupancy of 19 hydration-site centers detected in the apo
simulation, measured at those fixed centers in all four systems. Water oxygens
are assigned to their nearest apo center if it is within 1.0 Å. Red highlights
site 9 near W3; blue highlights site 2 near W249; dashed green highlights site 4,
the nearest detected apo site to W286. Gray lines show the other 16 apo sites.
The horizontal axis lists separate simulations, not time.

**What we are seeing.** This plot asks whether sites identified without a bound
ligand remain hydrated when each ligand is present. The site centers are held
fixed; the plotted lines do not connect separately detected holo sites. Site
numbers are identifiers, not ranking positions.

The red site is 0.32 Å from W3: its occupancy rises from about 0.41 in APO to
nearly 1.0 in 3RLP, then falls to zero in 3RLQ and 3RLR. The blue site is 0.46 Å
from W249: it remains substantially occupied through 3RLQ and is empty in 3RLR.
Gray lines show that other regions can remain hydrated or change occupancy as
the ligand changes.

**Why the green line needs care.** Apo site 4 lies approximately 1.50 Å from
W286, outside the 1.4 Å crystal-to-site matching threshold. It is a nearby site,
not a validated W286 correspondence. Its high APO/3RLP occupancy describes
water around the site-4 center, not around the experimental W286 position. The
dashed line and open circles mark this distinction. Its zero values overlap
the red line in 3RLQ and 3RLR.

**Why this differs from Figure 1.** Figure 1 uses a 1.4 Å sphere centered on
each experimental position. This figure uses a 1.0 Å assignment radius around
detected apo centers, with each oxygen assigned to at most one center. Changing
the center and radius changes the measured occupancy. For example, high
occupancy at site 4 can coexist with much lower occupancy around W286.

The distances in the legend describe spatial correspondence; they are not bond
lengths, ligand distances, or measures of water residence time. Occupancy alone
also does not determine the HydraRank displacement ranking.

Source: [`common_apo_sites.csv`](../common_apo_sites.csv).

## 4. `sensitivity_site_count.png`

![Number of detected apo sites under alternative clustering settings](sensitivity_site_count.png)

**Caption.** Number of apo hydration sites detected when one analysis setting
is varied at a time. Blue bars vary the site radius, orange bars vary the
density factor, and green bars discard the first 5 or 10 ns. All other settings
remain at their defaults: radius 1.0 Å, density factor 2.0, and no discarded
time. The dashed horizontal line marks the default result of 19 sites.

**What we are seeing.** This figure reruns site detection on the apo data;
it does not compare the four ligand systems.

| Setting varied | Tested values | Detected sites |
|---|---|---|
| Site radius | 0.8, 1.0, 1.2 Å | 26, 19, 10 |
| Density factor | 1.5, 2.5 | 24, 13 |
| Initial time discarded | 5, 10 ns | 19, 19 |

The radius controls the spatial scale of clustering. A smaller radius resolves
more separate sites here; a larger radius gives a coarser map with fewer sites.
The density factor sets the required observation count relative to the count
expected for bulk-density water in a sphere of that radius. Increasing it makes
site detection stricter. Discarding early frames means detecting sites again
using only the remaining trajectory, rather than deleting existing sites.

**How to interpret it.** The map's level of detail depends strongly on the
clustering settings. More sites do not automatically mean a better map. At
1.2 Å, W3 is no longer recovered as a distinct matching site, as recorded in
the sensitivity table. However, changing clustering settings does not change
the fixed experimental-position occupancy measurement in Figure 1.

Keeping 19 sites after discarding early frames establishes stability of the
count, not identical centers, occupancies, entropies, or ranks. In particular,
it does not resolve the apo W3 variation shown in Figure 2. The hydrogen-bond
penalty sensitivity listed in the CSV concerns ranking and is not plotted here.

Source: the `clustering` rows in [`sensitivity.csv`](../sensitivity.csv).

## Further context

The figures are generated by [`run_validation.py`](../run_validation.py).
[`summary.md`](../summary.md) connects them to crystal-water recovery, apo
ranking, and ligand coverage. Together they test structural correspondence,
occupancy changes, and sensitivity to analysis choices. They do not establish
prospective binding-affinity predictions; the apo system was also initialized
with crystallographic waters, which limits how independently its map tests
recovery of those locations.
