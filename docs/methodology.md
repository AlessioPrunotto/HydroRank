# Methodology and interpretation

This guide explains the choices behind HydraRank's current implementation and how
they affect interpretation. For commands and input examples, see the
[README](../README.md); for the mathematical description and references, see
[the manuscript](../paper/manuscript.md).

## Inputs and atom selections

Complete analysis requires atom identities, masses, molecular bonds, explicit
water hydrogens, and periodic-box information. Water residues are selected using
a configurable residue-name list, including OPC names, so selection does not
depend on a library's built-in water-name convention. An empty selection reports
the residue names present to help diagnose configuration errors.

Atom roles are assigned primarily by mass because oxygen and hydrogen names vary
between force fields. Three-, four-, and five-site water models are supported;
massless virtual sites are excluded from molecular geometry. Names provide a
fallback for selection-only operations when masses are unavailable, but analyses
that require masses reject missing mass information.

An explicit ligand selection is preferred. Automatic detection excludes recognized
waters and ions and applies a heavy-atom count threshold. Multiple remaining
candidates raise an error rather than selecting one arbitrarily. This guards
against defining the hydration region around the wrong molecule.

## Periodic boundaries and the local reference frame

The implemented transformation stack first makes the solute whole and tracks its
periodic images, then centres the ligand or apo pocket anchor, and finally wraps
the other molecules by residue. Residue wrapping translates each water together,
rather than moving its atoms independently. Correct topology connectivity and
intact input water geometry remain necessary.

The first solute frame is reconstructed from connectivity. Forward frames use a
minimum-image update relative to the previous processed coordinates. This assumes
that motion between processed frames can be resolved by that convention; widely
spaced frames can violate the assumption.

Spatial cutoffs are checked against half the shortest box-vector length. This is
an input safeguard, not an exhaustive proof that every periodic image is correct
for every pocket geometry or changing cell. Inspect the prepared structure and
diagnostics when applying the workflow to a new system.

Frames are rigidly fitted to a local protein reference, normally using pocket
Cα atoms, with other backbone atoms as a fallback. Local fitting reduces blurring
from distant protein motions. Fitting the ligand alone can give an unstable frame
for small or symmetric ligands and can transfer ligand motion into the protein
coordinates. Local fitting cannot remove genuine changes in pocket shape.

### What preprocessing QC establishes

QC collects the longest solute bond, local-fit RMSD, and number of shell waters
for every selected frame. The longest bond is measured without periodic-image
correction after the transformation stack. A maximum greater than 3.0 Å flags a
potential split solute and stops the unified analysis by default. A value of
exactly 3.0 Å passes. This uses the topology's recorded bonds: it does not validate
all expected chemical connectivity or independently prove that the solute is whole.

A maximum fit RMSD above 2.5 Å produces a warning. Residual motion can broaden
water distributions and affect inferred ordering, but high RMSD can also reflect
real pocket flexibility. The warning requires structural inspection rather than
automatically rejecting the trajectory. Shell-water counts are descriptive;
there is no acceptance threshold for them.

Passing these checks supports interpretation of preprocessing, but does not
establish sampling convergence, correctness of the force field, or accuracy of
the downstream thermodynamic interpretation.

## A fixed hydration region, including apo analysis

The hydration shell is defined around reference ligand heavy-atom coordinates in
the aligned frame, with a default cutoff of 5.0 Å. It remains fixed while the
trajectory is analysed. A region following the instantaneous ligand would mix
changes in water behaviour with changes in the measurement region.

For apo trajectories, an external bound structure supplies the reference pocket
and ligand pose. Pocket atoms must map uniquely onto the apo topology by residue
identity and atom name. The reference ligand defines a spatial region only: it is
not added to the simulation and contributes neither hydrogen bonds nor enclosure.
The map remains conditional on that reference pose and on the apo conformations
actually sampled. A holo map alone cannot reconstruct waters excluded by its ligand.

## Density-peak site detection

HydraRank uses deterministic density-peak clustering with spatial-bin candidates
and KD-tree neighbour searches. At each step it chooses the candidate with the
most remaining observations inside the site radius, assigns those observations,
and reports their centroid as the site center. Assigned observations are removed
before selecting the next peak. Observations outside accepted sites remain
unassigned. A configured site-count limit can also end detection.

The default radius is 1.0 Å. The default minimum observation count is:

```text
max(1, ceil(density_factor × bulk_density × (4π/3) × radius³ × n_frames))
```

Here `density_factor` defaults to 2.0 and `bulk_density` to 0.0333 molecules/Å³.
Thus the default threshold is twice the expected bulk population, rather than
simply bulk density. At a 1.0 Å radius, bulk corresponds to approximately 0.14
oxygen observations per frame and the threshold to approximately 0.28. If at most
one oxygen occupies the sphere at a time, this is roughly 28% occupancy. Different
water molecules can contribute over time; the threshold does not imply packing
multiple molecules into a single water-sized region.

The threshold selects localized enrichment. It is a methodological choice, not a
significance test, and regions below it may still have relevant solvent behaviour.
Radius and density factor affect the number and identity of detected sites. Site
numbers identify clusters, not rank positions or persistent water molecules.

## Structural and dynamical descriptors

Occupancy is the fraction of frames containing assigned water observations at a
site. It does not establish that one molecule persists throughout those frames.
Residence episodes instead follow unique topology water identities, joining visits
according to the configured maximum frame gap. Episodes touching trajectory
boundaries are tracked as censored. Mean episode duration is not a kinetic rate;
exchange faster than the saved frame spacing is unresolved.

Water–solute hydrogen bonds use topology-based donor hydrogens, candidate polar
atoms classified by mass, a default heavy-atom distance of 3.5 Å, and a minimum
donor–hydrogen–acceptor angle of 130°. This economical chemical model does not
replace protonation-aware perception or interaction-energy calculations.

Enclosure counts nearby solute heavy atoms within 5.0 Å by default. It is a
geometric descriptor, not an interaction energy, and is not used in the current
ranking score or category rules.

## Excess-entropy proxies

Translational and orientational entropy are estimated separately using
nearest-neighbour estimators. Values are expressed per water in units of the gas
constant relative to a bulk reference: uniform number density for translation
and uniformly random orientations for rotation. Negative values indicate ordering
relative to that reference. Their sum is converted to `−TΔS` in kcal/mol; positive
values describe an estimated entropic ordering cost.

The orientational estimator treats water's two hydrogens as indistinguishable.
A half turn about the molecular dipole bisector exchanges them. Both neighbour
distances and the reference rotational measure account for this symmetry, avoiding
dependence on arbitrary hydrogen labels.

These are first-order proxies. Higher-order water correlations are omitted, and
the observations pooled from a trajectory can be temporally correlated. A value
near zero indicates little estimated excess ordering under this model; it does
not establish that water displacement has no thermodynamic benefit or cost.
Interaction energies and other binding contributions remain relevant.

An estimator returns `nan` with fewer than ten observations or fewer than ten
valid positive nearest-neighbour distances. This prevents estimates from extremely
small or degenerate samples, but ten valid observations do not establish accuracy,
statistical independence, or convergence.

## Heuristic ranking and categories

The implemented score is:

```text
score = −TΔS − hbond_penalty × mean_water_solute_hydrogen_bonds
```

The default penalty is 1.0 kcal/mol per mean hydrogen bond. It is an uncalibrated
weight representing the need to consider contacts that a displacing ligand may
have to replace. It is not a measured bond-breaking free energy. Enclosure is
reported alongside the score but does not enter it.

Categories use separate thresholds. With the defaults, a site with an ordering
cost of at least 2.0 kcal/mol is `displaceable` if its mean solute hydrogen-bond
count is below 1.0, and `replace-hbonds` otherwise. Remaining sites, including
those with undefined ordering cost, are labelled `bulk-like`. That last label
can therefore mean insufficient evidence rather than demonstrated bulk behaviour.
Changing the penalty changes scores and rank order, but not category labels.

Categories provide an interpretation aid, and scores prioritize inspection.
Neither establishes the net free energy of displacing a water or changing ligand
affinity. The model omits explicit interaction energies, solvent correlations,
ligand strain, protein reorganization, and the interactions formed by a proposed
substituent. A favourable score is not sufficient evidence that a chemical change
will improve binding.

## Sampling and current validation

Input warnings identify fewer than 500 selected frames or effective frame spacing
greater than 20 ps; unknown time spacing is also reported. These are screening
rules, not system-specific convergence assessments. Increasing frame count by
saving more frequently does not necessarily increase independent information.
There is no universal trajectory length that guarantees reliable site statistics.

Inspect trajectory blocks for changes in occupancy, entropy, and site identity,
and compare independent simulations where available. Stable block averages do
not by themselves establish equilibrium or replica reproducibility.

The [HSP90 validation](../validation/hsp90/summary.md) provides retrospective
structural evidence for W3 and W249 displacement, while documenting incomplete
W286 sampling, drifting apo W3 occupancy, and parameter sensitivity. It uses one
trajectory per state, and the apo system began with crystallographic waters.
The [SSTMap comparison](../comparison/3rlp/report.md) assesses agreement on a shared
trajectory. These results support the tested analysis use case, while leaving
independent-replica reproducibility and prospective predictive performance open.
