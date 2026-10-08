# Changelog

All notable changes to this project are documented here. The format follows Keep a
Changelog, and releases use semantic versioning.

## Unreleased

### Added

- Apo-trajectory analysis using an external bound ligand pose as a non-interacting spatial
  reference for pocket definition and alignment.
- Journal-neutral scientific manuscript draft with a BibTeX library and structured
  placeholders for the planned HSP90 validation.
- Unified `hydrarank analyse` workflow with automatic cache validation, QC, ranking,
  human-readable output, provenance, and a standard results bundle.
- Heuristic hydration-site ranking using entropy and water-solute hydrogen bonds.
- CSV and PDB/mmCIF site exports.
- Optional occupancy, residence-distribution, convergence, and ranked-site plots.
- CLI progress reporting with an opt-out flag.
- Format-3 observation caches with globally unique water identities.
- End-to-end, CLI, hydrogen-bond, ranking, export, and plotting tests.

### Changed

- README status now reflects retrospective validation; detailed design rationale and
  interpretation limits are documented in `docs/methodology.md`.

- Orientational entropy now uses an exact symmetry-aware KD-tree search, reducing its
  memory scaling from quadratic to linear for long trajectories.
- Periodic unwrapping and residue wrapping are vectorised, allowing multi-gigabyte
  trajectories to remain streamed while avoiding per-frame Python loops.
- Hydration-site clustering uses compact spatial-bin candidates and length-only
  neighbour queries instead of materialising quadratic neighbour lists.
- Renamed the distribution, Python package, command, exception hierarchy, and generated
  format labels from `water-entropy` / `water_entropy` to HydraRank / `hydrarank`.
- Command-line arguments now override YAML configuration values.
- Replaced the obsolete demonstration dataset with a compact 200 ps 3RLP HSP90 fixture;
  full-length HSP90 trajectories remain the basis for scientific validation.
- Data contracts, cutoff checks, time sampling, and configuration validation are stricter.

### Fixed

- Ligand oxygen atoms are included in the heavy-atom hydration-shell reference.
- The bundled observation cache is compatible with the current loader.

## 0.1.0

- Initial preprocessing, hydration-site clustering, entropy analysis, and command-line tools.
