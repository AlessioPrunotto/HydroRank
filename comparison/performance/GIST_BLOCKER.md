# SSTMap GIST preflight blocker

The installed SSTMap 1.1.4 environment imports both analysis classes, but both
`gist_full` and `gist_entropy` terminate with SIGSEGV on the first analysed frame
of the prepared 3RLP input. They are excluded from repeated performance measurements.
Raw exit codes, logs, input hashes, and compiled-extension hashes are retained in
the pilot evidence bundle.

A separate one-frame run with Python's fault handler locates the crash at
`GridWaterAnalysis._process_frame`, line 315 of the installed module, while calling
`_sstmap_ext.assign_voxels`. This reproduces in entropy-only mode, so it does not
require hydrogen-bond or interaction-energy calculations. Python exception handling
cannot capture this native crash; absence of `adapter.json` is expected for this failure.

The retained extension source uses `PyArray_FromDims` to create the per-water voxel
record and dereferences its return value without checking for an allocation/API
error. This is a candidate cause, not a confirmed diagnosis of the installed binary.
A correction would need an isolated rebuild, explicit NumPy C-API error handling,
and validation against independently computed voxel membership.

The source also hardcodes 0.5 Å in its voxel-index calculation. Supporting other
resolutions requires a documented API/source change and boundary tests; passing
another Python grid resolution alone would not make assignment correct.

Before running GIST's energy/hydrogen-bond mode, also validate its neighbor-index
logic for OPC: it compares the oxygen-neighbor vector with the two-element voxel
record and subsequently indexes all atom IDs with neighbor-list positions. A
working import or repaired allocation is not sufficient validation of this path.

No SSTMap installation, extension, or scientific result has been patched by this
benchmark. Resolving these issues is a prerequisite for completing the planned
GIST comparison. The HSA pilot can still establish measured HSA workflow costs,
without supporting a claim about GIST performance.
