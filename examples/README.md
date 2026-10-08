# Software demonstration

Run the included HSP90 3RLP demonstration from the repository root:

```bash
uv run hydrarank analyse -c examples/hsp90_demo.yaml
```

`hsp90_demo.yaml` uses the topology and 101-frame trajectory in
[`data/hsp90_3rlp_demo/`](data/hsp90_3rlp_demo/README.md) and writes to `output/demo/`.
Input and output paths are resolved relative to the configuration file.

The 200 ps slice exercises the real-file software workflow. It is too short for
scientific interpretation of residence times, entropy, or displacement rankings.

The four full-trajectory HSP90 configurations are under
[`validation/hsp90/`](../validation/hsp90/README.md). They require production data
that must be supplied separately.
