# Contributing

Contributions are welcome. Please open an issue before making a scientific or file-format
change so assumptions and compatibility can be discussed explicitly.

## Development setup

```bash
uv sync --python 3.11 --all-extras
```

Run the same checks as CI before submitting a change:

```bash
uv run ruff check src tests
uv run ruff format --check src tests
uv run pytest
```

Tests should use the compact in-memory systems in `tests/conftest.py`; do not add large
trajectory files. A bug fix should include a regression test, and new CLI behavior should
cover both human-readable and machine-readable or file output where relevant.

Scientific-method changes should document their assumptions, cite the underlying method,
include synthetic limiting-case tests, and report comparisons against an established tool
or reference calculation when possible.

## Style and compatibility

- Support the Python versions declared in `pyproject.toml`.
- Keep distances in ångström, times in picoseconds, and energies in kcal/mol.
- Preserve saved-observation compatibility or add an explicit migration.
- Update `CHANGELOG.md`, the README, and API documentation for user-visible changes.
