# degrseq

Subcommands for the [rnadegradome](../nf-core-rnadegradome) Nextflow pipeline: mapping RNA
5′ ends, classifying them as capped or noncapped, and detecting small-RNA-directed cleavage
sites. Reimplements the custom logic of nanoPARE (Schon et al. 2018, Genome Res
28:1931–1942) from behavioural specifications.

**Status:** pre-alpha. Package skeleton only; no subcommands yet. See `specs/` and `AGENTS.md`.

## Development

```bash
uv lock --check && uv sync --frozen
uv run pre-commit install
uv run pre-commit run --all-files
uv run pytest
```
