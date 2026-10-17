# AGENTS.md

Notes for coding agents and contributors. Structure adapted from
[mgnifam](https://github.com/vagkaratzas/mgnifam) (MIT). This file covers the rules that will bite you.

## Ground rules

- **Clean-room: write from `specs/`, never from nanoPARE source.** nanoPARE
  (https://github.com/Gregor-Mendel-Institute/nanoPARE, commit `6264d6f`) is GPL-3. This
  package is MIT because it is written from behavioural specifications, not translated line by
  line. So, unlike mgnifam's `reference/`, **the original scripts are not vendored here**. The
  behavioural oracle is the set of golden outputs produced by running the original (see
  Testing). Each module docstring states which spec it implements.
- **Never add a dependency without updating the design-side dependency list.** Allowed:
  numpy, pysam, scipy (upper-bounded on purpose: they decide numeric output, so parity
  holds only for the resolved set in `uv.lock`). Interval operations (sort, merge, closest,
  subtract) belong in bedtools processes in the pipeline, not here.
- Use `uv lock --check && uv sync --frozen` before trusting a test result.
- Commits and pushes are the human's call.

## Naming contract

Subcommand `end-counts` ↔ module `src/degrseq/end_counts.py` ↔ spec
`specs/end-counts.md` ↔ pipeline module `modules/local/end_counts/` ↔ process `END_COUNTS`.
`tests/test_cli.py` enforces the first two. Each module exposes `main(argv)` and owns its
parser, with `prog="degrseq <name>"`.

Flags are `--snake_case` (`--n_shuffles`, `--nanopare_compat`), matching Nextflow `params`,
so a pipeline parameter and the flag it feeds have the same name.

## Behaviour vs nanoPARE

- **Default = fixed behaviour; `--nanopare_compat` = the original's quirks.** Every upstream
  bug a subcommand touches is listed in its spec under "Quirks" with a decision (fix /
  reproduce / fix + compat). Golden tests run in compat mode and must match. Default-mode
  tests must differ from golden output by exactly what each fix predicts.
- Deliberate differences go in `CHANGELOG.md` and in the spec, never only in code.

## Reproducibility

- Outputs must be byte-identical across runs, `--cpus` values and `PYTHONHASHSEED`. Sort
  explicitly before writing; never depend on dict or set iteration order for output.
- Anything random takes a `--seed` (e.g. `srna-shuffle`), and the seed is written to the stats file.

## Stats files (MultiQC contract)

Each subcommand that produces QC numbers writes `<prefix>.<subcommand>.stats.json` with a
`schema_version` field. A MultiQC module will parse it, so its shape is a contract: bump
`schema_version` on any breaking change. Record scientific parameters only (no paths, no `--cpus`).

## Testing

- `uv run pytest`. Unit tests use small handmade inputs placed in `tests/fixtures/`, each
  built to isolate one rule from the spec.
- Golden tests compare against outputs of the original nanoPARE run on its own test set
  (Arabidopsis chr2 window). Only the *outputs* are committed as fixtures (keep each under
  8 MB; `pre-commit` enforces it), plus a README recording the command, interpreter and
  commit that produced them.
- When you change scientific behaviour, the honest check is a diff against golden output,
  not a green test suite.

## Adding a subcommand

1. Approve its spec (`specs/<name>.md`, from `specs/TEMPLATE.md`).
2. Add `src/degrseq/<name_with_underscores>.py` with `main(argv)`.
3. Register it in `cli.COMMANDS`.
4. Add unit tests for each rule and edge case in the spec, then the golden test.
5. Update `CHANGELOG.md` and the pipeline's local module (replace its stub).
