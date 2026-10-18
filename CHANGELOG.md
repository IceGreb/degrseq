# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Package skeleton: `degrseq` command dispatcher, spec template, CI, pre-commit.
- Release workflow: tags `v*` run CI, build, `twine check`, smoke-test the wheel and publish to
  TestPyPI via trusted publishing.
- `end-counts` (default mode; `--nanopare_compat` not yet implemented): per-nucleotide 5′ ends,
  cap track (5P), optional coverage, stats JSON. Implements `specs/end-counts.md`.

### Differences from nanoPARE (`end-counts`, see the spec's quirk table)

- Cap counted once per read at its 5′ end, weighted, instead of per untemplated G upstream.
- Same minimum length (`--min_aligned`, default 16) for unique and multimapped reads,
  reclaimed bases included.
- Reads with untemplated bases only at the 3′ side are kept.
- Loci keyed by (reference, strand, 5′ end), not start position only.
- Values at 4 decimals, one row per nucleotide (no merged runs).
- Ties between multimapper blocks broken by sequence, then locus list; identical reads
  with different locus sets are shared separately.
