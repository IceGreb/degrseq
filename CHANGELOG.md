# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Package skeleton: `degrseq` command dispatcher, spec template, CI, pre-commit.
- Release workflow: tags `v*` run CI, build, `twine check`, smoke-test the wheel and publish to
  TestPyPI via trusted publishing.
- `end-counts`: per-nucleotide 5′ ends, cap track (5P), optional coverage, stats JSON.
  Implements `specs/end-counts.md`.
- `end-counts --nanopare_compat`: reproduces the original's quirks (spec quirk table, DL-045).
  Matches the golden fixture byte for byte (5P ends, cap, BODY ends, BODY coverage); golden test
  `tests/test_end_counts_golden.py`. Requires `--multimapper_mode rescue`; values fixed at 2
  decimals, so `--decimals` is an error (DL-047). New read fate `multimapped_lost_compat` (always 0 in default mode).
  Single-end only: paired-end input is an error (DL-046). With `--write_coverage` it warns about
  memory and builds coverage one chromosome and strand at a time (DL-047), skipping positions past
  the chromosome end (DL-048).
- `end-counts` golden fixture (`tests/fixtures/end_counts_golden/`): nanoPARE test-set subsample and the
  original's outputs on it. Pre-commit whitespace hooks skip `tests/fixtures/` so fixtures stay byte-exact;
  `.gitignore` ignores `golden/` and `scratch/` only at the repo root, so fixture folders are tracked.

- `srna-shuffle`: shuffled sRNA sets for the cleavage null model (`--n_shuffles`, `--seed`,
  `--layout combined|per_set`), stats JSON. Implements `specs/srna-shuffle.md`.

### Differences from nanoPARE (`srna-shuffle`)

- Seeded (`--seed`, default 0): the same input and seed give identical output. The original was
  unseeded. Shuffles are statistically indistinguishable from the original's.
- Default output is one combined FASTA (`<prefix>.shuffled.fa`); `--layout per_set` writes one
  file per set, as the original.
- Multi-line FASTA accepted; header descriptions dropped from IDs; duplicate IDs and empty
  sequences are errors.

### Differences from nanoPARE (`end-counts`, see the spec's quirk table)

- Cap counted once per read at its 5′ end, weighted, instead of per untemplated G upstream.
- Same minimum length (`--min_aligned`, default 16) for unique and multimapped reads,
  reclaimed bases included.
- Clipped bases are taken back from the alignment edge outwards on both ends, ignoring genome
  case; `N` never matches; bases past a chromosome end count as untemplated. None of these
  differs from nanoPARE on STAR output (no clipped base next to an alignment matches the genome).
- BODY multimapped reads are grouped by sequence and locus set; nanoPARE groups them by CIGAR
  and counts every read of a group at the last read's loci (21–26 % of BODY multimapped reads in
  the nanoPARE test set). `--nanopare_compat` reproduces it.
- Paired-end: the 5′ end always comes from mate 1, and a non-G clip on mate 1 discards the read.
  nanoPARE takes the outermost position over both mates and keeps mate 2; this can't occur with
  its published STAR settings. `--nanopare_compat` refuses paired-end input.
- Bedgraphs: no `0.0` rows in the cap track; a 5′ end at chromosome position 0 is kept; BODY
  coverage sums unrounded strand values.
- Only the 5′ clip is examined: reads with untemplated bases at the 3′ side (alone or together
  with a 5′ clip) are kept, and 3′ bases never count as cap signal. The original drops them.
- Loci keyed by (reference, strand, 5′ end), not start position only.
- Values at 4 decimals, one row per nucleotide (no merged runs).
- Ties between multimapper blocks broken by sequence, then locus list; identical reads
  with different locus sets are shared separately.
- Multimapped reads: alignments with a non-G extra 5′ base are removed before the read is
  shared, so a read with one usable alignment is never lost. The original computes the shares
  first and can drop such reads silently.
- Cap signal comes only from reads that are counted, and only from the 5′ clip. The original
  also adds it for reads too short to be counted, and for G's in a 3′ clip.
