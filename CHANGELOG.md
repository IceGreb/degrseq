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

- `end-counts`: paired-end `--min_aligned` counts both mates' positions, as nanoPARE does (DL-062); test added.

- `srna-shuffle`: shuffled sRNA sets for the cleavage null model (`--n_shuffles`, `--seed`,
  `--layout combined|per_set`), stats JSON. Implements `specs/srna-shuffle.md`.

- `cleavage-stats`: significant sRNA cleavage sites from real vs shuffled-sRNA step-1 tables
  (empirical p-values, Fisher, BH, window join), `--nanopare_compat`, stats JSON.
  Implements `specs/cleavage-stats.md`.

- `scale-factor`: the 5P-vs-BODY scale factor for EndGraph (`factor × 5P − BODY`), from gene regions
  (BED6) and ends bedgraphs; `--fraglen`, `--body_strand`, `--max_scale_factor`, `--write_genes`, stats JSON.
  Implements `specs/scale-factor.md` (DL-063). Deliberate differences from nanoPARE: no compat mode
  (DL-058); per-gene BODY counts are counts with shared positions split, not density × length (known
  issues 44, 45); last gene of a chromosome included (#43); gene length = region size (DL-059); empty
  totals are an error (DL-060); factor printed with 4 fixed decimals. On nanoPARE's test set:
  0.3871 / 0.1903 / 0.2103 vs the original's 0.3924 / 0.2100 / 0.2121. Fixture
  `tests/fixtures/scale_factor_test_set/` (1.5 MB, gzipped).
- `scale-factor --body_strand unstranded`: BODY plus + minus summed per nt as background for genes on both strands,
  for unstranded BODY libraries (DL-069). Default `reverse` unchanged (= nanoPARE's hard-coded flip).
- `kde-peaks`: 5′-end features of a 5P library (EndGraph): scale × 5P − BODY, Laplace kernel density, runs of
  positive density ≥ `--min_length`, density peaks, dominant 5P position, RPM filter; `--body_strand`
  (`reverse`/`forward`/`unstranded`), `--write_density`, `--nanopare_compat`, stats JSON. Implements
  `specs/kde-peaks.md` (DL-072, DL-073). Compat matches nanoPARE's golden EndGraph outputs on the test set
  (density byte-identical; features identical up to the order of rows with equal start).

- `cap-class`: capped / noncapped 5′-end features of a replicate group (EndClass, numeric part): per-nt
  replicate sums of 5P and cap tracks, trim to signal, rightmost peak, uG fraction split; `--min_uug`,
  `--min_reads`, `--decimals`, `--nanopare_compat`, stats JSON. Implements `specs/cap-class.md` (DL-087).
  Compat matches nanoPARE's EndClass outputs on the test set (sum bedgraphs byte-identical; peaks and
  capped/noncapped sets identical; 323 / 279). Fixture `tests/fixtures/cap_class_golden/` (415 kB).

### Differences from nanoPARE (`cap-class`, see the spec's quirk table)

- uG counted over the exact feature, no "+1 buffer"; ends and peak offsets 1 nt shorter (DL-074, known issue 15).
- `--min_uug 0` calls every feature capped; compat keeps the automatic cutoff (DL-080, known issue 50).
- Features keyed by chromosome, start, end and strand (known issues 52, 53).
- Replicate sums exact, written with 4 decimals (DL-085); uG sums always written (DL-086).
- Default mode stops when uG exceeds 5P in a feature (per-G cap track fed to default mode).
- Not ported: gene counts tables (DL-076), cap-mask bedgraphs (EndMask).

### Documentation

- `end-counts` spec: limitation added, a cap G over a genomic G looks templated (DL-079).

### Fixed (external code review of `cb08476`, 30/10)
- `end-counts`: an alignment past the reference end is an error, not a crash in
  `--write_coverage`; `XW` must be finite and ≥ 0; paired-end reads with `HI` on only some
  records pair mates by position; `kept_3prime_clip_only` no longer counts reads dropped by
  `--multimapper_mode unique`; mate sequences are kept apart in block identity; clean errors
  for missing inputs and output folders (also `srna-shuffle`); separate messages for
  `--min_aligned` and `--decimals`; the FASTA is closed.
- `end-counts`: read names that come back after another name are counted in stats
  (`read_names_repeated`) and warned about; counting is unchanged.
- `end-counts --write_coverage` (default) uses memory proportional to the intervals, not the
  chromosome length; multimapper rescue and coverage no longer scan every interval per
  reference. Output unchanged.
- All outputs are written as UTF-8.

### Differences from nanoPARE (`srna-shuffle`)

- Seeded (`--seed`, default 0): the same input and seed give identical output. The original was
  unseeded. Shuffles are statistically indistinguishable from the original's.
- Default output is one combined FASTA (`<prefix>.shuffled.fa`); `--layout per_set` writes one
  file per set, as the original.
- Multi-line FASTA accepted; header descriptions dropped from IDs; duplicate IDs and empty
  sequences are errors.

### Differences from nanoPARE (`kde-peaks`, see the spec's quirk table)

- Full kernel (the original drops the +W weight), density stops at the chromosome end, feature end exclusive
  (the original's is 1 nt short and leaves the last nt out of RPM and the dominant position), RPM over both
  strands' total (the original uses twice the plus total), no 2-decimal rounding of the difference before smoothing.
  Test set: replicable features 602 → 608 (DL-065, DL-066).
- A chromosome missing from `--lengths` is an error (the original crashes or drops features silently).
- No Gaussian kernel (DL-068); `--bandwidth` above 30 is an error (the original caps it silently).
- Smoothed density only with `--write_density` (DL-067).

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

### Differences from nanoPARE (`cleavage-stats`)
- Inputs are arguments, not fixed `Ath_annotations/` paths; class and family can come from
  `--srna_table` (name rule as fallback); unclassified sRNAs are skipped with a warning.
- Shuffled sites come as one pooled table per window (DL-050), so any number of sets works.
- An sRNA with no shuffled sites gets p = 1 and a warning; the original stops with an error.
- An sRNA may have several precursors (self:self removal checks all); the original stops on
  duplicate rows. sRNAs with real sites missing from `--srna_table` are counted and warned.
- A site found only in W50, or a window without sites, no longer stops the run (Q10).
- Labels: a site found in one window keeps its common name; `verified.` is `yes` if either
  window says yes; `interaction` uses the reported sRNA. Compat keeps R's `min()` and the
  W20 sRNA (Q3, Q9).
- Output order breaks adj-p ties by fold change, sRNA, transcript and site; compat keeps R's
  order (by `Transcript_Slice.site` text, then adj p).
- New per class: `significant_sites.tsv`, `window_tests.tsv` (every tested site per window with
  its p-values) and `shuffled_values.tsv` (the shuffled values behind them). Per-sRNA ECDF PDFs
  are not written; plotting is left to the pipeline, from these tables (DL-054, DL-055).
- A validated target also matches by its own sRNA name, so explicit families work
  (`miR2934` validates `miR2934-3p`); no change with the name rule alone.
- p-values are kept as the original (no +1; DL-051).
- adj p can differ from R in the 15th significant digit (scipy vs R chi-square tail).
