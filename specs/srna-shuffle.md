# `srna-shuffle`: behavioural spec

Status: draft (13/10/2026) · Replaces (nanoPARE 6264d6f): `EndCut.pipeline/scripts/sRNA_shuffler.py`
Pipeline module: `modules/local/srna_shuffle/` · Design-note section: §6 EndCut

## Purpose
Builds the null model for cleavage calling. Each small RNA (miRNA, tasiRNA) is shuffled N
times, keeping its length and nucleotide composition. GSTAr then predicts target sites for
every shuffled set exactly as for the real set. For each real sRNA, `cleavage-stats` compares
the real sites' fold change and Allen score with the distribution from that sRNA's own shuffles.
A real target site only counts if it beats sequences of identical composition that cannot pair
specifically.

## Inputs
| Argument | Format | Required | Notes |
|---|---|---|---|
| `--fasta` | FASTA of sRNAs (DNA alphabet; the bundled `anno.mir.tas.fa` has 587 records, 19–24 nt, no spaces in headers) | yes | Multi-line FASTA accepted (original required strict 2-line records) |
| `--n_shuffles` | integer ≥ 1 | no, default 1000 | Original: 2nd positional argument |
| `--seed` | integer | no, default 0 | Original was unseeded |

## Outputs
| File | Format | Consumer |
|---|---|---|
| `<prefix>.shuffled.fa` (default layout `combined`) | All N sets in one FASTA, set-major (set 0 for every sRNA in input order, then set 1, …). Header `><id>_shuffled_<j>` | GSTAr (chunked by the pipeline), then the pred-sites conversion, which strips `_shuffled…` but must keep `j` as a column |
| `<prefix>.<j>.shuffled.fa` × N (layout `per_set`) | One file per set, as in the original (`<stem>.<j>.shuffled.fa`) | For compat / golden-style runs |
| `<prefix>.srna-shuffle.stats.json` | see Stats | MultiQC |

## Rules
1. Records are read in input order; the header ID is the text after `>` up to the first whitespace.
2. For each set `j` in `0 … n_shuffles−1`, and each record in input order, output one record:
   header `>` + ID + `_shuffled_` + `j`, sequence = a uniformly random permutation of the
   input sequence's characters (Fisher–Yates).
3. Length and per-character counts of every shuffled sequence equal those of its input.
   Case is preserved character by character.
4. Sequences are written on one line.
5. All randomness comes from one `random.Random(seed)` stream, consumed in the order of rule 2.
   Same input + same seed + same `--n_shuffles` → byte-identical output, in either layout.
6. A shuffle identical to the input sequence is kept, not redrawn. This matches the original
   and leaves the null distribution unbiased.

## Edge cases
- Empty sequence → error naming the record (the original wrote an empty sequence line).
- Duplicate IDs → error (shuffled IDs would collide downstream).
- Header with a description (`>id desc`) → ID only. The original appended `_shuffled_j` after
  the whole header line, so the suffix landed after the description.
- Characters outside `ACGTUN` (any case) → warning, shuffled as-is.
- `--n_shuffles 0` or negative → error.

## Quirks of the original
| # | Original behaviour | Decision | Expected diff vs original |
|---|---|---|---|
| Q1 | Unseeded `random.shuffle`: the null set is not reproducible | **fix**: `--seed` | Shuffles differ; only statistical equivalence is testable |
| Q2 | Python 2 only (`print >>`) | **fix** (Python 3) | — |
| Q3 | Requires 2-line FASTA; multi-line FASTA silently produces garbage | **fix**: real FASTA parser | None on the bundled FASTA |
| Q4 | Output stem = input name split at the first `.fa` (`anno.mir.tas.fa` → `anno.mir.tas`), which breaks for names containing `.fa` earlier (e.g. `my.fasta.fa`) | **fix**: explicit `--prefix` | File names only |
| Q5 | Only the last output file handle is explicitly closed | **fix** | — |

No `--nanopare_compat` mode: an unseeded shuffle can't be reproduced, so there is nothing to match.

## Parameters
| Flag | Default | Original equivalent |
|---|---|---|
| `--fasta` | — | argv[1] |
| `--n_shuffles` | 1000 | argv[2] |
| `--seed` | 0 | — |
| `--layout` | `combined` | `per_set` behaviour |
| `--prefix` | input stem | split at `.fa` |

## Stats file
`schema_version`, `n_srnas`, `n_shuffles`, `seed`, `n_identical_to_input` (count of shuffles equal to their input; useful for short or low-complexity sRNAs).

## Tests
- Unit: composition and length preserved for every record; header format; input order and set-major order; same seed → identical bytes in both layouts; different seeds → different output; multi-line FASTA; description in header; duplicate IDs and empty sequence raise; `n_identical_to_input` counted.
- Golden: none byte-level (Q1). Statistical check: on the bundled `anno.mir.tas.fa`, per-position nucleotide frequencies of shuffled sets match those of the input pooled over each sRNA (χ² test, informational only).
- Downstream compatibility: GSTAr runs on a combined-layout chunk, and the pred-sites conversion recovers sRNA name and set `j` for every hit.

## Open
- Whether the pipeline chunks the combined FASTA for GSTAr (design proposal) or keeps per-set files. That only affects `--layout`'s default.
