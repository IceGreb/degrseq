# `end-counts`: behavioural spec

Status: **approved** (18/10/2026) · Replaces (nanoPARE 6264d6f): `scripts/python_scripts/readmapIO.py`, called at `endMap.sh` l.472–509; BODY coverage merge (`bedgraph_combine.py`, l.531)
Pipeline module: `modules/local/end_counts/` · Design-note section: §2.5 · Decisions: `package_build/degrseq_decision_log.md` DL-007–DL-017, DL-029

## Purpose
Turns aligned reads into per-nucleotide counts of RNA 5′ ends. For 5P libraries it also counts
reads carrying the cap signature (an untemplated G before the 5′ end). Reads that align equally
well to several places are shared out between them. The ends tracks feed peak calling (EndGraph),
cap classification (EndClass), isoform choice and masking (EndMask) and the scale factor. The cap
track feeds cap classification and masking.

## Inputs
| Argument | Format | Required | Notes |
|---|---|---|---|
| `--bam` | BAM from STAR, all alignments of a read adjacent (STAR `Unsorted`, or `samtools sort -n`) | yes | Header `SO:coordinate` → error (rule 1). Uses tags `NH` (if present), `XW` (if present). Read names are opaque grouping keys; a UMI suffix from umi_tools needs no handling (DL-029) |
| `--fasta` | genome FASTA with `.fai` | yes (5P) | read via `pysam.FastaFile`; every BAM reference must exist with the same length, else error |
| `--library` | `5P` or `BODY` | yes | switches rules 6–9 (5P only) and rule 17 (BODY only) |
| `--prefix` | string | yes | output file prefix |

## Outputs
All bedgraphs are 0-based half-open: `chrom  start  end  value`. Sorted by chrom (bytewise, as `LC_ALL=C sort`) then start, which is the order `bedtools sort` produces.

| File | Content | Consumer |
|---|---|---|
| `<prefix>.ends.plus.bedgraph`, `<prefix>.ends.minus.bedgraph` | 5′ ends (rule 21); one row per nucleotide, never merged | EndGraph, EndClass, EndMask, scale-factor |
| `<prefix>.cap.plus.bedgraph`, `<prefix>.cap.minus.bedgraph` | 5P only: cap-positive reads at their 5′ end (rule 22); one row per nucleotide | EndClass (`cap-class`), cap masking |
| `<prefix>.cov.plus.bedgraph`, `<prefix>.cov.minus.bedgraph` | only with `--write_coverage`, 5P: per-nt coverage; runs of equal value merged | none in the pipeline (QC) |
| `<prefix>.cov.bedgraph` | only with `--write_coverage`, BODY: coverage of both strands summed; runs merged | none in the pipeline (QC) |
| `<prefix>.end-counts.stats.json` | see Stats file | MultiQC |

Ends and cap rows are one nucleotide each, because downstream `bedtools map -o sum` sums rows, not per-nt values.

## Rules

### Reading the BAM
1. **Grouping.** Consecutive records with the same query name form one read. A header with `SO:coordinate` is an error. If `NH` is present and a group's number of distinct loci (after rule 4) differs from `NH`, the read is still processed with the loci present, and counted as `nh_mismatch`. This is expected after UMI dedup (pipeline G1).
2. **Records used.** Unmapped (0x4), supplementary (0x800), QC-fail (0x200) and duplicate-flagged (0x400) records are ignored. Primary and secondary (0x100) records are both loci (DL-013).
3. **Strand.** SE: a reverse (0x10) record is on `−`, else `+`. PE: the locus strand is mate 1's strand (rule 17).
4. **Locus identity.** A locus is keyed by (reference, strand, 5′-end position after rule 7). Records with the same key within one read count once (DL-011). Compat: keyed by alignment start only, keeping the first record (#24).

### Aligned positions and soft-clips
5. **CIGAR positions.** Aligned positions are those consumed by `M`, `=`, `X`. `D` and `N` advance the reference without adding positions. `I`, `S`, `H` add none.
6. **Reclaim (5P only).** For each soft-clip adjacent to the alignment, walk outward from the alignment edge base by base. A clipped base is reclaimed (treated as templated and aligned) while it equals the genome base at the corresponding position (case-insensitive; `N` never matches). Stop at the first non-matching base or the reference end. The bases left over, from that point to the read end, are **untemplated**.
7. **5′ end.** `+`: the leftmost aligned-or-reclaimed position. `−`: the rightmost. The 5′-side untemplated bases are those beyond it (read start for `+`, read end for `−`, in read orientation).
8. **Untemplated non-G at the 5′ side → locus unusable (5P).** If any 5′-side untemplated base is not `G` (after reverse-complementing for `−`), the locus is unusable (DL-009). Counted as `loci_untemplated_non_g`; a read with no usable locus left is counted under `discarded_untemplated_non_g`.
9. **Untemplated bases at the 3′ side only (5P).** These are ignored; the locus stays usable (DL-010). Reads with such a locus are counted as `kept_3prime_clip_only`. Compat: the locus is unusable.
10. **Cap signature (5P).** A usable locus is cap-positive if it has ≥ 1 untemplated 5′-side `G` (DL-008). The run length (number of untemplated G's) goes into the stats histogram (rule 26).
11. **Minimum length.** A locus is usable only if its aligned positions (rule 5) plus reclaimed bases (rule 6) total ≥ `--min_aligned` (default 16), for unique and multimapped reads alike (DL-007). Counted as `loci_too_short`; a read with no usable locus is `discarded_too_short`. Compat: unique reads need `> 16` (reclaim-aware); multimapper loci need `≥ 16` CIGAR positions only.
12. **Discard precedence.** A read with no usable locus is counted under exactly one reason, checked in this order: `too_short`, `untemplated_non_g`, `3prime_clip_only` (compat only).

### Weights
13. **Read weight.** `w` = the `XW` tag (float) if present on the read's records, else 1. The same `w` applies to ends, cap and coverage (DL-016). Records of one read with different `XW` values → error.

### Unique reads
14. **Unique read** = exactly one usable locus. It adds `w` to the ends track at its 5′ end, `w` to the cap track at its 5′ end if cap-positive (5P), and `w` to coverage at every aligned and reclaimed position.

### Multimappers
15. **Multimapped read** = ≥ 2 usable loci. Processed after all unique reads.
16. **Mode** (`--multimapper_mode`, DL-012):
    - `rescue` (default): rules 18–20;
    - `uniform`: share `1/k` to each of the k loci;
    - `unique`: multimapped reads are dropped (counted `multimapped_dropped`).
17. **Paired-end.** A locus is a mate-1 record plus the mate-2 record with the same `HI` tag. The 5′ end, strand and cap come from mate 1. Coverage is added from both mates. A locus without an aligned mate 1 is unusable (DL-017).
18. **Rescue weight of a locus** = the sum of current same-strand coverage over the locus's aligned and reclaimed positions.
    - BODY only: if every locus of the read has same-strand weight 0, opposite-strand coverage over the same positions is used.
    - If all weights are still 0: equal split (counted `multimapped_naive`), else counted `multimapped_rescued`.
19. **Shares.** `share_i = weight_i / Σ weight_j`. A locus with weight 0 gets share 0 unless the split is equal (rule 18). Each locus receives `w × share_i` on the ends track (5′ end), cap track (if that locus is cap-positive) and coverage (its positions).
20. **Order (D1 decided, DL-031: the original's order).**
    - **Tier** = number of loci the read has before the length filter (rule 11). Tiers are processed in ascending order.
    - **Block** = reads of one tier with the same leftmost alignment start (integer, chromosome not included) and the same sequence (mate 1 + mate 2). A block is shared once, with `w` = the sum of its reads' weights. Default: reads join a block only if their locus sets are also identical (D3, DL-033). Compat: the block takes the loci of its **last** read (G1, verified 17/10).
    - **Passes:** visit blocks by leftmost start, ascending on odd passes and descending on even passes. Each placed block updates coverage before the next one is visited. A block whose loci all have weight 0 (rule 18) stays for the next pass. Repeat until a pass places nothing; then split the remaining blocks equally (rule 18).
    - **Ties** (blocks with the same leftmost start): bytewise order of the block sequence, then of the sorted locus list (chromosome, start, strand); reversed on descending passes (D2, DL-032; D3, DL-033). Compat: first-seen order in the BAM.
    - Compat: identical, plus the G1 details as audited.

### Writing
21. **Ends track:** the per-nt sums of rules 14 and 19.
22. **Cap track (5P):** the per-nt sums of the cap contributions, at the 5′-end coordinate (DL-008). Compat: per untemplated G, +1 unweighted for unique reads, at each G's own position (1…k nt upstream); multimappers weighted by share (#8, #26).
23. **Rounding.** Values are rounded to `--decimals` (default 4; compat 2) with Python `round()`. Rows whose rounded value is 0 are not written. Values are printed without trailing zeros (`1`, `0.5`, `0.3333`) (DL-014). Compat format: 2 decimals, Python `round()` repr (`1.0`), equal adjacent values merged into one row (G4).
24. **Determinism.** Output is byte-identical across runs and `PYTHONHASHSEED`. Accumulation is single-threaded, in the order defined by rules 14 and 20, in float64.
25. **Coverage files:** only with `--write_coverage` (DL-015).

### Stats
26. See Stats file. Every read in the BAM is counted under exactly one fate.

## Edge cases
- **Empty BAM, or no usable reads:** empty bedgraphs (0 bytes) and a stats file with zeros; exit 0.
- **Reference in the BAM but not in the FASTA, or a length mismatch:** error before reading any records.
- **Soft-clip reaching the reference start or end:** reclaim stops there (rule 6); the remaining bases are untemplated.
- **5′ end at position 0:** valid. The cap is at the 5′-end coordinate, so it never goes negative (in compat, a −1 position is dropped and counted).
- **Read entirely soft-clipped after alignment (`0M`):** cannot occur from STAR; error if seen.
- **`H` (hard clip):** ignored (no positions, no reclaim).
- **`NH` missing (non-STAR aligner):** grouping by name only; `nh_mismatch` not computed.
- **Multimapper whose loci all share one key after rule 4:** treated as unique.
- **All loci of a multimapper unusable except one:** treated as unique (rule 14) and counted in `unique_after_filter`.
- **`XW` = 0:** the read contributes 0 everywhere but is still counted under its fate.

## Quirks of the original
| # (known_issues) | Original behaviour | Decision | Expected diff vs golden |
|---|---|---|---|
| #7 | Unique reads need > 16 positions (reclaim-aware); multimapper loci ≥ 16 CIGAR-only | fix + compat (DL-007) | TG3: up to 825 more unique alignments; multimapper loci passing only through reclaim |
| #8, #26 | Cap counted per untemplated G, +1 unweighted for uniques, at each G's position | fix + compat (DL-008) | cap moved to the 5′-end coordinate; multi-G reads count once; uniques weighted (identical when no `XW`) |
| #24 | Loci keyed by start position only | fix + compat (DL-011) | honeybee 5P sample, 5 % subsample: 101/282,142 reads (0.04 %) change from unique to multimapped |
| — | Reads with 3′-only untemplated bases discarded | fix + compat (DL-010) | TG3: ~0.6 % of unique alignments kept |
| — | Values rounded to 2 dp; shares < 0.005 vanish | fix + compat (DL-014) | more rows; values at 4 dp |
| — | Multimapper order depends on BAM order and integer positions across chromosomes, in alternating passes | fix + compat (DL-012, G1) | small redistribution among multimapper loci |
| #25 | Splice-junction tally uses a stale chromosome | not ported (option unused, DL-015) | — |
| #9 | PE 5P aligned end-to-end | pipeline (STAR args), not here | — |

## Parameters
| Flag | Default | Original equivalent |
|---|---|---|
| `--library` | — (required) | `-R 5P` / `-R BODY` plus `--softclip_type 5p --untemp_out G` (5P), `--allow_nonstranded` (BODY) |
| `--min_aligned` | 16 | `--minmatch 16` |
| `--multimapper_mode` | `rescue` | always rescue, `--allow_naive` |
| `--decimals` | 4 | 2 (hard-coded) |
| `--write_coverage` | off | coverage always written |
| `--nanopare_compat` | off | — |

Not ported: `--secondary` (always on), `--weighted` (automatic when `XW` present), size-class, BED, junction and multimapper-dump outputs (DL-015).

## Stats file
`<prefix>.end-counts.stats.json`, `schema_version: 1`.
- **Parameters:** `library`, `min_aligned`, `multimapper_mode`, `decimals`, `nanopare_compat`.
- **Read fates** (they sum to `reads_total`):
  - `unique`, `unique_after_filter`;
  - `multimapped_rescued`, `multimapped_naive`, `multimapped_uniform`, `multimapped_dropped`;
  - `discarded_too_short`, `discarded_untemplated_non_g`, `discarded_3prime_clip_only` (compat);
  - `unmapped` (reads whose records are all unmapped).
- **Loci:** `loci_too_short`, `loci_untemplated_non_g`, `loci_key_collisions` (compat).
- **Other:** `kept_3prime_clip_only`, `nh_mismatch`.
- **5P only:** `cap_positive_reads` (weighted sum) and `untemplated_g_run_length` histogram `{0, 1, 2, 3, ≥4}` over usable loci of unique reads. For GenXPro reads after normalisation this checks that the TSO GGG was removed.

## Tests
- **Unit**, small hand-made BAM + FASTA per rule:
  - rule 1: coordinate-sorted header rejected; NH mismatch counted;
  - rule 2: flag filtering;
  - rule 4: two loci with the same position on different chromosomes (default 2 loci; compat collision);
  - rule 6: reclaim stops at the first mismatch, at the reference end, and on `N`;
  - rule 7: 5′ end on both strands, with and without reclaim;
  - rule 8: `GA` 5′ clip → unusable; `GGG` clip → cap-positive;
  - rule 9: 3′-only clip kept (compat: discarded);
  - rule 11:
    - a 16-nt unique read (default kept, compat discarded);
    - a two-locus 16-nt read (kept in both modes);
    - a 15 + 2-reclaimed multimapper locus (default kept, compat dropped);
  - rule 13: `XW` 0.5 applied to all tracks;
  - rules 18–20:
    - a two-locus read with 3:1 coverage → 0.75 / 0.25;
    - no coverage → 0.5 / 0.5;
    - BODY opposite-strand fallback;
    - two reads in one tier at different leftmost starts: the result follows the ascending/descending pass order;
    - a zero-coverage read is placed in pass 2 after a neighbour adds coverage;
    - identical reads form one block;
    - two blocks with the same leftmost start: same output in both BAM orders (compat: follows BAM order);
  - rule 17: paired-end 5′ end from mate 1;
  - rule 22: GGG read → cap 1 at the 5′ end (compat: 1 at each of −1, −2, −3);
  - rule 23: 1/3 → `0.3333`, 1/300 kept at 4 dp, dropped at 2 dp;
  - rule 24: byte-identical output under two `PYTHONHASHSEED` values;
  - every edge case above.
- **Golden (compat):**
  - a collaborator's honeybee 5P/BODY run (STAR BAMs and EndMap bedgraphs of the original); a patched original reproduced them byte-for-byte on one sample (13/10). Fixtures: an 8 MB subsample.
  - nanoPARE test set after Phase 5b.
  - Expected: exact match in compat mode.
- **Default vs golden:** the difference must be explained entirely by the quirk table.

## Open
- **D1: multimapper order: decided 17/10 (DL-031), the original's order** (rule 20). Background: `package_build/17_10_2026_end-counts_D1_multimapper_order.md`.
- **D2: tie order between blocks with the same leftmost start: decided 17/10 (DL-032), bytewise sequence order** (rule 20). The toy run showed the original's tie order is random between runs under its pinned Python 3.5.
- **D3: identical reads whose locus sets differ: decided 18/10 (DL-033), option (a)** (rule 20). In the original, reads with the same sequence and leftmost start form one block, and the block takes the **last** read's loci. An earlier read's other loci are silently replaced (toy run: a read aligned to chr3:101/501 was counted at chr3:101/901).
  - **When it happens:** STAR gives the same loci to the same sequence, so normally never. Per-alignment UMI dedup (pipeline, UMI G1) can remove a read at one locus but not another, so it can happen after dedup [I].
  - **Options:**
    - (a) a block requires the same sequence **and** the same locus set (recommended). Identical to the original whenever loci agree, which is the normal case; otherwise no read is moved to loci it never aligned to;
    - (b) reproduce the original (last read's loci).
  - Compat: (b) either way.
  - **Blocks:** rule 20's block definition and one unit test.
- **Gaps (G):** behaviour of the original needed only for `--nanopare_compat` and the golden test. Each needs a targeted audit (nanopare-port-audit) of `readmapIO.py` before compat is implemented. Default mode doesn't depend on them.
  - **G1: closed 17/10 [V]** by a toy-SAM run of the original (`answers/17_10_2026_multimapper_rescue_audit.md` §6): sum of same-strand coverage over the footprint; ascending first pass, alternating; immediate update; zero-coverage blocks deferred; chromosome-blind integer `minpos`; block uses its last read's loci; ties follow first-seen order on Python ≥ 3.6 and the hash seed on Python 3.5 (pinned in `OPTIONS_lmod`), so golden runs need Python ≥ 3.6 or a fixed `PYTHONHASHSEED`. Not run, read from code only: tier = locus count before the length filter (l.284 vs l.593) [I].
  - **G2:** whether a multimapper locus with a non-G untemplated 5′ base is dropped alone or discards the whole read.
  - **G3:** whether reclaim is contiguous from the alignment edge exactly as rule 6 states, including on `N`.
  - **G4 (partly closed 17/10 [V]):** the 5′-end and coverage files **merge runs of equal value** (3 adjacent ends of 1 → one row `chr1 1000 1003 1.0`); values printed as Python `round(x, 2)` (`1.0`, `2.67`); the untemp file has a `#chrom start end G` header. Still open: zero rows, and the sort order before `bedtools sort`.
  - **G5:** BODY soft-clips (STAR Local in the collaborator's run): ignored, or reclaimed?
  - **G6:** paired-end loci with only one mate aligned.
  - **G7:** the BODY coverage merge (`bedgraph_combine.py`, l.531): strand sum or other.
- **Implementation notes (not rules):**
  - ends and cap are kept sparse;
  - coverage is needed only where multimapper loci fall, so it can be restricted to those positions after a first pass;
  - multimapper records are held in memory until uniques are done (a large fraction of 5P reads in typical runs).
