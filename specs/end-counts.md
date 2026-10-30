# `end-counts`: behavioural spec

Status: **approved** (18/10/2026); default mode implemented 18/10, `--nanopare_compat` 29/10 (golden fixture matched byte for byte) · Replaces (nanoPARE 6264d6f): `scripts/python_scripts/readmapIO.py`, called at `endMap.sh` l.472–509; BODY coverage merge (`bedgraph_combine.py`, l.531)
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
4. **Locus identity.** A locus is keyed by (reference, strand, 5′-end position after rule 7). Records with the same key within one read count once (DL-011). Compat: keyed by alignment start only (#24), separately for primary and secondary records, primary loci listed first (G9). Colliding records (counted in `loci_key_collisions`) form one locus with the first record's reference and strand and the union of their positions; for a unique read each record passes or fails the clip checks on its own and only passing records add positions. Golden fixture: 127 collisions in 5P, mostly one read's alignments at one start with different introns.

### Aligned positions and soft-clips
5. **CIGAR positions.** Aligned positions are those consumed by `M`, `=`, `X`. `D` and `N` advance the reference without adding positions. `I`, `S`, `H` add none.
6. **Reclaim (5P only).** For each soft-clip adjacent to the alignment, walk outward from the alignment edge base by base. A clipped base is reclaimed (treated as templated and aligned) while it equals the genome base at the corresponding position (case-insensitive; `N` never matches). Stop at the first non-matching base or the reference end. The bases left over, from that point to the read end, are **untemplated**.
7. **5′ end.** `+`: the leftmost aligned-or-reclaimed position. `−`: the rightmost. The 5′-side untemplated bases are those beyond it (read start for `+`, read end for `−`, in read orientation).
8. **Untemplated non-G at the 5′ side → locus unusable (5P).** If any 5′-side untemplated base is not `G` (after reverse-complementing for `−`), the locus is unusable (DL-009; kept after G8, DL-035). Counted as `loci_untemplated_non_g`; a read with no usable locus left is counted under `discarded_untemplated_non_g`.
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
    - Compat: identical, plus the G1 and G9 details as audited:
      - a tier-≥ 2 read stays a multimapper even if only one locus passes rule 11 (it is shared with the others in its tier);
      - block key: 5P (tier, leftmost start, sequence); BODY (tier, leftmost start, CIGAR of the leftmost record). The block takes the last read's loci; its weight is the read count × the last read's `w`;
      - per block, in this order: drop loci with < `--min_aligned` CIGAR positions (union over the locus's records, before take-back; all dropped → read lost, counted `discarded_too_short`); rescue weights (rule 18); all 0 → deferred; shares; 5P: loci failing rule 8 or the compat clip checks (3′-only clip, or any non-G base in either clip) get share 0 (all 0 → `multimapped_lost_compat`, G2); a locus merged from several records (rule 4) gets share 0 if any one of them fails, unlike a unique read (`readmapIO.py` l.652–737); remaining shares renormalised;
      - a pass counts as having placed something only if every block of at least one leftmost start was placed;
      - the final pass runs once, ascending, ties first-seen: each block is shared by rule 18 if it now has coverage (earlier blocks of the same pass count), else equally.

### Writing
21. **Ends track:** the per-nt sums of rules 14 and 19.
22. **Cap track (5P):** the per-nt sums of the cap contributions, at the 5′-end coordinate (DL-008). Compat: per untemplated G, at each G's own coordinate, from both clips (#34): +1 unweighted for unique reads, recorded before the length check (#33); multimappers weighted by count × share × `w` (#8, #26).
23. **Rounding.** Values are rounded to `--decimals` (default 4; compat fixed at 2, `--decimals` is an error, DL-047) with Python `round()`. Rows whose rounded value is 0 are not written. Values are printed without trailing zeros (`1`, `0.5`, `0.3333`) (DL-014). Compat format: 2 decimals, Python `round()` repr (`1.0`), equal adjacent values merged into one row (G4).
24. **Determinism.** Output is byte-identical across runs and `PYTHONHASHSEED`. Accumulation is single-threaded, in the order defined by rules 14 and 20, in float64.
25. **Coverage files:** only with `--write_coverage` (DL-015). Compat: a memory warning is printed to stderr; coverage is built one reference and strand at a time, adding in read order (DL-047). Positions past the end of the locus's reference are skipped; they occur only when records on different references merge (rule 4, #24) (DL-048).

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
| — | Multimapper shares are computed from coverage before non-G loci are dropped; if all loci with coverage are dropped, the read is lost silently (G2) | fix + compat (DL-039) | reads with a usable locus that had no coverage are kept instead of lost (toy: 1 of 4); candidates: TG3 2,710 / 100,264 multimapped reads, test flower_1 196 / 2,271 |
| — | Reads clipped at both ends are discarded if the 3′ clip has any non-G base, even when the 5′ clip is all G | fix + compat (DL-040) | TG3: 4.46 % of unique reads kept (as cap-positive); test flower_1: 0.19 % |
| #33 | Untemplated G's are added to the cap track before the length check (l.449–453 vs l.460): reads too short to be counted still add cap signal | fix + compat | test flower_1: 2,255 too-short unique reads with an all-G 5′ clip add cap signal only in the original; honeybee TG3: 787 |
| #34 | A read clipped at both ends with G's on both sides passes the check, and the 3′ G's go into the cap track at the read's 3′ end (l.449–453) | fix + compat | TG3: 1,457 unique reads; test flower_1: 5 |
| #35 | Taking back clipped bases: case-sensitive (lowercase genome never matches), `N` matches `N`, both-end clips walked from one outer end, bases past the chromosome edge ignored | fix; **not reproduced in compat** (DL-041) | none on STAR output: 0 bases taken back in TG3 (430,097 reads) and test flower_1 (32,340) |
| #36 | BODY multimapper blocks keyed by CIGAR (no sequence kept without `--untemp_out`); all reads of a block take the last read's loci | fix + compat (DL-042) | test set BODY: 21–26 % of multimapped reads counted at loci they aligned to |
| #37 | Last bedgraph run dropped when a chromosome/strand's only value is at position 0 | fix + compat (DL-044) | toy only; one extra row per such chromosome/strand |
| #38 | Cap track compared as text: `0.0` rows written, `1` / `1.0` neighbours not merged, integer repr for unique-only positions | fix + compat (DL-044) | test set: one `0.0` row each in flower_2, flower_3 `uG.minus` |
| #39 | BODY coverage = sum of strand values already rounded to 2 dp | fix + compat (DL-044) | toy 0.66 vs 0.6667; file has no consumer |
| #40 | PE: 5′ end over both mates; mate 2 only with 0x2; non-G clip on mate 1 drops only mate 1 | fix; compat refuses paired-end input (DL-043, DL-046) | none with published STAR flags (PE `EndToEnd`, no protrusion) |
| #25 | Splice-junction tally uses a stale chromosome | not ported (option unused, DL-015) | — |
| #9 | PE 5P aligned end-to-end | pipeline (STAR args), not here | — |

## Parameters
| Flag | Default | Original equivalent |
|---|---|---|
| `--library` | — (required) | `-R 5P` / `-R BODY` plus `--softclip_type 5p --untemp_out G` (5P), `--allow_nonstranded` (BODY) |
| `--min_aligned` | 16 | `--minmatch 16` |
| `--multimapper_mode` | `rescue` | always rescue, `--allow_naive` |
| `--decimals` | 4 | 2 (hard-coded); not allowed with `--nanopare_compat` (DL-047) |
| `--write_coverage` | off | coverage always written |
| `--nanopare_compat` | off | single-end only: paired-end input is an error (DL-046) |

Not ported: `--secondary` (always on), `--weighted` (automatic when `XW` present), size-class, BED, junction and multimapper-dump outputs (DL-015).

## Stats file
`<prefix>.end-counts.stats.json`, `schema_version: 1`.
- **Parameters:** `library`, `min_aligned`, `multimapper_mode`, `decimals`, `nanopare_compat`.
- **Read fates** (they sum to `reads_total`):
  - `unique`, `unique_after_filter`;
  - `multimapped_rescued`, `multimapped_naive`, `multimapped_uniform`, `multimapped_dropped`;
  - `multimapped_lost_compat` (compat only, DL-039): reads the original loses because all its loci with coverage were non-G (G2);
  - `discarded_too_short`, `discarded_untemplated_non_g`, `discarded_3prime_clip_only` (compat);
  - `unmapped` (reads whose records are all unmapped);
  - `discarded_flagged` (approved 21/10, DL-037: reads whose records are all removed by rule 2 but not all unmapped, e.g. duplicate- or QC-fail-flagged);
  - `discarded_no_mate1` (approved 21/10, DL-037: PE reads whose only loci lack an aligned mate 1, rule 17).
- **Loci:** `loci_too_short`, `loci_untemplated_non_g`, `loci_no_mate1` (approved 21/10), `loci_key_collisions` (compat).
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
  - both-end clip: `GG` 5′ clip + `AA` 3′ clip → kept, cap-positive (compat: discarded);
  - G2 toy cases: both loci non-G → discarded; usable locus without coverage + non-G locus with coverage → unique (compat: lost); no coverage anywhere → usable locus;
  - every edge case above.
- **Golden (compat):**
  - a collaborator's honeybee 5P/BODY run (STAR BAMs and EndMap bedgraphs of the original); a patched original reproduced them byte-for-byte on one sample (13/10). Fixtures: an 8 MB subsample.
  - nanoPARE test set after Phase 5b: `tests/fixtures/end_counts_golden/` (flower_1 subsample, DL-038), `tests/test_end_counts_golden.py`. **Exact match 29/10** (all 7 files).
  - Expected: exact match in compat mode.
- **Default vs golden:** the difference must be explained entirely by the quirk table. **Done 29/10** on the golden fixture: fully attributed, no residual (`nanopare_nfcore_design/package_build/29_10_2026_end-counts_default_vs_golden.md`).

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
- **D4: decided 22/10 (DL-041):** the four take-back differences (case, `N`, walk order, chromosome edge) are not reproduced in compat; unreachable with STAR output (G3).
- **D5: decided 22/10 (DL-040):** only the 5′ clip is examined; reads clipped at both ends are kept (capped if the 5′ clip is all G); 3′ G's are never cap signal. Compat reproduces the original. To confirm with the authors (Q8).
- **D6: decided 23/10 (DL-042):** BODY multimapper blocks keyed by CIGAR is a bug of the original (known issue 36). Default unchanged (rule 20); compat reproduces it for BODY.
- **D7: decided 23/10 (DL-043), changed 30/10 (DL-046):** compat refuses paired-end input. Default mode handles it by rule 17; the original's PE differences (G6, known issue 40) are not reproduced. Paired-end behaviour is with the authors (Q10–Q12).
- **Gaps (G):** behaviour of the original needed only for `--nanopare_compat` and the golden test. Each needs a targeted audit (nanopare-port-audit) of `readmapIO.py` before compat is implemented. Default mode doesn't depend on them.
  - **G1: closed 17/10 [V]** by a toy-SAM run of the original (`answers/17_10_2026_multimapper_rescue_audit.md` §6): sum of same-strand coverage over the footprint; ascending first pass, alternating; immediate update; zero-coverage blocks deferred; chromosome-blind integer `minpos`; block uses its last read's loci; ties follow first-seen order on Python ≥ 3.6 and the hash seed on Python 3.5 (pinned in `OPTIONS_lmod`), so golden runs need Python ≥ 3.6 or a fixed `PYTHONHASHSEED`. Not run, read from code only: tier = locus count before the length filter (l.284 vs l.593) [I].
  - **G2: closed 22/10 [V]** (toy SAM, `nanopare_nfcore_design/answers/22_10_2026_end-counts_g2_g3.md`). The original computes each locus's share from coverage **before** it drops loci with a non-G untemplated 5′ base (l.606–633, then l.763–766). If every locus with coverage is dropped, the read is lost (`return`, l.770–772) and never retried, even when a usable locus remains; the run log still counts it as rescued. A read whose usable loci all lack coverage is deferred to the equal-split pass and lands on its usable locus. Default mode drops unusable loci first (rule 8), so a read with one usable locus is a unique read. Compat: reproduce the loss and count it as `multimapped_lost_compat` (decided 22/10, DL-039).
  - **G3: closed 22/10 [V]** (toy SAMs + both STAR BAMs, `nanopare_nfcore_design/answers/22_10_2026_end-counts_g2_g3.md`). How the original takes back clipped bases (l.382–419):
    - it compares case-sensitively (the genome is loaded as is, `fasta_utils.import_genome`), so lowercase genome bases never match; `N` in the read matches `N` in the genome;
    - for a read clipped at both ends, both clips are walked in descending position order, so one of them is walked from its outer end and is never taken back beyond its innermost base;
    - clipped bases beyond the chromosome start/end are skipped (neither taken back nor untemplated); a 5′ end at chromosome position 1 is then missing from the ends file, though coverage is written.
    - **None of this is reachable with STAR output:** on TG3 5P (430,097 primary records with a 5′ clip) and test flower_1 (32,340), no clipped base adjacent to the alignment ever matches the genome, so nothing is ever taken back. Not reproduced in compat (decided 22/10, DL-041).
  - **Reads clipped at both ends (found with G3) [V]:** the original's untemplated check (l.444–447) includes the 3′ clip, so a read with an all-G 5′ clip and any non-G 3′ base is discarded. Default keeps it (rules 8–9): TG3 4.46 % of unique reads, test flower_1 0.19 %. Compat: discard (quirk table). Decided 22/10: keep (DL-040).
  - **G4: closed 23/10 [V]** (toy SAMs + test-set outputs, `nanopare_nfcore_design/answers/23_10_2026_end-counts_g4_g7.md` §1). Ends and coverage: each value rounded to 2 dp **before** runs are merged (so `1.0028` merges with `1.0`); rows that round to 0 are not written; float repr (`1.0`, `0.33`). Rows are already in bytewise chromosome order, then start; `bedtools sort` changes nothing. Cap track: compared as text, so `0.0` rows **are** written (test set: one each in flower_2 and flower_3 `uG.minus`) and `1` / `1.0` side by side are not merged; integer repr (`1`) where only unique reads contributed; the `#chrom start end G` header is dropped by `bedtools sort`. On every track the last run is dropped when a chromosome/strand's only value is at position 0 (`if position`, l.1015).
  - **G5: closed 23/10 [V]** (`nanopare_nfcore_design/answers/23_10_2026_end-counts_g4_g7.md` §2). BODY soft-clipped bases are ignored (no `--untemp_out`, no `--map_softclip`): coverage and 5′ end come from aligned bases only. Same as default; no compat work.
  - **G6: closed 23/10 [V]** (toy SAM, `nanopare_nfcore_design/answers/23_10_2026_end-counts_g4_g7.md` §3). Mate 2 is used only with flag 0x2; a read with only mate 2 aligned is lost silently; the 5′ end is the leftmost/rightmost position over **both** mates; a non-G clip on mate 1 of a unique pair drops only mate 1, so the 5′ end moves to mate 2. The last two differ from rule 17 but are unreachable with the published flags (PE `EndToEnd`, no `--alignEndsProtrude`). Not reproduced in compat (D7, DL-043).
  - **G7: closed 23/10 [V]** (`nanopare_nfcore_design/answers/23_10_2026_end-counts_g4_g7.md` §4). `bedgraph_combine.py` sums the plus and minus coverage per nt from the files already rounded to 2 dp (toy: `0.33 + 0.33 = 0.66`, default `0.6667`), then writes with the same writer as G4. No later stage reads `BODY.cov.bedgraph`.
  - **BODY multimapper blocks (found with G7) [V]:** without `--untemp_out` the block identifier is the CIGAR of the leftmost alignment, not the sequence (l.285–307). BODY reads with different sequences but the same locus count, leftmost start and CIGAR form one block and all take the last read's loci. Test set BODY: 1,021–1,548 reads per sample (21–26 % of multimapped reads) moved to loci they never aligned to. Reproduced in compat; default unchanged (D6, DL-042).
  - **G9: closed 29/10 [V]** (read from code, `readmapIO.py` l.556–632, l.898–958; DL-045). Details of the original's multimapper passes needed for compat, now in rule 20: the length filter (CIGAR positions, ≥) runs when a block is shared, not when it is read, so a tier-≥ 2 read with one long-enough locus is still shared as a multimapper; a block with no long-enough locus is dropped silently; a rescue pass counts as successful only when every block of at least one leftmost start was placed, so a pass that places some blocks but leaves one at every start ends the loop; the final pass is one ascending pass in which a block is shared by coverage if any exists (including coverage added earlier in that pass) and equally otherwise.
  - **G8: closed 18/10 [V]** (`answers/18_10_2026_untemplated_non_g_audit.md`). The published original discards non-G untemplated reads (`--untemp_out G` acts as a filter, `readmapIO.py` l.444–447); rule 8 matches it. The collaborator's honeybee results differ because their `endMap.sh` uses `--untemp_out A C G N T` (and STAR `Local`); the unchanged original with those flags reproduces their TG3 5P files byte-for-byte. Golden data for compat: the original re-run with the published flags on their BAMs.
- **Implementation notes (not rules):**
  - ends and cap are kept sparse;
  - coverage is needed only where multimapper loci fall, so it can be restricted to those positions after a first pass;
  - multimapper records are held in memory until uniques are done (a large fraction of 5P reads in typical runs).
