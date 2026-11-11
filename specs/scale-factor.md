# `scale-factor`: behavioural spec

Status: **approved** (11/11/2026, DL-063; D1–D4 decided 10/11, DL-058–DL-061) · Replaces (nanoPARE 6264d6f): `scripts/python_scripts/gtf_quantify.py` (×2) + `endgraph_ratio.py`, called at `endGraph.sh` l.212–232; cap at l.243–247
Pipeline module: `modules/local/scale_factor/` · Audit: `nanopare_nfcore_design/answers/10_11_2026_scale_factor_audit.md`

## Purpose
EndGraph subtracts the BODY (RNA-seq) signal from the 5P signal before smoothing and peak calling, so that 5′ ends explained by ordinary fragmentation of the transcript body are removed. The two libraries have different depths, and one RNA molecule gives one 5P read but many BODY reads. This subcommand computes one number per sample, the factor that 5P values are multiplied by before BODY is subtracted (`factor × 5P − BODY`). Consumer: `kde-peaks`.

## Inputs
| Argument | Format | Required | Notes |
|---|---|---|---|
| `--ends_plus`, `--ends_minus` | bedgraph, 5P ends after the TSO mask | yes | per-nt rows (end-counts output); 0-based half-open |
| `--body_plus`, `--body_minus` | bedgraph, BODY ends | yes | as written by `end-counts --library BODY`; strand as aligned |
| `--regions` | BED6: chrom, start, end, gene_id, score (ignored), strand | yes | gene regions = exons of all the gene's transcripts ∪ `buffer` nt before the first and after the last exon of each transcript, merged per gene and strand. Built by the pipeline (gffread GTF + awk + bedtools merge). Several rows per gene allowed; rows of one gene must not overlap each other |

A bedgraph row spanning n nt contributes its value at each of its n positions (merged and per-nt bedgraphs give the same result). Several rows may cover the same position; their values add.

## Outputs
| File | Format | Consumer |
|---|---|---|
| `<prefix>.scale_factor.txt` | one line: the factor after the cap, 4 decimals (`0.2121`) | `kde-peaks` |
| `<prefix>.scale-factor.stats.json` | see Stats file | MultiQC |
| `<prefix>.scale-factor.genes.tsv` | only with `--write_genes`: header `gene_id length end body`; one row per gene in `--regions`, sorted by gene_id (bytewise); end and body split by rule 4, 4 fixed decimals | none (QC) |

## Rules
1. **Sense strand.** A gene on `+` is matched with `ends_plus`. With `--body_strand reverse` (default) it is matched with `body_minus`, and a `−` gene with `body_plus`; with `forward`, with `body_plus` / `body_minus` respectively.
2. **Ownership.** A position (chrom, strand, nt) is owned by every gene whose region covers it on that strand.
3. **Totals.** E = sum of 5P values at owned positions; B = sum of sense BODY values at owned positions. Each position counts once, however many genes own it. Values outside every region are ignored. Sums are exact (`math.fsum`), so they do not depend on row order.
4. **Per-gene BODY counts.** b(g) = sum over g's positions of the BODY value divided by the number of owners of that position.
5. **Gene length** len(g) = number of positions in g's region (DL-059).
6. **Mean length.** L̄ = Σ b(g) / Σ (b(g) / len(g)) over genes with b(g) > 0.
7. **Factor.** raw = (2 · F / L̄) · (B / E), F = `--fraglen`.
8. **Cap.** factor = min(raw, `--max_scale_factor`). Written rounded half-even to 4 decimals.
9. Output is the same for any row order of the inputs.

## Edge cases
- E = 0, B = 0, or no gene with b(g) > 0: no factor can be computed → exit non-zero with a message naming the empty quantity (DL-060).
- Bedgraph chromosome absent from `--regions`: its values are ignored, counted in stats.
- Negative or non-finite bedgraph value: error.
- Overlapping rows of the same gene in `--regions`: error (the pipeline merges them).
- Region on strand `.`: error (gene regions are stranded).

## Quirks of the original
| # (known_issues) | Original behaviour | Decision | Expected diff vs golden |
|---|---|---|---|
| 43 | last locus of a chromosome dropped | fix (no compat, DL-058) | test set: AT2G40740 included |
| 44 | residual reset inside the loop | fix by not having the step (D1) | 0.3924/0.2100/0.2121 → 0.3920/0.2001/0.2127 for this bug alone |
| 45 | per-gene reads = density × full length; totals inflated (BODY up to +22.6 %) | fix: counts, not densities (rules 3–4) | 0.3871 / 0.1903 / 0.2103 (−1.4 / −9.4 / −0.8 %) |
| 46 | gene = ID before `.`, isoforms by substring | not carried over: gene = `gene_id` from the regions BED | none on the test set |
| 47 | `flip_5p` hard-coded (BODY antisense) | parameter `--body_strand`, default `reverse` = original | none |
| — | gene length = longest isoform + 200 (inherited from the general gene-level TPM code, eac817b) | region size (DL-059) | ~1 %: longest + 200 would give 0.3904 / 0.1919 / 0.2121 |
| — | `round()` repr (`0.21`) | 4 fixed decimals (`0.2100`) | text only; value equal |

## Parameters
| Flag | Default | Original equivalent |
|---|---|---|
| `--fraglen` | 200 | `--fraglen` / `FRAGLEN` |
| `--body_strand` | `reverse` | `flip_5p="true"` (hard-coded) |
| `--max_scale_factor` | 10 | `SCALE_CAP=10` (hard-coded) |
| `--write_genes` | off | — (the original writes per-gene quant tables always) |
| `--prefix` | required | `$SAMPLE_NAME` |
| (pipeline) `params.scale_factor_buffer` | 100 | `gtf_quantify.py --buffer 100` |

## Stats file
`schema_version`, `fraglen`, `body_strand`, `max_scale_factor`, `end_total` (E), `body_total` (B), `mean_gene_length` (L̄), `genes_in_regions`, `genes_with_body`, `raw_factor`, `scale_factor`, `capped` (bool), `values_outside_regions` (5P and BODY, separately), `chromosomes_without_regions`.

## Tests
- Unit: one per rule; strand matching both ways; shared position split between 2 and 3 genes (rule 4) but counted once in E/B (rule 3); cap hit; input row order shuffled → same output; each edge case.
- Golden: none in compat (DL-058). Default checked against the audit's measured values on the test set: 0.3871 / 0.1903 / 0.2103, and E/B equal to 138,729 / 160,859 (flower_1; the audit's instrumented original saw 160,857 because it drops AT2G40740, 2 BODY reads, quirk 43), 214,506 / 121,465 (flower_2), 219,022 / 139,615 (flower_3), within 1 read (bedgraph values have 2 dp).
- Downstream golden tests (`kde-peaks`) take the golden factors 0.3924 / 0.2100 / 0.2121 as input.

## Open
- **G1:** where the factor 2 in 2F/L̄ comes from is not documented (audit §2). It appeared in 2f0dcd2 (08/08/2018), replacing a fixed baseline of 0.25; an earlier version (`RPM_TPM_normalize.py`, ed2ea74) used F/L̄ without the 2. Asked as Q15.

## Decided
- D1 → DL-058 (no compat), D2 → DL-059 (region size), D3 → DL-060 (error), D4 → DL-061 (`--write_genes`).
