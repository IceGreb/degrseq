# `<subcommand>`: behavioural spec

Status: draft | approved (DD/MM/YYYY) · Replaces (nanoPARE 6264d6f): `<scripts / shell lines>`
Pipeline module: `modules/local/<name>/` · Design-note section: `§<n>`

Write this from the audit (design note, answers/), not side by side with the source. It is
what the code is written from and what the tests check.

## Purpose
One paragraph: what the step is for, biologically and technically, and who consumes its output.

## Inputs
| Argument | Format | Required | Notes (coordinate base, sort order, strand) |
|---|---|---|---|

## Outputs
| File | Format / columns | Consumer |
|---|---|---|

## Rules
Numbered, exact, testable statements. Each rule gets at least one unit test.
1. …

## Edge cases
Empty input, ties, position 0, chromosome ends, missing IDs, … and the required behaviour for each.

## Quirks of the original
| # (known_issues) | Original behaviour | Decision: fix / reproduce / fix + `--nanopare_compat` | Expected diff vs golden |
|---|---|---|---|

## Parameters
| Flag | Default | Original equivalent |
|---|---|---|

## Stats file
Fields written to `<prefix>.<subcommand>.stats.json` (if any).

## Tests
- Unit: …
- Golden: input files, command, expected match (exact / numeric tolerance / compat only).
