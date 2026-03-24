# Design: Separate Annotation Data from Library Data

**Date:** 2026-03-19
**Issue:** [pantherdb/fullgo_paint_update#77](https://github.com/pantherdb/fullgo_paint_update/issues/77)

## Problem

The current `prepare` command mixes two unrelated operations — annotation processing and sequence normalization — and writes annotation output into the PANTHER library directory (`datadir`). Since annotations update monthly and the library updates yearly, this couples fast-moving data to a slow-moving directory, making it awkward to manage multiple annotation versions or keep the library directory read-only.

## Design Decisions

1. **Remove sequence normalization** from `prepare`. The U/O→X fasta rewriting belongs upstream in the PANTHER library build, not in TreeGrafter.
2. **`prepare` becomes annotation-only.** It takes an input annotation file and an output directory for per-family JSONs.
3. **`run` takes two separate directory arguments** — one for the library (trees, fastas, HMMs) and one for the annotation JSONs.
4. **Use named flags `-d` and `-a`** on `run` instead of positional arguments, so the two directories are self-documenting and order-independent.
5. **Clean break, no backwards compatibility shim.** The data layout is fundamentally changing; old invocations will get an argparse error.

## CLI Interface

### `prepare`

```
treegrafter.py prepare ANNOTATION_FILE OUTPUT_DIR
```

- `ANNOTATION_FILE` — path to `PAINT_Annotations_TOTAL.txt` (or equivalent)
- `OUTPUT_DIR` — directory where per-family `{family_id}.json` files are written; created if it doesn't exist

### `run`

```
treegrafter.py run FASTA HMMSEARCH -d LIBDIR -a ANNOTDIR [options]
```

- `FASTA` — query sequences (positional, unchanged)
- `HMMSEARCH` — hmmsearch output file (positional, unchanged)
- `-d LIBDIR` — **required** — PANTHER library directory containing `Tree_MSF/` and `famhmm/`
- `-a ANNOTDIR` — **required** — directory containing per-family annotation JSON files (output of `prepare`)
- Other flags unchanged: `-e`, `-o`, `--epa-ng`, `-t`, `-T`, `--keep`, `--print-go`

### `treegrafter.sh search`

```
treegrafter.sh search FASTA LIBDIR ANNOTDIR OUTPUT
```

- Runs `hmmsearch` using `LIBDIR/famhmm/binHmm`
- Passes `-d LIBDIR -a ANNOTDIR` to `treegrafter.py run`

## Data Flow

```
PAINT_Annotations_TOTAL.txt
         |
         v
  treegrafter.py prepare
         |
         v
  ANNOTDIR/              (per-family JSONs, updated monthly)
    PTHR10000.json
    PTHR10001.json
    ...

  LIBDIR/                (PANTHER library, updated yearly, read-only)
    Tree_MSF/
      *.fasta
      *.newick
    famhmm/
      binHmm

  query.fasta + hmmsearch.out
         |
         v
  treegrafter.py run -d LIBDIR -a ANNOTDIR
         |
         v
  predictions.tsv
```

## Changes by File

| File | Change |
|------|--------|
| `treegrafter.py` `main()` | `prepare` parser: replace `datadir` + `-a` with positional `annotation_file` + `output_dir` |
| `treegrafter.py` `main()` | `run` parser: remove positional `datadir`, add required `-d` (`libdir`) and `-a` (`annotdir`) |
| `treegrafter.py` `prepare()` | Remove sequence normalization (lines 547-563). Read from `args.annotation_file`, write JSONs to `args.output_dir` (create dir if needed). |
| `treegrafter.py` `run()` | Use `args.libdir` for tree/fasta/hmm paths, `args.annotdir` for JSON paths. Drop `annot_dir` override logic. |
| `treegrafter.py` `process_matches_epang()` | Split `datadir` param into `libdir` + `annot_dir` (both required, no fallback) |
| `treegrafter.py` `process_tree()` | Split `datadir` param into `libdir` + `annot_dir` (both required, no fallback) |
| `treegrafter.py` `_commonancestor()` | `datadir` → `libdir` |
| `treegrafter.py` `generate_fasta_for_panthr()` | `datadir` → `libdir` |
| `treegrafter.py` `_run_epang()` | `datadir` → `libdir` |
| `treegrafter.py` `align_length()` | `datadir` → `libdir` |
| `treegrafter.sh` | `search` signature: `FASTA LIBDIR ANNOTDIR OUTPUT`. Pass `-d`/`-a` to `run`. |

## Error Handling

- `prepare`: validate `ANNOTATION_FILE` exists (fail with message if not), create `OUTPUT_DIR` if it doesn't exist
- `run`: validate `FASTA` and `HMMSEARCH` are files, validate `-d` and `-a` are directories

## Breaking Changes

- `prepare DATADIR` → `prepare ANNOTATION_FILE OUTPUT_DIR`
- `run ... DATADIR` → `run ... -d LIBDIR -a ANNOTDIR`
- `treegrafter.sh search FASTA DATADIR OUTPUT` → `treegrafter.sh search FASTA LIBDIR ANNOTDIR OUTPUT`
