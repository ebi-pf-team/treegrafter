# Annotation File Option (`-a`) for treegrafter.py

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow users to specify the PAINT annotation file path independently from the `datadir` positional argument, for both `prepare` and `run` subcommands.

**Context:** This mirrors the Perl TreeGrafter implementation from [fullgo_paint_update#77](https://github.com/pantherdb/fullgo_paint_update/issues/77). In the Python version, the annotation file is currently hard-coded as `<datadir>/PAINT_Annotations/PAINT_Annotations_TOTAL.txt` in `prepare()` and per-family JSON files are read from `<datadir>/PAINT_Annotations/` in `process_tree()`.

**Architecture:** Add a new `-a` CLI option to both `prepare` and `run` subcommands. When provided, it overrides the default annotation file path. When omitted, behavior is unchanged (backwards-compatible).

**Tech Stack:** Python 3, argparse

---

## File Structure

- **Modify:** `treegrafter.py` — Add `-a` option to both subcommand parsers, thread through relevant functions.

---

### Task 1: Add `-a` option to `prepare` subcommand

The `prepare` subcommand reads `PAINT_Annotations_TOTAL.txt` and writes per-family `.json` files to the same `PAINT_Annotations/` directory. When `-a` is provided, the annotation source file changes but the JSON output directory must also be configurable (or remain in `datadir`).

**Key decision:** The `-a` flag only overrides the *input* annotation file. The per-family JSON output still goes to `<datadir>/PAINT_Annotations/` since that's where `run` reads them from.

**Files:**
- Modify: `treegrafter.py:621-623` (prepare subparser)
- Modify: `treegrafter.py:501-506` (prepare function)

- [ ] **Step 1: Add `-a` argument to `prepare` subparser**

At `treegrafter.py:622`, after the `datadir` argument, add:

```python
parser_pre.add_argument("-a", dest="annotation_file", metavar="FILE",
                        help="PAINT annotation file path "
                             "(default: <datadir>/PAINT_Annotations/"
                             "PAINT_Annotations_TOTAL.txt)")
```

- [ ] **Step 2: Use `-a` in `prepare()` with fallback to default**

At `treegrafter.py:501-506`, update the `prepare` function to use `args.annotation_file` when provided:

```python
def prepare(args):
    datadir = args.datadir

    sys.stderr.write("Loading PAINT annotations\n")
    paintdir = os.path.join(datadir, "PAINT_Annotations")

    if args.annotation_file:
        paintfile = args.annotation_file
    else:
        paintfile = os.path.join(paintdir, "PAINT_Annotations_TOTAL.txt")
```

No other changes needed in `prepare()` — the JSON output still writes to `paintdir` (line 540).

---

### Task 2: Add `-a` option to `run` subcommand

The `run` subcommand doesn't read the TOTAL annotation file directly — it reads per-family `.json` files from `<datadir>/PAINT_Annotations/`. The `-a` flag for `run` needs to override the *directory* where these JSON files are read from.

**Key decision:** For `run`, the `-a` flag specifies an alternative `PAINT_Annotations/` directory (containing the per-family `.json` files), not the TOTAL text file. This is the directory that `prepare` wrote the JSON files into.

Actually, re-reading the Perl plan: the `-a` flag in Perl points to the annotation *file* (`PAINT_Annotations_TOTAL.txt`), and the Perl code reads it at runtime. The Python version pre-processes this into per-family JSON files during `prepare`. So for `run`, the equivalent is to override the directory containing those JSON files.

**Revised decision:** Add `-A` (or `--annot-dir`) to `run` to specify an alternative directory containing the per-family annotation JSON files. This is the Python-specific equivalent of the Perl `-a` flag for the `run` phase. Alternatively, we can keep `-a` and have it point to the directory containing the `.json` files.

**Simplest approach:** Use `-a` on `run` to specify an alternative `PAINT_Annotations` directory (the one containing `*.json` files). This directory would be the parent dir of the JSON files that `prepare` wrote.

**Files:**
- Modify: `treegrafter.py:625-645` (run subparser)
- Modify: `treegrafter.py:562-601` (run function)
- Modify: `treegrafter.py:28` (process_matches_epang signature)
- Modify: `treegrafter.py:178` (process_tree signature)
- Modify: `treegrafter.py:216` (annotation file path in process_tree)
- Modify: `treegrafter.py:59-65` (process_tree call in process_matches_epang)

- [ ] **Step 3: Add `-a` argument to `run` subparser**

At `treegrafter.py:644` (before `set_defaults`), add:

```python
parser_run.add_argument("-a", dest="annot_dir", metavar="DIR",
                        help="directory containing per-family annotation "
                             "JSON files from 'prepare' step "
                             "(default: <datadir>/PAINT_Annotations)")
```

- [ ] **Step 4: Resolve `annot_dir` in `run()` with fallback**

At `treegrafter.py:562`, after validation, compute the effective annotation directory:

```python
annot_dir = args.annot_dir or os.path.join(args.datadir, "PAINT_Annotations")
```

Add validation that `annot_dir` exists:

```python
if not os.path.isdir(annot_dir):
    sys.stderr.write("Error: {}: "
                     "no such directory.\n".format(annot_dir))
    sys.exit(1)
```

- [ ] **Step 5: Thread `annot_dir` through `process_matches_epang()`**

Update the signature at line 28:

```python
def process_matches_epang(matches, datadir, tempdir, binary=None, threads=1, print_go=False, annot_dir=None):
```

Update the `process_tree` call at line 65:

```python
for result in process_tree(pthr, result_tree, matches[pthr], datadir, print_go=print_go, annot_dir=annot_dir):
```

Update the call in `run()` at line 598:

```python
results = process_matches_epang(matches, args.datadir, tempdir,
                                binary=args.epang,
                                threads=args.threads,
                                print_go=args.print_go,
                                annot_dir=annot_dir)
```

- [ ] **Step 6: Thread `annot_dir` through `process_tree()`**

Update the signature at line 178:

```python
def process_tree(pthr, result_tree, pthr_matches, datadir, print_go=False, annot_dir=None):
```

Update the annotation file path at line 216:

```python
effective_annot_dir = annot_dir or os.path.join(datadir, 'PAINT_Annotations')
annot_file = os.path.join(effective_annot_dir, pthr + '.json')
```

---

### Task 3: Test

- [ ] **Step 7: Test `prepare` with default behavior (no `-a`)**

```bash
python treegrafter.py prepare ./Test/PANTHER_mini
```

Expected: Works as before — reads from default annotation file location.

- [ ] **Step 8: Test `prepare` with explicit `-a` pointing to same file**

```bash
python treegrafter.py prepare -a ./Test/PANTHER_mini/PAINT_Annotations/PAINT_Annotations_TOTAL.txt ./Test/PANTHER_mini
```

Expected: Same result as default — backwards compatible.

- [ ] **Step 9: Test `run` with default behavior (no `-a`)**

```bash
python treegrafter.py run -o /dev/stdout --print-go ./Test/sample.fasta ./Test/sample.fasta.hmmsearch.out ./Test/PANTHER_mini
```

Expected: Works as before.

- [ ] **Step 10: Test `run` with explicit `-a` pointing to same dir**

```bash
python treegrafter.py run -a ./Test/PANTHER_mini/PAINT_Annotations -o /dev/stdout --print-go ./Test/sample.fasta ./Test/sample.fasta.hmmsearch.out ./Test/PANTHER_mini
```

Expected: Same result as default.

- [ ] **Step 11: Test `prepare` with nonexistent `-a` path**

```bash
python treegrafter.py prepare -a /nonexistent/file.txt ./Test/PANTHER_mini 2>&1
```

Expected: Python raises `FileNotFoundError` naturally when trying to open the file. Optionally add explicit validation with a user-friendly error message.

---

## Summary of Changes

| Location | Change |
|----------|--------|
| `main()` prepare subparser | Add `-a` / `--annotation-file` argument |
| `main()` run subparser | Add `-a` / `--annot-dir` argument |
| `prepare()` | Use `args.annotation_file` with fallback to default path |
| `run()` | Resolve `annot_dir`, validate, pass to `process_matches_epang()` |
| `process_matches_epang()` | Accept and forward `annot_dir` parameter |
| `process_tree()` | Accept `annot_dir` parameter, use it for JSON file lookup |

All changes are backwards-compatible — omitting `-a` preserves existing behavior.