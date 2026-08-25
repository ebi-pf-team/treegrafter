# Add `--print-go` Flag for GO/PC Annotations

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `--print-go` flag to the `run` subcommand that appends GO terms and protein class columns to the TSV output, sourced from the per-family PAINT annotation JSONs.

**Architecture:** The PAINT annotation data is already loaded into per-family JSON files during `prepare`. Each JSON maps `node_id` to `[subfam_id, go_terms, protein_class, graft_point]`. The `process_tree` function (line 214-215) already opens these JSONs but only extracts `pthrsf` (subfam_id), discarding GO/PC. We thread a `print_go` flag through to `process_tree` and `process_matches_epang`, and when enabled, include the GO and PC fields in the result rows and output header.

**Tech Stack:** Python 3, no new dependencies.

---

## Data Flow

Current result row (13 columns):
```
[query_id, panther_id, score, evalue, dom_score, dom_evalue,
 hmm_start, hmm_end, ali_start, ali_end, env_start, env_end, node_id]
```

With `--print-go` (15 columns):
```
[query_id, panther_id, score, evalue, dom_score, dom_evalue,
 hmm_start, hmm_end, ali_start, ali_end, env_start, env_end, node_id,
 go_terms, protein_class]
```

- `go_terms`: comma-separated GO IDs (e.g. `GO:0005634,GO:0006355`) or `-` if none
- `protein_class`: PC ID (e.g. `PC00216`) or `-` if none

---

## File Structure

All changes are in a single file:
- **Modify:** `treegrafter.py` — add flag, thread it through, update output

---

## Task 1: Add `--print-go` argument to the `run` subparser

**Files:**
- Modify: `treegrafter.py:630-632`

- [ ] **Step 1: Add the argument**

After the `--keep` argument (line 631), add:

```python
    parser_run.add_argument("--print-go", action="store_true",
                            help="include GO terms and protein class in output")
```

---

## Task 2: Update `process_tree` to optionally return GO/PC annotations

**Files:**
- Modify: `treegrafter.py:153-233`

- [ ] **Step 1: Add `print_go` parameter to `process_tree`**

Change the function signature at line 153 from:

```python
def process_tree(pthr, result_tree, pthr_matches, datadir):
```

to:

```python
def process_tree(pthr, result_tree, pthr_matches, datadir, print_go=False):
```

- [ ] **Step 2: Unpack GO/PC from the JSON**

At line 215, change:

```python
            pthrsf, _, _, _ = json.load(fh)[common_an]
```

to:

```python
            pthrsf, go_terms, protein_class, _ = json.load(fh)[common_an]
```

- [ ] **Step 3: Conditionally append GO/PC to the result row**

At line 231 (after `common_an` in the result list), add two extra columns when `print_go` is enabled. Change the result construction at lines 217-231 from:

```python
        results_pthr.append([
            query_id,
            pthrsf or pthr,
            pthr_matches[query_id]['score'][0],
            pthr_matches[query_id]['evalue'][0],
            pthr_matches[query_id]['domscore'][0],
            pthr_matches[query_id]['domevalue'][0],
            pthr_matches[query_id]['hmmstart'][0],
            pthr_matches[query_id]['hmmend'][0],
            pthr_matches[query_id]['alifrom'][0],
            pthr_matches[query_id]['alito'][0],
            pthr_matches[query_id]['envfrom'][0],
            pthr_matches[query_id]['envto'][0],
            common_an
        ])
```

to:

```python
        row = [
            query_id,
            pthrsf or pthr,
            pthr_matches[query_id]['score'][0],
            pthr_matches[query_id]['evalue'][0],
            pthr_matches[query_id]['domscore'][0],
            pthr_matches[query_id]['domevalue'][0],
            pthr_matches[query_id]['hmmstart'][0],
            pthr_matches[query_id]['hmmend'][0],
            pthr_matches[query_id]['alifrom'][0],
            pthr_matches[query_id]['alito'][0],
            pthr_matches[query_id]['envfrom'][0],
            pthr_matches[query_id]['envto'][0],
            common_an
        ]
        if print_go:
            row.append(go_terms or '-')
            row.append(protein_class or '-')
        results_pthr.append(row)
```

---

## Task 3: Update `process_matches_epang` to thread `print_go` through

**Files:**
- Modify: `treegrafter.py:28-66`

- [ ] **Step 1: Add `print_go` parameter**

Change the function signature at line 28 from:

```python
def process_matches_epang(matches, datadir, tempdir, binary=None, threads=1):
```

to:

```python
def process_matches_epang(matches, datadir, tempdir, binary=None, threads=1, print_go=False):
```

- [ ] **Step 2: Pass `print_go` to `process_tree`**

At line 62, change:

```python
        for result in process_tree(pthr, result_tree, matches[pthr], datadir):
```

to:

```python
        for result in process_tree(pthr, result_tree, matches[pthr], datadir, print_go=print_go):
```

---

## Task 4: Update `run` to use the flag

**Files:**
- Modify: `treegrafter.py:555-600`

- [ ] **Step 1: Update the header line**

At lines 583-585, change:

```python
    fh.write("query_id\tpanther_id\tscore\tevalue\tdom_score\tdom_evalue\t"
             "hmm_start\thmm_end\tali_start\tali_end\tenv_start\tenv_end\t"
             "node_id\n")
```

to:

```python
    header = ("query_id\tpanther_id\tscore\tevalue\tdom_score\tdom_evalue\t"
              "hmm_start\thmm_end\tali_start\tali_end\tenv_start\tenv_end\t"
              "node_id")
    if args.print_go:
        header += "\tgo_terms\tprotein_class"
    fh.write(header + "\n")
```

- [ ] **Step 2: Pass `print_go` to `process_matches_epang`**

At lines 588-590, change:

```python
        results = process_matches_epang(matches, args.datadir, tempdir,
                                        binary=args.epang,
                                        threads=args.threads)
```

to:

```python
        results = process_matches_epang(matches, args.datadir, tempdir,
                                        binary=args.epang,
                                        threads=args.threads,
                                        print_go=args.print_go)
```

---

## Task 5: Handle the fallback row (no EPA-ng result)

**Files:**
- Modify: `treegrafter.py:39-53`

When EPA-ng fails or a query has no tree placement, the fallback row at lines 40-54 has 13 columns with `"-"` as node_id. If `print_go` is enabled, this row also needs the extra columns (both `-` since there's no node to look up annotations from).

- [ ] **Step 1: Add `print_go` to the fallback row logic**

Change lines 39-54 from:

```python
        for query_id in matches[pthr]:
            results[query_id] = [
                query_id,
                pthr,
                matches[pthr][query_id]['score'][0],
                matches[pthr][query_id]['evalue'][0],
                matches[pthr][query_id]['domscore'][0],
                matches[pthr][query_id]['domevalue'][0],
                matches[pthr][query_id]['hmmstart'][0],
                matches[pthr][query_id]['hmmend'][0],
                matches[pthr][query_id]['alifrom'][0],
                matches[pthr][query_id]['alito'][0],
                matches[pthr][query_id]['envfrom'][0],
                matches[pthr][query_id]['envto'][0],
                "-"
            ]
```

to:

```python
        for query_id in matches[pthr]:
            row = [
                query_id,
                pthr,
                matches[pthr][query_id]['score'][0],
                matches[pthr][query_id]['evalue'][0],
                matches[pthr][query_id]['domscore'][0],
                matches[pthr][query_id]['domevalue'][0],
                matches[pthr][query_id]['hmmstart'][0],
                matches[pthr][query_id]['hmmend'][0],
                matches[pthr][query_id]['alifrom'][0],
                matches[pthr][query_id]['alito'][0],
                matches[pthr][query_id]['envfrom'][0],
                matches[pthr][query_id]['envto'][0],
                "-"
            ]
            if print_go:
                row.extend(['-', '-'])
            results[query_id] = row
```

---

## Task 6: Integration test

- [ ] **Step 1: Run treegrafter without `--print-go`**

Run the same command as before. Verify output is unchanged (13 columns, same header).

- [ ] **Step 2: Run treegrafter with `--print-go`**

Run with `--print-go` added. Verify:
- Header has 15 columns: the original 13 plus `go_terms` and `protein_class`
- Data rows have 15 tab-separated fields
- GO terms appear as comma-separated `GO:NNNNNNN` values or `-`
- Protein class appears as `PCNNNNN` or `-`

---

## Notes

- The `--print-go` flag is intentionally off by default to preserve backward compatibility with any downstream parsers that expect exactly 13 columns.
- The flag name `--print-go` covers both GO terms and protein class since they come from the same annotation source and are typically wanted together.
- No changes to `prepare` are needed — the per-family JSONs already contain all the data.