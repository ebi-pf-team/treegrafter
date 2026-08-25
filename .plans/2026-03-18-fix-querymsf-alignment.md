# Fix Query MSF Alignment Reconstruction

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the "length of query MSF longer than expected PANTHER alignment length" error by fixing multi-block alignment parsing and rewriting `_querymsf` with a position-based placement strategy.

**Architecture:** Two bugs cause incorrect query MSF lengths. First, `parsehmmsearch` only captures the first ~80-character block of HMMER domain alignments, truncating long alignments. Second, `_querymsf` uses fragile incremental string-building that breaks with overlapping domains. We fix parsing to accumulate all alignment blocks, then rewrite `_querymsf` to pre-allocate a fixed-length gap array and place characters by position.

**Tech Stack:** Python 3, no new dependencies.

---

## Bug Analysis

There are two distinct bugs in `treegrafter.py`, both contributing to wrong-length MSF output:

### Bug 1: Multi-block alignment parsing (primary cause of reported error)

In `parsehmmsearch()` (line 368-400), when a `== domain N` section is encountered, the parser reads exactly **one** alignment block (3 lines: HMM line, consensus line, query line). But HMMER3 splits long alignments into multiple blocks of ~80 characters. For alignments longer than ~80 HMM positions, subsequent blocks are silently skipped, producing truncated `hmmalign` and `matchalign` strings.

This directly explains the error: "expected 157, got 82" — only the first block (~82 positions) was captured from a 157-position alignment.

### Bug 2: Fragile incremental `_querymsf` construction

`_querymsf()` (line 105-149) builds the output string by concatenation, assuming domains are non-overlapping and sequential. If domains overlap or gap-bridging arithmetic is off, the string length is wrong.

---

## File Structure

All changes are in a single file:
- **Modify:** `treegrafter.py` — fix `parsehmmsearch()` and rewrite `_querymsf()`
- **Create:** `tests/test_querymsf.py` — unit tests for the fixed functions
- **Create:** `tests/test_parsehmmsearch.py` — unit tests for multi-block parsing

---

## Task 1: Add unit tests for `_querymsf` with single-domain, single-block input

**Files:**
- Create: `tests/test_querymsf.py`

- [ ] **Step 1: Write the failing test for single-domain single-block**

```python
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from treegrafter import _querymsf


def test_single_domain_single_block():
    """Single domain, alignment fits in one block (< 80 chars)."""
    match_data = {
        'hmmstart': ['3'],
        'hmmend': ['8'],
        'hmmalign': ['abcdef'],
        'matchalign': ['ABCDEF'],
        'domscore': ['100.0'],
    }
    pthr_align_length = 10
    result = _querymsf(match_data, pthr_align_length)
    # positions: 1-2 = '--', 3-8 = 'ABCDEF', 9-10 = '--'
    assert result == '--ABCDEF--'
    assert len(result) == pthr_align_length
```

- [ ] **Step 2: Run test to verify it passes with current code (baseline)**

Run: `cd /Users/ebertdu/panther/ebi-treegrafter && python -m pytest tests/test_querymsf.py::test_single_domain_single_block -v`

Expected: PASS (this case works with current code)

---

## Task 2: Add unit tests for `_querymsf` with insert states

**Files:**
- Modify: `tests/test_querymsf.py`

- [ ] **Step 1: Write test for alignment with insert states**

```python
def test_single_domain_with_inserts():
    """Insert states (dots in hmmalign) should be skipped."""
    match_data = {
        'hmmstart': ['1'],
        'hmmend': ['5'],
        'hmmalign': ['ab.c.de'],  # 5 non-insert + 2 inserts
        'matchalign': ['AB-C-DE'],  # corresponding query chars
        'domscore': ['100.0'],
    }
    pthr_align_length = 5
    result = _querymsf(match_data, pthr_align_length)
    # inserts at positions 2,4 skipped; remaining: A,B,C,D,E
    assert result == 'ABCDE'
    assert len(result) == pthr_align_length
```

- [ ] **Step 2: Run test to verify it passes (baseline)**

Run: `python -m pytest tests/test_querymsf.py::test_single_domain_with_inserts -v`

Expected: PASS

---

## Task 3: Add tests for multi-domain cases

**Files:**
- Modify: `tests/test_querymsf.py`

- [ ] **Step 1: Write test for non-overlapping multi-domain placement**

```python
def test_non_overlapping_multi_domain():
    """Two non-overlapping domains should both be placed correctly."""
    match_data = {
        'hmmstart': ['1', '6'],
        'hmmend': ['3', '8'],
        'hmmalign': ['abc', 'def'],
        'matchalign': ['ABC', 'DEF'],
        'domscore': ['100.0', '90.0'],
    }
    result = _querymsf(match_data, 10)
    assert result == 'ABC--DEF--'
    assert len(result) == 10
```

- [ ] **Step 2: Write test for overlapping domains**

```python
def test_overlapping_domains_higher_score_wins():
    """Two domains overlap; higher-scoring domain should take priority."""
    match_data = {
        'hmmstart': ['1', '3'],
        'hmmend': ['5', '8'],
        'hmmalign': ['abcde', 'fghijk'],  # domain 2 overlaps positions 3-5
        'matchalign': ['ABCDE', 'FGHIJK'],
        'domscore': ['100.0', '50.0'],
    }
    pthr_align_length = 8
    result = _querymsf(match_data, pthr_align_length)
    # Domain 1 (score 100): positions 1-5 = ABCDE
    # Domain 2 (score 50): positions 3-8, but 3-5 overlap -> skip domain 2
    #   (entire domain skipped because it overlaps a higher-scoring one)
    # Positions 6-8 remain gaps since domain 2 is skipped entirely
    assert len(result) == pthr_align_length
    # Domain 1 fills positions 1-5, rest are gaps
    assert result == 'ABCDE---'
```

- [ ] **Step 3: Run tests to verify the overlap test fails**

Run: `python -m pytest tests/test_querymsf.py -v`

Expected: `test_non_overlapping_multi_domain` may pass, `test_overlapping_domains_higher_score_wins` FAILS (current code crashes or produces wrong length with overlapping domains)

---

## Task 4: Rewrite `_querymsf` with position-based placement

**Files:**
- Modify: `treegrafter.py:105-149`

- [ ] **Step 1: Replace `_querymsf` with position-based implementation**

Replace lines 105-149 in `treegrafter.py` with:

```python
def _querymsf(match_data, pthr_align_length):
    # Pre-allocate output as fixed-length gap array
    querymsf = ['-'] * pthr_align_length

    # Track which HMM positions are already claimed
    used_positions = set()

    # Sort domain indices by score (highest first)
    num_domains = len(match_data['matchalign'])
    domain_indices = list(range(num_domains))
    domain_indices.sort(
        key=lambda i: float(match_data['domscore'][i]), reverse=True
    )

    for i in domain_indices:
        hmmstart = int(match_data['hmmstart'][i])
        hmmend = int(match_data['hmmend'][i])

        # Check for overlap with already-claimed positions
        domain_positions = set(range(hmmstart, hmmend + 1))
        if domain_positions & used_positions:
            sys.stderr.write(
                "Warning: domain {} (score {}) overlaps higher-scoring "
                "domain, skipping.\n".format(i, match_data['domscore'][i])
            )
            continue

        # Mark positions as used
        used_positions |= domain_positions

        # Place alignment characters at correct MSA positions
        hmmalign = match_data['hmmalign'][i]
        matchalign = match_data['matchalign'][i]
        msa_pos = hmmstart - 1  # convert to 0-based

        for j in range(len(hmmalign)):
            if hmmalign[j] == '.':
                # Insert state — skip, don't advance MSA position
                continue
            querymsf[msa_pos] = matchalign[j]
            msa_pos += 1

    return ''.join(querymsf).upper()
```

- [ ] **Step 2: Run all `_querymsf` tests**

Run: `python -m pytest tests/test_querymsf.py -v`

Expected: ALL PASS

---

## Task 5: Add test for multi-block alignment parsing

**Files:**
- Create: `tests/test_parsehmmsearch.py`

- [ ] **Step 1: Create a minimal multi-block HMMER output fixture**

Create `tests/fixtures/multiblock.hmmsearch` containing a synthetic hmmsearch output where one domain alignment spans two blocks. The fixture needs to follow HMMER3 text output format.

```python
import sys
import os
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from treegrafter import parsehmmsearch

# Minimal HMMER3 output with a single query matching PTHR00001,
# one domain whose alignment spans two blocks.
MULTIBLOCK_HMMER = """\
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# query HMM file: test.hmm
# target sequence database: test.fa
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Query:       PTHR00001  [M=10]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence  Description
    ------- ------ -----    ------- ------ -----   ---- --  --------  -----------
    1.0e-10  100.0   0.0    1.0e-10  100.0   0.0    1.0  1  query1    test query

Domain annotation for each sequence (and target coverage):
>> query1  test query
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  100.0   0.0   1.0e-10   1.0e-10       1      10 ..       1      10 ..       1      10 .. 0.99

  Alignments for each domain:
  == domain 1  score: 100.0 bits;  conditional E-value: 1.0e-10
  PTHR00001.1   1 abcde 5
                  ABCDE
       query1   1 ABCDE 5
                  99999

  PTHR00001.1   6 fghij 10
                  FGHIJ
       query1   6 FGHIJ 10
                  99999

//
"""


def test_multiblock_alignment_captured():
    """Alignment spanning two blocks should be fully captured."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.out', delete=False) as f:
        f.write(MULTIBLOCK_HMMER)
        f.flush()
        try:
            matches = parsehmmsearch(f.name)
        finally:
            os.unlink(f.name)

    assert 'PTHR00001' in matches
    assert 'query1' in matches['PTHR00001']
    align = matches['PTHR00001']['query1']
    # Full alignment should be 10 chars, not just 5 from first block
    assert align['hmmalign'][0] == 'abcdefghij'
    assert align['matchalign'][0] == 'ABCDEFGHIJ'
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_parsehmmsearch.py::test_multiblock_alignment_captured -v`

Expected: FAIL — current parser only captures first block: `hmmalign[0] == 'abcde'`

---

## Task 6: Fix `parsehmmsearch` to accumulate multi-block alignments

**Files:**
- Modify: `treegrafter.py:368-400`

- [ ] **Step 1: Replace the `== domain` parsing block**

The current code at lines 368-400 reads one block per domain. Replace it with a deterministic loop that reads exactly 4 lines per block (HMM, consensus, query, PP) plus a blank separator, then checks whether the next line starts a continuation block by validating the model name. This avoids the fragile peek-ahead heuristic.

**Key design decisions:**
- The inner `while True` loop reads the known 4-line HMMER3 block structure deterministically
- Continuation detection: after consuming a full block + blank separator, the next iteration reads a line and checks if its first token matches the model name
- When the loop breaks, `line` holds the first non-block line — the outer loop's `continue` skips its bottom `line = fp.readline()`, so this line gets processed next
- When `domain_num not in store_domain`, we fall through to only the outer loop's bottom `line = fp.readline()` (one read, not two)

Replace this code (lines 368-400):

```python
            elif m.match(r'\s+==\sdomain\s(\d+)') and store_align:
                domain_num = m.group(1)

                if domain_num in store_domain:
                    line = fp.readline()
                    hmmalign_array = line.split()

                    hmmalign_model = hmmalign_array[0]
                    hmmalign_model = re.sub(r'\..+', '', hmmalign_model)
                    hmmalign_seq = hmmalign_array[2]

                    line = fp.readline()
                    line = fp.readline()

                    matchalign_array = line.split()

                    matchalign_query = matchalign_array[0]
                    matchalign_query = stringify(matchalign_query)

                    matchlign_seq = matchalign_array[2]

                    if matchpthr == hmmalign_model and query_id == matchalign_query:
                        if len(current_match['align']['hmmalign']) >= len(current_match['align']['hmmstart']):
                            sys.stderr.write("Trying to add alignment sequence"
                                             " without additional data.\n")
                            sys.exit(1)

                        current_match['align']['hmmalign'].append(hmmalign_seq)
                        current_match['align']['matchalign'].append(matchlign_seq)

                    match_store[query_id] = current_match

                line = fp.readline()
```

With this new code:

```python
            elif m.match(r'\s+==\sdomain\s(\d+)') and store_align:
                domain_num = m.group(1)

                if domain_num in store_domain:
                    hmmalign_seq = ''
                    matchlign_seq = ''
                    hmmalign_model = None
                    matchalign_query = None

                    # Read all alignment blocks for this domain.
                    # HMMER3 outputs 4 lines per block:
                    #   1. HMM sequence line (e.g. "PTHR00001.1  1 abcde 5")
                    #   2. Consensus/match line
                    #   3. Query sequence line (e.g. "query1  1 ABCDE 5")
                    #   4. PP (posterior probability) line
                    # Blocks are separated by a blank line.
                    while True:
                        # Read HMM line (line 1 of block)
                        line = fp.readline()
                        if not line or not line.strip():
                            break
                        hmmalign_array = line.split()
                        if len(hmmalign_array) < 3:
                            break

                        # On continuation blocks, verify model name matches
                        block_model = re.sub(r'\..+', '', hmmalign_array[0])
                        if hmmalign_model is not None:
                            if block_model != hmmalign_model:
                                break  # not a continuation — stop
                        else:
                            hmmalign_model = block_model
                        hmmalign_seq += hmmalign_array[2]

                        fp.readline()            # consensus line (line 2)
                        line = fp.readline()      # query line (line 3)

                        matchalign_array = line.split()
                        if matchalign_query is None:
                            matchalign_query = stringify(matchalign_array[0])
                        matchlign_seq += matchalign_array[2]

                        fp.readline()            # PP line (line 4)
                        fp.readline()            # blank separator between blocks

                    # `line` now holds the first non-block line (next section
                    # or blank). The outer while-loop will process it via
                    # m = ReMatcher(line) at the top, so we must NOT read
                    # another line at the bottom of the outer loop.

                    if hmmalign_model and matchalign_query:
                        if matchpthr == hmmalign_model and query_id == matchalign_query:
                            if len(current_match['align']['hmmalign']) >= len(current_match['align']['hmmstart']):
                                sys.stderr.write("Trying to add alignment sequence"
                                                 " without additional data.\n")
                                sys.exit(1)

                            current_match['align']['hmmalign'].append(hmmalign_seq)
                            current_match['align']['matchalign'].append(matchlign_seq)

                    match_store[query_id] = current_match

                    continue  # skip bottom-of-loop readline; line already set

                # domain_num not in store_domain — just let the outer loop
                # advance via its bottom-of-loop readline (line 408)
```

- [ ] **Step 2: Run parsing tests**

Run: `python -m pytest tests/test_parsehmmsearch.py -v`

Expected: PASS

- [ ] **Step 3: Run all tests**

Run: `python -m pytest tests/ -v`

Expected: ALL PASS

---

## Task 7: Remove `filter_best_domain` single-domain restriction (optional)

**Decision: Keep `filter_best_domain` as-is. No changes needed.**

**Reasoning:** Downstream processing matches queries to a single PANTHER family for tree grafting (EPA-ng placement). Multiple domains from the same HMM don't change which family the query gets grafted to, so keeping only the best-scoring domain gives the strongest placement signal. The overlap handling in the rewritten `_querymsf` still acts as a safety net if `filter_best_domain` is ever bypassed or if edge cases slip through.

- [x] **Step 1: Discussed — single-domain behavior is intentional for downstream tree grafting**
- [x] **Step 2: N/A — keeping `filter_best_domain`**
- [x] **Step 3: N/A**

---

## Task 8: Integration test with real data

- [ ] **Step 1: Run treegrafter on the input that previously produced the error**

Run the same command that produced: `Error: length of query MSF longer than expected PANTHER alignment length: expected 157, got 82`

- [ ] **Step 2: Verify no length mismatch errors**

Expected: No "Error: length of query MSF" messages on stderr.

- [ ] **Step 3: Verify output is reasonable**

Compare output against a known-good reference if available. At minimum, confirm the output file is non-empty and contains expected columns.

---

## Notes

- The error message at line 145 says "longer" but the actual querymsf (82) is **shorter** than expected (157). Consider fixing the message to say "different from" instead.
- The `domscore` field is already collected during parsing (line 352), so no additional parsing changes are needed to support score-based domain ranking in `_querymsf`.
- The HMMER3 alignment block structure is: (1) HMM line, (2) consensus/match line, (3) query line, (4) PP line — with a blank line between blocks. The parser needs to handle this 4-line-per-block format.