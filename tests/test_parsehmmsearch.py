import sys
import os
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from treegrafter import parsehmmsearch

# Minimal HMMER3 output with a single query matching PTHR00001,
# one domain whose alignment spans two blocks.
# Format closely mirrors real hmmsearch output (see homo_sapiens_19_sample.hits.out).
MULTIBLOCK_HMMER = """\
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# query HMM file: test.hmm
# target sequence database: test.fa
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Query:       PTHR00001.orig.30.pir  [M=10]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence  Description
    ------- ------ -----    ------- ------ -----   ---- --  --------  -----------
    1.0e-10  100.0   0.0    1.0e-10  100.0   0.0    1.0  1  query1    test query


Domain annotation for each sequence (and alignments):
>> query1  test query
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  100.0   0.0   1.0e-10   1.0e-10       1      10 ..       1      12 ..       1      12 .. 0.99

  Alignments for each domain:
  == domain 1  score: 100.0 bits;  conditional E-value: 1.0e-10
                  PTHR00001.orig.30.pir  1 abcde 5
                                           ABCDE
                             query1  1 ABCDE 5
                                       99999 PP

                  PTHR00001.orig.30.pir  6 fghij 10
                                           FGHIJ
                             query1  6 FGHIJ 10
                                       99999 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (10 nodes)
Target sequences:                            1  (12 residues searched)
Passed MSV filter:                         1  (1.0); expected 0.0 (0.02)
Passed bias filter:                        1  (1.0); expected 0.0 (0.02)
Passed Vit filter:                         1  (1.0); expected 0.0 (0.001)
Passed Fwd filter:                         1  (1.0); expected 0.0 (1e-05)
Initial search space (Z):                  1  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: 0.00
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
