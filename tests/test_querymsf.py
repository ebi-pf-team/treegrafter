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
    assert len(result) == pthr_align_length
    # Domain 1 fills positions 1-5, rest are gaps
    assert result == 'ABCDE---'