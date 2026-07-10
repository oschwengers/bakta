import pytest
import numpy as np

from bakta import plot


def test_calc_gc_content_short_sequence():
    """Test calc_gc_content with very short sequences that trigger step_size = 0 bug"""
    # Test with single base sequence
    seq = "G"
    pos_list, gc_content_list = plot.calc_gc_content(seq)
    assert len(pos_list) > 0
    assert len(gc_content_list) > 0
    assert isinstance(pos_list, np.ndarray)
    assert isinstance(gc_content_list, np.ndarray)


def test_calc_gc_content_two_bases():
    """Test calc_gc_content with two base sequence"""
    seq = "GC"
    pos_list, gc_content_list = plot.calc_gc_content(seq)
    assert len(pos_list) > 0
    assert len(gc_content_list) > 0
    assert isinstance(pos_list, np.ndarray)
    assert isinstance(gc_content_list, np.ndarray)


def test_calc_gc_content_normal_sequence():
    """Test calc_gc_content with a normal length sequence"""
    seq = "ATGCATGCATGCATGC" * 100  # 1600 bases
    pos_list, gc_content_list = plot.calc_gc_content(seq)
    assert len(pos_list) > 0
    assert len(gc_content_list) > 0
    assert isinstance(pos_list, np.ndarray)
    assert isinstance(gc_content_list, np.ndarray)
    # Verify that step_size is not zero for this length
    # With len(seq) = 1600, step_size = int(1600/1000) = 1
    assert len(pos_list) > 1


def test_calc_gc_skew_short_sequence():
    """Test calc_gc_skew with very short sequences that trigger step_size = 0 bug"""
    # Test with single base sequence
    seq = "G"
    pos_list, gc_skew_list = plot.calc_gc_skew(seq)
    assert len(pos_list) > 0
    assert len(gc_skew_list) > 0
    assert isinstance(pos_list, np.ndarray)
    assert isinstance(gc_skew_list, np.ndarray)


def test_calc_gc_skew_two_bases():
    """Test calc_gc_skew with two base sequence"""
    seq = "GC"
    pos_list, gc_skew_list = plot.calc_gc_skew(seq)
    assert len(pos_list) > 0
    assert len(gc_skew_list) > 0
    assert isinstance(pos_list, np.ndarray)
    assert isinstance(gc_skew_list, np.ndarray)


def test_calc_gc_skew_normal_sequence():
    """Test calc_gc_skew with a normal length sequence"""
    seq = "ATGCATGCATGCATGC" * 100  # 1600 bases
    pos_list, gc_skew_list = plot.calc_gc_skew(seq)
    assert len(pos_list) > 0
    assert len(gc_skew_list) > 0
    assert isinstance(pos_list, np.ndarray)
    assert isinstance(gc_skew_list, np.ndarray)
    # Verify that step_size is not zero for this length
    assert len(pos_list) > 1
