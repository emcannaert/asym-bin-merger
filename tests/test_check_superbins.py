import pytest
import numpy as np
from asym_bin_merger import AsymBinMerger  

def create_merger(hist_shape, superbin_indices):
    dummy_hist = np.ones(hist_shape)
    merger = AsymBinMerger(hist=dummy_hist, max_stat_uncert=0.1, output_dir=".", debug=True)
    merger.superbin_indices = superbin_indices
    return merger


@pytest.mark.parametrize("superbins", [
    [[(0, 0)]],                                         # single superbin with one bin
    [[(0, 0), (0, 1), (1, 1)]],                         # single superbin with multiple bins
    [[(0, 0)], [(0, 1,)], [(1, 0), (1, 1), (2, 2)]],    # multiple superbins with varied lengths
    [[(0, 0)], [(1, 1), (1, 2)], [(2, 0), (2, 1)]],     # multiple superbins with mixed lengths
])
def test_check_superbins_valid_varied_lengths(superbins):
    merger = create_merger(hist_shape=(3, 3), superbin_indices=superbins)
    # Should not raise exception
    merger._check_superbins()
    # Check that none of the superbins are empty
    assert all(len(sb) > 0 for sb in merger.superbin_indices)


@pytest.mark.parametrize("superbins", [
    [[], [(0, 1)]],                   # one empty superbin (should print warning, not fail)
    [[(0, 0)], "not a list"],         # invalid superbin structure (not list)
])
def test_check_superbins_invalid_structure(superbins):
    if not all(isinstance(sb, list) for sb in superbins):
        with pytest.raises(TypeError, match="superbin_indices must be a list of lists."):
            merger = create_merger(hist_shape=(2, 2), superbin_indices=superbins)
            merger._check_superbins()
    else:
        merger = create_merger(hist_shape=(2, 2), superbin_indices=superbins)
        merger._check_superbins()
        # After cleanup, empty superbins should be removed
        assert all(len(sb) > 0 for sb in merger.superbin_indices)


@pytest.mark.parametrize("bad_superbins", [
    [[(5, 5)]],               # bin index out of bounds (too large)
    [[(-1, 0)]],              # negative index
    [[(0, 3)]],               # y-index out of bounds
    [[(1, 1), (4, 4)]],       # mixed valid and invalid bins
])
def test_check_superbins_out_of_bounds(bad_superbins):
    merger = create_merger(hist_shape=(3, 3), superbin_indices=bad_superbins)
    with pytest.raises(Exception, match="Invalid bin index"):
        merger._check_superbins()
