import numpy as np
import pytest
from asym_bin_merger import AsymBinMerger 


superbin_cases = {
    "TC1_one_bad_bin_center": {
        "bin_contents": np.array([
            [20, 20, 20],
            [20,  1, 20],
            [20, 20, 20]
        ]),
        "merged_contents": [
            [(0,0)], [(0,1)], [(0,2)],
            [(1,0)], [(1,1)], [(1,2)],
            [(2,0)], [(2,1)], [(2,2)]
        ]
    },
    "TC2_no_bad_bin_center": {
        "bin_contents": np.array([
            [20, 20, 20],
            [20, 20, 20],
            [20, 20, 20]
        ]),
        "merged_contents": [
            [(0,0)], [(0,1)], [(0,2)],
            [(1,0)], [(1,1)], [(1,2)],
            [(2,0)], [(2,1)], [(2,2)]
        ]
    },
    "TC3_all_bad_bins": {
        "bin_contents": np.array([
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1]
        ]),
        "merged_contents": [
            [(0,0)], [(0,1)], [(0,2)],
            [(1,0)], [(1,1)], [(1,2)],
            [(2,0)], [(2,1)], [(2,2)]
        ]
    },
    "TC4_random_bad_bins": {
        "bin_contents": np.array([
            [20, 1, 20],
            [.5, 1, 3],
            [2, 50, 1]
        ]),
        "merged_contents": [
            [(0,0)], [(0,1)], [(0,2)],
            [(1,0)], [(1,1)], [(1,2)],
            [(2,0)], [(2,1)], [(2,2)]
        ]
    },
}


badbin_cases = {
    "TC1_one_bad_bin_center": {
        "bin_contents": np.array([
            [20, 20, 20],
            [20,  1, 20],
            [20, 20, 20]
        ]),
        "merged_contents": [
            [(1,1)]
        ]
    },
    "TC2_no_bad_bin_center": {
        "bin_contents": np.array([
            [20, 20, 20],
            [20, 20, 20],
            [20, 20, 20]
        ]),
        "merged_contents": []
    },
    "TC3_all_bad_bin_center": {
        "bin_contents": np.array([
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1]
        ]),
        "merged_contents": [
            [(1,1)], [(2,1)], [(1,2)],
            [(1,0)], [(0,1)], [(2,2)],
            [(2,0)], [(0,2)], [(0,0)]
        ]
    },
    "TC4_random_bad_bins": {
        "bin_contents": np.array([
            [20, 1, 20],
            [.5, 1, 3],
            [2, 50, 1]
        ]),
        "merged_contents": [
            [(1,0)], [(1,1)], [(0,1)], [(2,2)], [(2,0)], [(1,2)]
        ]
    },
}



#Test the "_init_superbin_indices" method and "get_bad_bins" method
@pytest.mark.parametrize("name, case", superbin_cases.items())
def test_init_superbin_indices(name, case):
    merger = AsymBinMerger(hist=case["bin_contents"], max_stat_uncert=0.25, output_dir="bin_maps/", debug = True)
    output = merger._init_superbin_indices()
    assert output == case["merged_contents"]

@pytest.mark.parametrize("name, case", badbin_cases.items())
def test_get_bad_bins(name, case):
    merger = AsymBinMerger(hist=case["bin_contents"], max_stat_uncert=0.25, output_dir="bin_maps/", debug = True)
    merger._init_superbin_indices()
    bad_bin_indices, bad_bin_nums = merger._get_bad_bins()
    assert bad_bin_indices == case["merged_contents"]