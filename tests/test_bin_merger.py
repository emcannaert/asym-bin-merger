import numpy as np
import pytest
from asym_bin_merger import AsymBinMerger 

test_cases_with_bottom_left_origin = {
    "TC1_one_bad_bin_center": {
        "bin_contents": np.array([
            [20, 20, 20],
            [20,  1, 20],
            [20, 20, 20]
        ]),  # row 0 = bottom
        "superbin_labels": np.array([
            [0, 1, 2],
            [3, 6, 4],
            [5, 6, 7]
        ])
    },
    "TC5_snake_bad_path": {
        "bin_contents": np.array([
            [50,  1, 50],
            [ 1,  1,  1],
            [50,  1, 50]
        ]),
        "superbin_labels": np.array([
            [0, 3, 1],
            [3, 3, 3],
            [2, 3, 3]
        ])
    },
    "TC6_corner_spill": {
        "bin_contents": np.array([
            [25,  2, 25],
            [ 2,  1,  2],
            [25,  2, 25]
        ]),
        "superbin_labels": np.array([
            [0, 3, 1],
            [3, 3, 3],
            [2, 3, 3]
        ])
    },
    "TC7_uniform_low_3x3": {
        "bin_contents": np.array([
            [10, 10, 10],
            [10, 10, 10],
            [10, 10, 10]
        ]),
        "superbin_labels": np.array([
            [1, 0, 0],
            [1, 2, 2],
            [1, 3, 3]
        ])
    },
    "TC8_bad_ring_around_good_center": {
        "bin_contents": np.array([
            [ 2,  2,  2],
            [ 2, 20,  2],
            [ 2,  2,  2]
        ]),
        "superbin_labels": np.array([
            [0, 0, 0],
            [0, 1, 0],
            [0, 0, 0]
        ])
    },
}


@pytest.mark.parametrize("name, case", test_cases_with_bottom_left_origin.items())
def test_bin_merging(name, case):
    merger = AsymBinMerger(hist=case["bin_contents"], max_stat_uncert=0.25, output_dir="bin_maps/", debug = True)
    merger._run()
    assert np.array_equal(merger._get_merged_hist(), case["merged_contents"]) == True
    np.testing.assert_array_equal(merger._get_merged_hist(), case["merged_contents"])
