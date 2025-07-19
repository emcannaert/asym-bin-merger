import numpy as np
import pytest

from asym_bin_merger import AsymBinMerger

test_cases_with_bottom_left_origin = {
    "test1": {
        "bin_contents": np.array(
            [[100, 200, 300], [400, 1, 500], [600, 700, 800]]
        ),  # row 0 = bottom
        "merged_content": np.array([[100, 201, 300], [400, 201, 500], [600, 700, 800]]),
    },
    "test2": {
        "bin_contents": np.array([[1, 100, 150], [50, 1000, 300], [200, 250, 350]]),
        "merged_content": np.array([[51, 100, 150], [51, 1000, 300], [200, 250, 350]]),
    },
    "test3": {
        "bin_contents": np.array([[1, 50, 150], [50, 1000, 300], [200, 250, 350]]),
        "merged_content": np.array([[51, 51, 150], [50, 1000, 300], [200, 250, 350]]),
    },
    "test4": {
        "bin_contents": np.array([[999, 999, 999], [999, 1, 999], [999, 999, 999]]),
        "merged_content": np.array(
            [[999, 1000, 999], [999, 1000, 999], [999, 999, 999]]
        ),
    },
    "test5": {
        "bin_contents": np.array([[999, 999, 999], [999, 1, 999], [999, 999, 999]]),
        "merged_content": np.array(
            [[999, 1000, 999], [999, 1000, 999], [999, 999, 999]]
        ),
    },
    "test6": {
        "bin_contents": np.array([[1000, 2, 50], [1000, 5, 10], [1000, 1000, 1000]]),
        "merged_content": np.array(
            [[1000, 17, 50], [1000, 17, 17], [1000, 1000, 1000]]
        ),
    },
    "test7": {
        "bin_contents": np.array([[1000, 2, 50], [1000, 5, 5], [1000, 1000, 1000]]),
        "merged_content": np.array(
            [[1000, 62, 62], [1000, 62, 62], [1000, 1000, 1000]]
        ),
    },
    "test8": {
        "bin_contents": np.array([[1, 2, 3], [5, 1000, 6], [4, 3, 2]]),
        "merged_content": np.array([[26, 26, 26], [26, 1000, 26], [26, 26, 26]]),
    },
}


@pytest.mark.parametrize("name, case", test_cases_with_bottom_left_origin.items())
def test_bin_merging(name, case):
    merger = AsymBinMerger(
        hist=case["bin_contents"],
        max_stat_uncert=0.25,
        output_dir="bin_maps/",
        debug=True,
        verbose=True,
    )
    merger.run()
    np.testing.assert_array_equal(merger._get_merged_hist(), case["merged_content"])
