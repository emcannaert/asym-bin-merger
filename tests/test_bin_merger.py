import numpy as np
import pytest
from asym_bin_merger import AsymBinMerger 

test_cases = {
    "TC1_one_bad_bin_center": {
        "bin_contents": np.array([
            [20, 20, 20],
            [20,  1, 20],
            [20, 20, 20]
        ]),
        "merged_contents": np.array([
            [20, 20, 20],
            [20, 21, 20],
            [20, 21, 20]
        ])
    },
}

more_test_cases = {
	"TC5_snake_bad_path": {
    	# Bad bins form a zig-zag "snake" that should cluster together
    	"bin_contents": np.array([
        	[50,  1, 50],
        	[ 1,  1,  1],
        	[50,  1, 50]
    	]),
    	"merged_contents": np.array([
        	[50, 17, 50],
        	[17, 17, 17],
        	[50, 17, 50]
    	])
	},
	"TC6_corner_spill": {
    	# Central cluster is bad, merges into adjacent good corner
    	"bin_contents": np.array([
        	[25,  2, 25],
        	[ 2,  1,  2],
        	[25,  2, 25]
    	]),
    	"merged_contents": np.array([
        	[25, 32, 25],
        	[32, 32, 32],
        	[25, 32, 25]
    	])
	},
	"TC7_uniform_low_3x3": {
    	# Uniform values below threshold, expect merging into blocks
    	"bin_contents": np.array([
        	[10, 10, 10],
        	[10, 10, 10],
        	[10, 10, 10]
    	]),
    	"merged_contents": np.array([
        	[20, 20, 20],
        	[20, 20, 20],
        	[20, 20, 20]
    	])
	},
	"TC8_bad_ring_around_good_center": {
    	# Good bin in center, bad ring around it should not disturb center
    	"bin_contents": np.array([
        	[ 2,  2,  2],
        	[ 2, 20,  2],
        	[ 2,  2,  2]
    	]),
    	"merged_contents": np.array([
        	[10, 10, 10],
        	[10, 20, 10],
        	[10, 10, 10]
    	])
	},
}


@pytest.mark.parametrize("name, case", test_cases.items())
def test_bin_merging(name, case):
    merger = AsymBinMerger(hist=case["bin_contents"], max_stat_uncert=0.25, output_dir="bin_maps/", debug = True)
    merger._run()
    assert np.array_equal(merger._get_merged_hist(), case["merged_contents"]) == True
    np.testing.assert_array_equal(merger._get_merged_hist(), case["merged_contents"])
