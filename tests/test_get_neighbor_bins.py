from unittest import case
import numpy as np
import pytest

from asym_bin_merger import AsymBinMerger

# Define test cases to validate _get_neighbor_bins
simple_neighbor_cases = {
    "TC1_simple_neighbors": {
        "bin_contents": np.array([[10, 1, 10], [1, 1, 1], [10, 1, 10]]),
        "expected_neighbors_for_bin_4": [
            # Simple case: Index 4 = (1,1), should return its 8-connected surrounding bins
            [(0, 1)],
            [(1, 0)],
            [(1, 2)],
            [(2, 1)],
            [(0, 0)],
            [(0, 2)],
            [(2, 0)],
            [(2, 2)],
        ],
    },
    "TC2_neighbors_higher": {
        "bin_contents": np.array([[10, 10, 10], [10, 1, 10], [10, 10, 10]]),
        "expected_neighbors_for_bin_4": [
            # Target Bin (1,1) is surrounded by same-valued bins, each in separate superbin
            [(0,0)], [(0,1)], [(0,2)], [(1,0)],
            [(1,2)], [(2,0)], [(2,1)], [(2,2)]
        ]

    },
    "TC3_neighbors_lower": {
        "bin_contents": np.array([[1, 1, 1], [1, 10, 1], [1, 1, 1]]),
        "expected_neighbors_for_bin_4": [
            # Target Bin (1,1) is surrounded by same-valued bins, but high values this time
            [(0,0)], [(0,1)], [(0,2)], [(1,0)],
            [(1,2)], [(2,0)], [(2,1)], [(2,2)]
        ]
    }

}

edge_cases = {
    "TC4_corner_bin_neighbors": {
        #Target bin is in corner
        "bin_contents": np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        "target_bin": (0, 0),
        "expected_neighbors": [[(0, 1)], [(1, 0)], [(1, 1)]],
    },
    "TC5_edge_bin_neighbors": {
        # Target bin is on non-corner edge
        "bin_contents": np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        "target_bin": (1, 2),
        "expected_neighbors": [[(0, 1)], [(0, 2)], [(1, 1)], [(2, 1)], [(2, 2)]],
    },
}

stress_cases = {
    "TC6_large_histogram_center": {
        #Target bin is in the center of a large histogram (10x10)
        "bin_contents": np.pad(np.ones((3, 3)), pad_width=3, mode='constant', constant_values=5),
        "target_bin": (4, 4),
        "expected_neighbors": [
            [(3, 3)], [(3, 4)], [(3, 5)],
            [(4, 3)],           [(4, 5)],
            [(5, 3)], [(5, 4)], [(5, 5)],
        ],
    },
    "TC7_high_values_3x3": {
        # Target bin is surrounded by high values 
        "bin_contents": np.array([
            [100, 100, 100],
            [100,   1, 100],
            [100, 100, 100]
        ]),
        "target_bin": (1, 1),
        "expected_neighbors": [
            [(0, 0)], [(0, 1)], [(0, 2)],
            [(1, 0)],           [(1, 2)],
            [(2, 0)], [(2, 1)], [(2, 2)],
        ],
    },
    "TC8_isolated_high_center": {
        #Target bin contains high values and is surrounded by low ones
        "bin_contents": np.array([
            [1, 1, 1],
            [1, 100, 1],
            [1, 1, 1]
        ]),
        "target_bin": (1, 1),
        "expected_neighbors": [
            [(0, 0)], [(0, 1)], [(0, 2)],
            [(1, 0)],           [(1, 2)],
            [(2, 0)], [(2, 1)], [(2, 2)],
        ],
    }
}

@pytest.mark.parametrize("name, case", simple_neighbor_cases.items())
def test_get_neighbor_bins(name, case):
    # Create a merger with debug True so we trigger both code paths
    merger = AsymBinMerger(
        hist=case["bin_contents"],
        max_stat_uncert=0.25,
        output_dir="bin_maps/",
        debug=True,
    )

    # Initialize the superbins and find bad bins
    merger._init_superbin_indices()
    bad_bins = merger._get_bad_bins()

    # Make sure we actually got the center bin (1,1),
    # i.e., index 4 in a 3x3 flattened grid
    center_coord = (1, 1)
    try:
        bad_bin_index = merger.superbin_indices.index([center_coord])
    except ValueError:
        pytest.fail(f"Center bin {center_coord} not found in bad bins for {name}")

    neighbors = merger._get_neighbor_bins(bad_bin_index)
    
    assert neighbors == case["expected_neighbors_for_bin_4"], \
    f"Failed {name}: Expected {case['expected_neighbors_for_bin_4']}, got {neighbors}"


@pytest.mark.parametrize("name, case", edge_cases.items())
def test_edge_and_corner_neighbors(name, case):
    merger = AsymBinMerger(
        hist=case["bin_contents"],
        max_stat_uncert=0.25,
        output_dir="bin_maps/",
        debug=True,
    )
    merger._init_superbin_indices()

    # Insert the corner/edge bin directly as a "bad bin" to simulate isolation
    bin_coord = case["target_bin"]
    try:
        index = merger.superbin_indices.index([bin_coord])
    except ValueError:
        merger.superbin_indices.append([bin_coord])
        index = len(merger.superbin_indices) - 1

    neighbors = merger._get_neighbor_bins(index)
    
    assert neighbors == case["expected_neighbors"], \
    f"Failed {name}: Expected {case['expected_neighbors']}, got {neighbors}"


@pytest.mark.parametrize("name, case", stress_cases.items())
def test_stress_neighbors_known(name, case):
    merger = AsymBinMerger(
        hist=case["bin_contents"],
        max_stat_uncert=0.25,
        output_dir="bin_maps/",
        debug=False,
    )
    merger._init_superbin_indices()

    coord = case["target_bin"]
    try:
        index = merger.superbin_indices.index([coord])
    except ValueError:
        merger.superbin_indices.append([coord])
        index = len(merger.superbin_indices) - 1

    neighbors = merger._get_neighbor_bins(index)

    # Sort each list of coordinates so that order doesn't affect comparison
    sorted_actual = [sorted(group) for group in neighbors]
    sorted_expected = [sorted(group) for group in case["expected_neighbors"]]

    assert sorted_actual == sorted_expected, \
        f"{name} failed: Expected {sorted_expected}, got {sorted_actual}"

