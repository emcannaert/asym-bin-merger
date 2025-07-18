import numpy as np
import pytest
from asym_bin_merger import AsymBinMerger

# Define test cases to validate _get_neighbor_bins
neighbor_bin_cases = {
    "TC1_simple_neighbors": {
        "bin_contents": np.array([
            [10,  1, 10],
            [ 1,  1,  1],
            [10,  1, 10]
        ]),
        "expected_neighbors_for_bin_4": [
            # Index 4 = (1,1), should return its 8-connected surrounding bins
            [(0,1)], [(1,0)], [(1,2)], [(2,1)],
            [(0,0)], [(0,2)], [(2,0)], [(2,2)]
        ]
    },
    "TC2_isolated_bad_bin": {
        "bin_contents": np.array([
            [10, 10, 10],
            [10,  1, 10],
            [10, 10, 10]
        ]),
        "expected_neighbors_for_bin_4": [
            # Bin (1,1) is surrounded by high-content bins, but each is in a separate superbin
            [(0,1)], [(1,0)], [(1,2)], [(2,1)],
            [(0,0)], [(0,2)], [(2,0)], [(2,2)]
        ]
    },
    "TC3_only_edge_neighbors": {
        "bin_contents": np.array([
            [1, 1, 1],
            [1, 10, 1],
            [1, 1, 1]
        ]),
        "expected_neighbors_for_bin_4": [
            # Neighbors should be bins adjacent to (1,1), all with low counts
            [(0,1)], [(1,0)], [(1,2)], [(2,1)],
            [(0,0)], [(0,2)], [(2,0)], [(2,2)]
        ]
    }
}

edge_cases = {
    "TC4_corner_bin_neighbors": {
        "bin_contents": np.array([
            [ 1, 2, 3],
            [ 4, 5, 6],
            [ 7, 8, 9]
        ]),
        "target_bin": (0, 0),
        "expected_neighbors": [
            [(0,1)], [(1,0)], [(1,1)]
        ]
    },
    "TC5_edge_bin_neighbors": {
        "bin_contents": np.array([
            [ 1, 2, 3],
            [ 4, 5, 6],
            [ 7, 8, 9]
        ]),
        "target_bin": (1, 2),
        "expected_neighbors": [
            [(0,1)], [(0,2)], [(1,1)], [(2,1)], [(2,2)]
        ]
    }
}

@pytest.mark.parametrize("name, case", neighbor_bin_cases.items())
def test_get_neighbor_bins(name, case):
    # Create a merger with debug True so we trigger both code paths
    merger = AsymBinMerger(hist=case["bin_contents"], max_stat_uncert=0.25, output_dir="bin_maps/", debug=True)
    
    # Initialize the superbins and find bad bins
    merger._init_superbin_indices()
    bad_bins = merger._get_bad_bins()
    
    # Make sure we actually got the center bin (1,1), i.e., index 4 in a 3x3 flattened grid
    center_coord = (1, 1)
    try:
        bad_bin_index = merger.superbin_indices.index([center_coord])
    except ValueError:
        pytest.fail(f"Center bin {center_coord} not found in bad bins for {name}")
    
    neighbors = merger._get_neighbor_bins(bad_bin_index)
    
    # Flatten neighbor coords for easier comparison
    flat_neighbors = [sorted(list(sb)) for sb in neighbors]
    expected_neighbors = [sorted([coord]) for coord in case["expected_neighbors_for_bin_4"]]
    
    assert flat_neighbors == expected_neighbors, f"Failed {name}: Expected {expected_neighbors}, got {flat_neighbors}"

@pytest.mark.parametrize("name, case", edge_cases.items())
def test_edge_and_corner_neighbors(name, case):
    merger = AsymBinMerger(hist=case["bin_contents"], max_stat_uncert=0.25, output_dir="bin_maps/", debug=True)
    merger._init_superbin_indices()
    
    # Insert the corner/edge bin directly as a "bad bin" to simulate isolation
    bin_coord = case["target_bin"]
    try:
        index = merger.superbin_indices.index([bin_coord])
    except ValueError:
        merger.superbin_indices.append([bin_coord])
        index = len(merger.superbin_indices) - 1
    
    neighbors = merger._get_neighbor_bins(index)
    flat_neighbors = [sorted(list(sb)) for sb in neighbors]
    expected_neighbors = [sorted([coord]) for coord in case["expected_neighbors"]]
    
    assert flat_neighbors == expected_neighbors, f"Failed {name}: Expected {expected_neighbors}, got {flat_neighbors}"
