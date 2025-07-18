import pytest
import numpy as np
from asym_bin_merger import AsymBinMerger

@pytest.mark.parametrize("hist, is_valid", [
    # Valid arguments

    # Regular 2D NumPy array
    (np.array([[1, 2], [3, 4]]), True),
    # Jagged NumPy array
    (np.array([np.array([1, 2]), np.array([3])], dtype=object), True),

    # Invalid arguments

    # List of lists (not a NumPy array)
    ([[1, 2], [3]], False),
    # Scalar input
    (42, False),
    # String input
    ("not an array", False),
    # Flat NumPy array
    (np.array([1, 2, 3]), False),
    # Empty NumPy array
    (np.array([]), False),
    # 3D NumPy array
    (np.zeros((2, 2, 2)), False),
])
def test_get_merged_hist_accepts_valid_inputs(hist, is_valid):
    try:
        merger = AsymBinMerger(hist=hist, max_stat_uncert=0.1, output_dir=".", debug=True)
        # Dummy superbin_indices and final_hist for testing
        merger.superbin_indices = [[(0, 0)]]
        
        if isinstance(hist, np.ndarray) and hist.dtype == object:
            merger.final_hist = np.array([np.zeros_like(row) for row in hist], dtype=object)
        elif isinstance(hist, np.ndarray) and hist.ndim == 2:
            merger.final_hist = np.zeros_like(hist)
        else:
            # Just to avoid attribute errors in the test; will error out before use
            merger.final_hist = hist

        if is_valid:
            result = merger._get_merged_hist()
            assert isinstance(result, np.ndarray)
        else:
            with pytest.raises((TypeError, ValueError, AttributeError)):
                _ = merger._get_merged_hist()

    except Exception as e:
        if is_valid:
            pytest.fail(f"Unexpected error for valid input: {e}")
