import pytest
import numpy as np
from asym_bin_merger import AsymBinMerger 

@pytest.mark.parametrize("input_array", [
    np.array([[0, 0], [0, 0]]),                    # all zeros
    np.array([[1, 2], [3, 4]]),                    # simple numbers
    np.ones((3, 3)),                               # all ones
    np.random.randint(1, 10, size=(4, 2)),         # rectangular random
])
def test_convert_hist_with_valid_numpy_arrays(input_array):
    merger = AsymBinMerger(hist=input_array, max_stat_uncert=0.1, output_dir=".", debug=True)
    result = merger._convert_hist()

    assert isinstance(result, np.ndarray)
    assert result.shape == input_array.shape
    np.testing.assert_array_equal(result, input_array)


@pytest.mark.parametrize("invalid_input", [
    np.array([1, 2, 3]),      # 1D array
    "not an array",           # string
    [[1, 2], [3, 4]],         # list of lists
    123,                      # int
    None,                     # NoneType
])
def test_convert_hist_raises_typeerror_on_invalid_input(invalid_input):
    with pytest.raises(TypeError, match="In debug mode, `hist` must be a numpy array."):
        AsymBinMerger(hist=invalid_input, max_stat_uncert=0.1, output_dir=".", debug=True)._convert_hist()

