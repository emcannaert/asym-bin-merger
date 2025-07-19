# AsymBinMerger

**AsymBinMerger** is a Python plugin/class designed to process 2D ROOT histograms (`TH2F`) by merging input TH2F bins into "superbins" such that **no resulting output superbin has a fractional statistical uncertainty above a specified threshold**, while **maximizing the number of bins retained**.

This is particularly useful for analyses that require high-statistics binning in 2D phase space where this process is non-trivial.

---

## Features

- Merges bins based on **fractional statistical uncertainty** (`1/√counts`)
- Respects 2D histogram topology using **adjacent neighbor merging**
- "x-biased" tie-breaking ensures **deterministic, reproducible** results
- Outputs:
  - Superbin-to-bin map
  - Diagnostic plot

---

## Algorithm Overview

### Inputs

- `ROOT.TH2F` histogram
- `float`: maximum fractional statistical uncertainty (`max_stat_uncert`)
- `str` *(optional)*: output directory for files and plots, default is "merged_bins/"

### Merging Procedure

1. **Initialize** each bin as its own superbin.
2. **Identify bad bins**: superbins whose uncertainty exceeds the threshold.
3. **Choose the "worst" bad superbin using**:
     1. Highest fractional uncertainty
     2. Superbin centroid furthest from (0,0)
     3. Larger `x` coordinate, then larger `y` if tied
4. **Merge** the bad superbin into a **neighboring superbin**:
   - Select the neighbor with:
     1. Highest fractional uncertainty
     2. Furthest centroid from origin (0,0), favoring larger `x`, then `y`
5. **Update**: Recalculate list of "bad" superbins and repeat from 3 until no bad bins remain.

### Notes

- Neighbors are defined as superbins that have at least one constituent that is adjacent (up/down/left/right) to at least one constituent of the candidate superbin.
- Tie-breaking rules ensure **deterministic x-biased merging**.

---

## Example Usage

```python
from asym_bin_merger import AsymBinMerger
import ROOT

# Open ROOT file and retrieve histogram
r_file = ROOT.TFile.Open("path/to/your/hist.root", "READ")
hist = r_file.Get("my_hist")  # Ensure it's a TH2F

# Merge with max stat uncertainty of 12%
bin_merger = AsymBinMerger(hist, 0.12, "my_hist")

# Retrieve outputs
bin_merger.write_output()
bin_merger.merged_hist_to_image()

```

---

## Output Files

Stored in the `output_dir` (e.g. `"merged_hist"`):

- `bin_map.txt`: A list of lists with format:

  ```python
  [
    [(i1, j1), (i2, j2), ...],  # Superbin 0
    [(k1, l1), (k2, l2), ...],  # Superbin 1
    ...
  ]
  ```

- `merged_hist.png`: An image containing the merged histogram with statistical uncertainty labeled for each superbin

---

## Installation

This package uses [`uv`](https://github.com/astral-sh/uv) for environment and dependency management.

> **Note:** `ROOT` must be installed **separately** on your system (e.g., via [conda](https://root.cern/install/) or prebuilt binaries). This tool **requires** a working Python-compatible ROOT environment (`import ROOT` must work).

### Dependencies

- `numpy`
- `os` (standard library)
- `ROOT` *(not installable via `uv`)*
- `matplotlib` (for merged_hist_to_image)

### Setup Instructions

```bash
# Clone the repository
git clone https://github.com/your-org/asym-bin-merger.git
cd asym-bin-merger

# Create and enter a uv virtual environment
uv venv
source .venv/bin/activate

# Install Python dependencies
uv sync
```


---

## API Reference

### `AsymBinMerger(hist: ROOT.TH2F, max_stat_uncert: float, output_dir: str = None)`

Initializes the merger and performs bin merging.

### `.get_merged_hist() -> ROOT.TH2F`

Returns the post-merged histogram.

### `.get_bin_map() -> list[list[tuple[int, int]]]`

Returns the bin grouping (superbin) structure.

### `.write_output()`

Generates bin map and saves it in the output directory.

### `.merged_hist_to_image()`

Generates plot and saves it in the output directory.


---

## FAQ

**Q:** What happens if two bins have equal uncertainty and distance from the origin?
**A:** The algorithm prefers the bin with **larger x**, then **larger y** coordinate to ensure deterministic merging. This algorithm is therefore "x-biased".

**Q:** Can I use this with histograms that contain empty bins?
**A:** Yes. Empty bins (with zero entries) are treated as having infinite uncertainty and are prioritized for early merging.

---

## License

MIT License

---

## 🤝 Acknowledgments

Developed by Ethan Cannaert, Mayuri Kawale, Elias Mettner, and Ashling Quinn for use in high-statistics ROOT-based 2D analyses.
