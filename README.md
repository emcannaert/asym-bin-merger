# asym-bin-merger

#[![PyPI - Version](https://img.shields.io/pypi/v/asym-bin-merger.svg)](https://pypi.org/project/asym-bin-merger)
#[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/asym-bin-merger.svg)](https://pypi.org/project/asym-bin-merger)

-----

## Table of Contents
- [Description]

Inputs:       numpy arrays → flipped relative to normal lists, stretch: root histograms
Outputs:   text file: [ [(i1,j1), (i2,j2,), ... ], [], []    ]
can be written as [ superbin1 = ( bin1, bin2, bin3, ... ), superbin2 = (binx, biny, binz, ...), ... ]
→ bin1 = (x1,y1)
stretch: add conversion to colored histogram, JSON?

Workflow: 
STRETCH: compare different algorithms
STRETCH: implement for multiple backgrounds

1. accept ROOT histograms,  old: (Accept numpy arrays as pkl)
2. Run bin merging scheme
    1. while len(bad_bins) > 0:
    2. identify "bad" bins = bins with stat uncertainty > some threshold

        * loop over all bins
            * calculate the stat uncertainty in each bin
    * choose a bad bin
        * choose randomly, from highest stat uncertainty to lowest, or lowest to highest
    * look at the neighbors of bad bin
        * STRETCH: add in some type of weight based on distance to neighbor centroids
    * add neighbor bin with the highest stat uncertainty with bad bin 
    * remove neighbor bin w/ highest stat uncertainty
        * STRETCH: try to "complete" a bad bin before moving on to another
        * STRETCH: fill in holes

1. Output bin map 


Classes: 
bin_calculator → inputs: stat uncertainty threshold, histograms, location of output text file 
hist_container → histograms & weigh

List of methods: 
* help
    * STRETCH: text interface 
* convert_root_to_list: takes in 2D histogram, returns corresponding numpy list
* init superbin: loop over all bins in 2D bins and create 1D list of superbins from these:
    * superbins_indices = [ [1], [2], [3]  ]   → bin1.extend(bin2)
* run_bin_merging
* identify bad bins: returns SORTED python list of superbin indices of "bad" bins
* get list of neighbor bins: for given superbin, get SORTED (by stat uncertainty) list of neighbor bins → do we want all superbins that border any bin in the current superbin?
* output final bin map: create a text file at the specified output location



- [Installation](#installation)
- [License](#license)

## Installation

```console
pip install asym-bin-merger
```

## License

`asym-bin-merger` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
