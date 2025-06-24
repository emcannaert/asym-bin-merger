"""
asym_bin_merger: a utility that finds the optimized (= largest number of 
bins with statistical uncertainty under some threshold) bin mappings for a given 
ROOT histogram. Bin mappings are saved in a text file with binning format
    [ [ (x1,y1),(x2,y2),... ], [(xi,yi),(xj,yj),...], ...  ] 
where entries in the outer list represent post-merged bins (="superbins") in terms 
of their constituents written as tuples containing their x- and y-coordinates
in the original 2D distribution.

Inputs:
    hist: a ROOT histogram over which to calculate binning
    max_stat_uncert: the maximum statistical uncertainty threshold for any bin 
    output_dir (optional): the location where the output text file is saved, default=
    debug (optional): testing par that bypasses main workflow, default=False
Outputs:
    Bin mappings are stored as bin_maps/<some name>_bin_map.txt by default, and 
    can also be obtained inline using the .get_bin_map() method.


"""
import ROOT  ## should be in documentation that this can only be run in a ROOT environment 
import os
import numpy as np

class AsymBinMerger:
    def __init__(self, hist=None, max_stat_uncert=None, output_dir=None, debug: bool = False):
        self.hist = hist
        self.max_stat_uncert = max_stat_uncert
        self.output_dir = output_dir
        self.debug = debug
        self.file_name = "bin_map.txt"
        self.output_file = os.path.join(self.output_dir, self.file_name)
        self.final_hist = np.zeros((5,5)) # placeholder for final hist
        self.debug = debug

        #check inputs
        self._validate_inputs()

        #run main workflow
        if not debug: 
            self._run()
            self._write_output()
        else:
            print("Running in debug mode.")




    ## Main methods
    def _convert_hist(self) -> list:  # convert hist to 2D array for easier handling, also set final_hist
        # TODO
        return []
    
    def _init_superbin_indices(self) -> list: # initialize the superbin index list
        # TODO
        return []
    
    def _check_superbins(self): # check whether there are problems with final superbins
        # TODO
        print("Checking bin map for issues.")
        return
    
    def _run_bin_merging(self): # run main bin merging scheme
        # TODO
        print("Run bin merging scheme.")
        return

    def _write_output(self): # write bin map to self.output_file
        # TODO
        print("Writing output to %s."%(self.output_file))
        return
    
    def _run(self): # do main bin-merging scheme
            self.converted_hist = self._convert_hist()
            self.superbin_indices = self._init_superbin_indices()
            self._run_bin_merging()
            self._check_superbins()


    ## Helper functions

    def _print_help(self): # give information about usage
        # TODO
        return
    
    def _get_bad_bins(self) -> list:# return list (by decreasing stat uncert) of current "bad" bins 
        # TODO
        return []

    def _get_neighbor_bins(self) -> list: # return list (by decreasing stat uncert) of neighbor superbins
        # TODO
        return []
    
    def print_bin_map(self): # helper for testing and checking merging process
        # TODO
        print("Printing out formatted bin map")
        return


    ## Getters and setters

    def get_bin_map(self) -> list: # external getter for superbin indices list 
        # TODO
        return []


    ## testing and meta functions
    def _validate_inputs(self):
        if self.debug:
            if not isinstance(self.hist, list) and not isinstance(self.hist, np.array):
                raise TypeError("In debug mode `hist` must be a np.array or list object.")
        else:
            if not isinstance(self.hist, ROOT.TH2):
                raise TypeError("`hist` must be a ROOT.TH2 object.")
        if not isinstance(self.max_stat_uncert, float):
            raise TypeError("`max_stat_uncert` must be a float.")
        if not isinstance(self.output_dir, str):
            raise TypeError("`output_dir` must be a string.")
    
    def _get_merged_hist(self): # testing: return np.array version of post-merge hist
        # TODO
        print("Returning post-merged hist as numpy array.")
        return []