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
        self.superbin_indices = []  # list of lists, where each inner list contains indices of bins that are merged together
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
        # Identify bins that are already under the max_stat_uncert threshold
        # and set them as their own superbins. Superbins are initialized as a list of lists, 
        # where each inner list contains the indices of bins that are merged together.
        if self.debug:
            if isinstance(self.hist, np.ndarray):
                # For debug mode, assume hist is a numpy array
                superbins = []
                for i in range(self.hist.shape[0]):
                    for j in range(self.hist.shape[1]):
                        superbins.append([(i, j)])
            else:
                raise TypeError("In debug mode, `hist` must be a numpy array.")
        elif isinstance(self.hist, ROOT.TH2):
            # For ROOT histograms, we need to iterate through the bins
            superbins = []
            for i in range(1, self.hist.GetNbinsX() + 1):
                for j in range(1, self.hist.GetNbinsY() + 1):
                    superbins.append([(i, j)])
        else:
            raise TypeError("`hist` must be a ROOT.TH2 object or a numpy array in debug mode.")
        # If no bins are under the threshold, return an empty list
        if not superbins:
            print("No bins under the max_stat_uncert threshold. Returning empty superbins list.")
            self.superbin_indices = []
            return []
        else:
            print(f"Initialized superbins with {len(superbins)} entries.")
        self.superbin_indices = superbins
        return superbins

    def _check_superbins(self): # check whether there are problems with final superbins
        # TODO
        print("Checking bin map for issues.")
        return
    
    def _run_bin_merging(self): # run main bin merging scheme
        
        # init bad_bins
        bad_bin_nums,_ = self._get_bad_bins()
        it_num = 0

        print("Running bin merging sequence.")

        while len(bad_bins) > 0:

            # get largest stat uncert superbin number
            bad_bin_num = bad_bin_nums[0]

            # get list of neighbors 
            bad_neighbor_nums = self._get_neighbor_bins(bad_bin)
            if not bad_neighbor_nums:
                raise ValueError("(Iteration %s) No neighbor superbins found for superbin number %s with superbins indices %s."%(it_num, bad_bin_num, self.superbin_indices))

            bad_neighbor_num = bad_neighbor_nums[0]

            # add bad bin into 'worst' neighbor
            self.superbin_indices[bad_neighbor_num].extend(self.superbin_indices[bad_bin_num])

            # remove bad bin from superbin_indices
            self.superbin_indices.pop(bad_neighbor_num)

            # update bad bins list for next iteration
            bad_bin_nums,_ = self._get_bad_bins()
            it_num+=1
        print("Finished bin merging sequence - took %s iterations."%(it_num))

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
        print("Hello world.")
        return
    
    def _get_bad_bins(self) -> list:
    # return list of bad bins (those that are above the max_stat_uncert threshold)
    # step one: access current superbin list and check each bin's stat uncertainty. 
    # Note that one superbin can contain multiple bins. Must find total uncertainty for 
    # each superbin. A superbin containing multiple bins can be considered bad if the
    # total uncertainty is above the threshold. Bins are sorted by decreasing stat uncertainty.
        if not self.superbin_indices:
            print("No superbins initialized. Please run _init_superbin_indices() first.")
            return []
        bad_bins = []
        for superbin in self.superbin_indices:
            total_superbin_value = 0
            for bin_index in superbin:
                if self.debug:
                    # In debug mode, assume hist is a numpy array
                    bin_value = self.hist[bin_index[0], bin_index[1]]
                else:
                    # For ROOT histograms, use GetBinContent
                    bin_value = self.hist.GetBinContent(bin_index[0], bin_index[1])
                # Calculate statistical uncertainty (for simplicity, assume Poisson statistics)
                total_superbin_value += bin_value
            total_stat_uncert = np.sqrt(total_superbin_value) / total_superbin_value if total_superbin_value > 0 else 0
            if total_stat_uncert > self.max_stat_uncert:
                bad_bins.append(superbin)
        if self.debug:
            print(f"Found {len(bad_bins)} bad bins exceeding the max_stat_uncert threshold.")
            # sort bad bins by decreasing stat uncertainty
            bad_bins.sort(key=lambda x: -np.sqrt(np.sum([self.hist[bin_index[0], bin_index[1]] for bin_index in x])) / np.sum([self.hist[bin_index[0], bin_index[1]] for bin_index in x]))
        else:
            # sort bad bins by decreasing stat uncertainty
            bad_bins.sort(key=lambda x: -np.sqrt(np.sum([self.hist.GetBinContent(bin_index[0], bin_index[1]) for bin_index in x])) / np.sum([self.hist.GetBinContent(bin_index[0], bin_index[1]) for bin_index in x]))
            print(f"Found {len(bad_bins)} bad bins exceeding the max_stat_uncert threshold in ROOT histogram.")
        return bad_bins  # return list of bad bins, each represented as a list of indices

    def _get_neighbor_bins(self, bad_index) -> list: # return list (by decreasing stat uncert) of neighbor superbin indices
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
            if not isinstance(self.hist, list) and not isinstance(self.hist, np.ndarray):
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