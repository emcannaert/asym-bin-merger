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
        # Convert the ROOT histogram to a 2D numpy array for easier handling
        if self.debug:
            if isinstance(self.hist, np.ndarray):
                # If hist is already a numpy array, just return it
                return self.hist
            else:
                raise TypeError("In debug mode, `hist` must be a numpy array.")
        elif isinstance(self.hist, ROOT.TH2):
            # If hist is a ROOT histogram, convert it to a numpy array
            n_bins_x = self.hist.GetNbinsX()
            n_bins_y = self.hist.GetNbinsY()
            hist_array = np.zeros((n_bins_x, n_bins_y))
            for i in range(1, n_bins_x + 1):
                for j in range(1, n_bins_y + 1):
                    hist_array[i - 1, j - 1] = self.hist.GetBinContent(i, j)
            self.final_hist = hist_array
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
        bad_bin_indices,bad_bin_nums = self._get_bad_bins()
        it_num = 0

        print("Running bin merging sequence.")
        
        while len(bad_bin_indices) > 0:

            # get largest stat uncert superbin number
            bad_bin_num = bad_bin_nums[0]
            bad_bin = bad_bin_indices[0]

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
            bad_bin_indices, bad_bin_nums = self._get_bad_bins()
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
        return
    
    def _get_bad_bins(self) -> list:
    # return list of bad bin indices (those that are above the max_stat_uncert threshold), 
    # and a list of bad bin numbers (the position the bad bin occupies in the superbin_indices list)
    # step one: access current superbin list and check each bin's stat uncertainty. 
    # Note that one superbin can contain multiple bins. Must find total uncertainty for 
    # each superbin. A superbin containing multiple bins can be considered bad if the
    # total uncertainty is above the threshold. Bins are sorted by decreasing stat uncertainty.
        if not self.superbin_indices:
            print("No superbins initialized. Please run _init_superbin_indices() first.")
            return []
        bad_bin_indices = []
        bad_bin_numbers = []
        # Check each superbin for statistical uncertainty
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
            # calculate total statistical uncertainty for the superbin
            total_stat_uncert = np.sqrt(total_superbin_value) / total_superbin_value if total_superbin_value > 0 else 0
            if total_stat_uncert > self.max_stat_uncert:
                bad_bin_indices.append(superbin)
                bad_bin_numbers.append(self.superbin_indices.index(superbin))
        if self.debug:
            # In debug mode, assume hist is a numpy array
            bad_bins = []
            for index in bad_bin_indices:
                total_superbin_value = 0
                for bin_index in index:
                    bin_value = self.hist[bin_index[0], bin_index[1]]
                    total_superbin_value += bin_value
                total_stat_uncert = np.sqrt(total_superbin_value) / total_superbin_value if total_superbin_value > 0 else 0
                centroid_x = np.mean([bin_index[0] for bin_index in index])
                centroid_y = np.mean([bin_index[1] for bin_index in index])
                # Distance from center of histogram (x,y) value:
                center_x = self.hist.shape[0] // 2
                center_y = self.hist.shape[1] // 2
                # Calculate central distance from the center of the histogram
                central_distance = np.sqrt((centroid_x - center_x) ** 2 + (centroid_y - center_y) ** 2)
                # Append to bad_bins list
                bad_bins.append((total_stat_uncert, central_distance, centroid_x, centroid_y, index))
            # Sort by total_stat_uncert (descending), then central distance (closest to center of
            # histogram), then by centroid_x (highest), then by centroid_y (highest)
            bad_bins.sort(key=lambda x: (-x[0], x[1], -x[2], -x[3]))
            # Extract bad bin indices (a doubly nested list containing indices of bins in each superbin)
            bad_bin_indices = [bin_info[4] for bin_info in bad_bins]
            bad_bin_numbers = [self.superbin_indices.index(bin_info[4]) for bin_info in bad_bins]
            print(f"Found {len(bad_bin_indices)} bad bins exceeding the max_stat_uncert threshold.")
        else:
            for index in bad_bin_indices:
                superbin = self.superbin_indices[index]
                total_superbin_value = 0
                for bin_index in superbin:
                    bin_value = self.hist.GetBinContent(bin_index[0], bin_index[1])
                    total_superbin_value += bin_value
                total_stat_uncert = np.sqrt(total_superbin_value) / total_superbin_value if total_superbin_value > 0 else 0
                # Use get_superbin_centroids method to get centroid coordinates
                centroids_x, centroids_y = self.get_superbin_centroids()
                centroid_x = centroids_x[index]
                centroid_y = centroids_y[index]
                center_x = (self.hist.GetNbinsX() + 1) // 2  # ROOT bins are 1-based
                center_y = (self.hist.GetNbinsY() + 1) // 2
                central_distance = np.sqrt((centroid_x - center_x) ** 2 + (centroid_y - center_y) ** 2)
                bad_bins.append((total_stat_uncert, central_distance, centroid_x, centroid_y, superbin))
            # Sort by total_stat_uncert (descending), then by central_distance (ascending), then by centroid_x (descending), then by centroid_y
            bad_bins.sort(key=lambda x: (-x[0], x[1], -x[2], -x[3]))
            bad_bin_indices = [bin_info[4] for bin_info in bad_bins]
            bad_bin_numbers = [self.superbin_indices.index(bin_info[4]) for bin_info in bad_bins]
            print(f"Sorted bad bins by decreasing statistical uncertainty. Total bad bins: {len(bad_bins)}")
        if not bad_bin_indices:
            print("No bad bins found exceeding the max_stat_uncert threshold.")
            return [], [] 
        return bad_bin_indices, bad_bin_numbers  # return list of bad bins, each represented as a list of indices

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

    def get_superbin_centroids(self) -> list: # return list of superbin centroids
        centroids = []
        for superbin in self.superbin_indices:
            if self.debug:
                # In debug mode, assume hist is a numpy array
                x_coords = [bin_index[0] for bin_index in superbin]
                y_coords = [bin_index[1] for bin_index in superbin]
            else:
                # For ROOT histograms, use GetBinCenter
                x_coords = [self.hist.GetXaxis().GetBinCenter(bin_index[0]) for bin_index in superbin]
                y_coords = [self.hist.GetYaxis().GetBinCenter(bin_index[1]) for bin_index in superbin]
            centroid_x = np.mean(x_coords)
            centroid_y = np.mean(y_coords)
            centroids.append((centroid_x, centroid_y))
        return centroids

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