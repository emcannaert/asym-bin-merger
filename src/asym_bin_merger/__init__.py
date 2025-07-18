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
Usage:
- Create an instance of AsymBinMerger.
- Call `.print_help()` for usage instructions.

"""

import os

import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class AsymBinMerger:
    def __init__(
        self, hist=None, max_stat_uncert=None, output_dir=None, debug: bool = False
    ):
        self.hist = hist
        self.max_stat_uncert = max_stat_uncert
        self.output_dir = output_dir
        self.debug = debug
        self.file_name = "bin_map.txt"
        self.output_file = os.path.join(self.output_dir, self.file_name)
        self.final_hist = np.zeros((5, 5))  # placeholder for final hist
        self.superbin_indices = []
        # check inputs
        self._validate_inputs()

        # run main workflow
        if not debug:
            print("Running AsymBinMerger in normal mode.")
            #self._run()
            #self._write_output()
        else:
            print("Running in debug mode.")

    ## Main methods
    def _convert_hist(self) -> list:
        """
        Convert the ROOT histogram to a 2D NumPy array for easier handling.
        Also sets self.final_hist.

        Returns:
            A 2D NumPy array.

        Raises:
            TypeError: If `hist` is not a NumPy array in debug mode
            or not a ROOT.TH2 histogram otherwise.
        """
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
            return self.final_hist
        else:

            raise TypeError(
                "`hist` must be a ROOT.TH2 object or a numpy array in debug mode."
            )

    def _init_superbin_indices(self) -> list:
        # Initialize superbin_indices s.t. every 2D bin is a superbin
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
            raise TypeError(
                "`hist` must be a ROOT.TH2 object or a numpy array in debug mode."
            )

        # If no bins are under the threshold, raise error
        if not superbins:
            raise ValueError(
                "Bins were not correctly initialized; check the input histogram."
            )

        self.superbin_indices = superbins
        return superbins

    def _check_superbins(self): # check whether there are problems with final superbins
        print("Checking bin map for issues.")

        # Check if superbin_indices is a list of lists
        if not all(isinstance(superbin, list) for superbin in self.superbin_indices):
            raise TypeError("superbin_indices must be a list of lists.")
        
        # Check if there are any empty superbins    
        # This is a simple check, but it can be extended to more complex checks if needed
        empty_superbins = [superbin for superbin in self.superbin_indices if not superbin]
        if empty_superbins:
            print(f"Warning: Found {len(empty_superbins)} empty superbins. These will be ignored in the merging process.")
            # Optionally, remove empty superbins from the list
            self.superbin_indices = [superbin for superbin in self.superbin_indices if superbin]
        
        # check if there are any issues with the superbin indices
        if not self.superbin_indices:
            raise Exception("No superbins initialized. Please run _init_superbin_indices() first.")
        
        # Check if all bins in superbin_indices are valid
        for superbin in self.superbin_indices:
            for bin_index in superbin:
                if self.debug:
                    # In debug mode, assume hist is a numpy array
                    if not (0 <= bin_index[0] < self.hist.shape[0] and 0 <= bin_index[1] < self.hist.shape[1]):
                        raise Exception(f"Invalid bin index {bin_index} in superbin {superbin}.")
                else:
                    # For ROOT histograms, check if the bin index is valid
                    if not (1 <= bin_index[0] <= self.hist.GetNbinsX() and 1 <= bin_index[1] <= self.hist.GetNbinsY()):
                        raise Exception(f"Invalid bin index {bin_index} in superbin {superbin}.")
                    
        print("Superbins check completed successfully. No issues found.")
        # If all checks pass, return successfully
        return

    def _run_bin_merging(self):  # run main bin merging scheme
        # init bad_bins
        bad_bin_indices, bad_bin_nums = self._get_bad_bins()
        it_num = 0

        print("Running bin merging sequence.")
        while len(bad_bin_indices) > 0:
            # get largest stat uncert superbin number
            bad_bin_num = bad_bin_nums[0]
            bad_bin = bad_bin_indices[0]

            # get list of neighbors 
            bad_neighbor_nums = self._get_neighbor_bins(bad_bin_num)

            if not bad_neighbor_nums:
                raise ValueError(
                    f"(Iteration {it_num}) No neighbor superbins \
                    found for superbin number {bad_bin_num} \
                    with superbins indices {self.superbin_indices}."
                )

            
            bad_neighbor_superbin = bad_neighbor_nums[0]  # This is a list of coordinates
            try:
                bad_neighbor_num = self.superbin_indices.index(bad_neighbor_superbin)
            except ValueError:
                # Fallback: find the superbin that contains all coordinates of bad_neighbor_superbin
                found = False
                for idx, superbin in enumerate(self.superbin_indices):
                    if all(coord in superbin for coord in bad_neighbor_superbin):
                        bad_neighbor_num = idx
                        found = True
                        break
                if not found:
                    print(f"Warning: Neighbor superbin {bad_neighbor_superbin} not found in current superbin_indices. Skipping this merge.")
                    it_num += 1
                    continue  # Skip this iteration
            self.superbin_indices[bad_neighbor_num].extend(self.superbin_indices[bad_bin_num])
            
            # remove bad bin from superbin indices
            self.superbin_indices.pop(bad_bin_num)

            # update bad bins list for next iteration
            bad_bin_indices, bad_bin_nums = self._get_bad_bins()
            it_num += 1
        print(f"Finished bin merging sequence - took {it_num} iterations.")

        return

    def _write_output(self): 
        """
        Writes bin map to self.output_file
        Raises:
            Exception: If no superbins are initialized or if the output directory does not exist.
        """
        # Check if output directory is set
        if self.output_dir is None:
            raise ValueError("Output directory is not set. Please provide a valid output directory.")
        
        # Check if superbin_indices is empty
        if not self.superbin_indices:
            raise Exception("No superbins to write. Please run _run_bin_merging() first.")
        # Check if output directory exists, create it if not
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Write the superbin indices to the output file
        with open(self.output_file, 'w') as f:
            f.write(str(self.superbin_indices))
        print("Writing output to %s."%(self.output_file))

        return

    def run(self):  # do main bin-merging scheme
        self.converted_hist = self._convert_hist()
        self.superbin_indices = self._init_superbin_indices()
        self._run_bin_merging()
        self._check_superbins()
        self._write_output()
        print("Bin merging completed successfully.")
        self._merged_hist_to_image()

    ## Helper functions

    def _print_help(self): # give information about usage
        """
        Print usage instructions for the AsymBinMerger utility.
        """
        print("=" * 60)
        print("AsymBinMerger: Optimize binning for a 2D ROOT histogram")
        print("=" * 60)
        print("Description:")
        print("  Merges bins in a 2D histogram such that the statistical uncertainty")
        print("  in each merged bin (superbin) remains below a given threshold.")
        print()
        print("Inputs:")
        print("  hist             : Required. A ROOT.TH2 histogram object")
        print("                     (or a 2D NumPy array in debug mode).")
        print("  max_stat_uncert  : Required. A float indicating the maximum allowed")
        print("                     relative statistical uncertainty per bin.")
        print("  output_dir       : Optional. Directory where output bin map is saved.")
        print("                     Defaults to current working directory if not provided.")
        print("  debug            : Optional. If True, skips main execution for testing (default: False).")
        print()
        print("Output:")
        print("  - A text file 'bin_map.txt' is written to output_dir.")
        print("  - It contains the final binning map in the form:")
        print("      [ [ (x1, y1), (x2, y2), ... ], [ ... ], ... ]")
        print("    where each sublist represents a merged 'superbin'.")
        print()
        print("Usage Example:")
        print("  merger = AsymBinMerger(hist=my_hist, max_stat_uncert=0.1,")
        print("                         output_dir='bin_maps', debug=False)")
        print("  bin_map = merger.get_bin_map()")
        print()
        print("Note:")
        print("  - Requires ROOT to be available in the runtime environment.")
        print("  - Use debug=True to test with NumPy arrays directly.")
        print("=" * 60)

    def _get_bad_bins(self) -> list:
        # return bad bin indices (those above the max_stat_uncert),
        # and bad bin numbers (index of bad bin in superbin_indices list)
        # step one: get stat uncertainty of all bins

        # Note: superbins can contain multiple bins.
        #  Must find total uncertainty for each superbin.
        # A superbin containing multiple bins can be considered bad if the
        # total uncertainty is above the threshold.
        # Bins are sorted by decreasing stat uncertainty.

        if not self.superbin_indices:
            print(
                "No superbins initialized. Please run _init_superbin_indices() first."
            )
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
                # Calculate Poissonian stat uncertainty
                total_superbin_value += bin_value
            # calculate total statistical uncertainty for the superbin
            total_stat_uncert = (
                np.sqrt(total_superbin_value) / total_superbin_value
                if total_superbin_value > 0
                else 0
            )
            if total_stat_uncert > self.max_stat_uncert:
                bad_bin_indices.append(superbin)
                bad_bin_numbers.append(self.superbin_indices.index(superbin))
        if self.debug:
            # In debug mode, assume hist is a numpy array
            bad_bins = []
            total_superbin_value = 0
            for index in bad_bin_numbers:
                superbin = self.superbin_indices[index]
                total_superbin_value = 0
                for bin_index in superbin:
                    bin_value = self.hist[bin_index[0], bin_index[1]]
                    total_superbin_value += bin_value
                total_stat_uncert = (
                    np.sqrt(total_superbin_value) / total_superbin_value
                    if total_superbin_value > 0
                    else 0
                )
                centroid_x = np.mean([bin_index[0] for bin_index in self.superbin_indices[index]])
                centroid_y = np.mean([bin_index[1] for bin_index in self.superbin_indices[index]])
                # Distance from center of histogram (x,y) value:
                center_x = self.hist.shape[0] // 2
                center_y = self.hist.shape[1] // 2
                # Calculate central distance from the center of the histogram
                central_distance = np.sqrt(
                    (centroid_x - center_x) ** 2 + (centroid_y - center_y) ** 2
                )
                # Append to bad_bins list
                bad_bins.append(
                    (total_stat_uncert, central_distance, centroid_x, centroid_y, superbin)
                )
            # Sort by total_stat_uncert (descending)
            # then central distance (closest to center of histogram),
            # then by centroid_x (highest),
            # then by centroid_y (highest)
            bad_bins.sort(key=lambda x: (-x[0], x[1], -x[2], -x[3]))
            # Extract bad bin indices (a doubly nested list
            # containing indices of bins in each superbin)
            bad_bin_indices = [bin_info[4] for bin_info in bad_bins]
            bad_bin_numbers = [
                self.superbin_indices.index(bin_info[4]) for bin_info in bad_bins
            ]
            print(
                f"Found {len(bad_bin_indices)} bad bins exceeding  \
                the max_stat_uncert threshold."
            )
        else:
            bad_bins = []
            for index in bad_bin_numbers:
                superbin = self.superbin_indices[index]
                total_superbin_value = 0
                for bin_index in superbin:
                    bin_value = self.hist.GetBinContent(bin_index[0], bin_index[1])
                    total_superbin_value += bin_value
                total_stat_uncert = (
                    np.sqrt(total_superbin_value) / total_superbin_value
                    if total_superbin_value > 0
                    else 0
                )
                # Use get_superbin_centroids method to get centroid coordinates
                centroids = self.get_superbin_centroids()
                centroids_x = [c[0] for c in centroids]
                centroids_y = [c[1] for c in centroids]
                centroid_x = centroids_x[index]
                centroid_y = centroids_y[index]
                center_x = (self.hist.GetNbinsX() + 1) // 2  # ROOT bins are 1-based
                center_y = (self.hist.GetNbinsY() + 1) // 2
                central_distance = np.sqrt(
                    (centroid_x - center_x) ** 2 + (centroid_y - center_y) ** 2
                )
                bad_bins.append(
                    (
                        total_stat_uncert,
                        central_distance,
                        centroid_x,
                        centroid_y,
                        superbin,
                    )
                )
            # Sort by total_stat_uncert (descending),
            # then by central_distance (ascending),
            # then by centroid_x (descending),
            # then by centroid_y
            bad_bins.sort(key=lambda x: (-x[0], x[1], -x[2], -x[3]))
            bad_bin_indices = [bin_info[4] for bin_info in bad_bins]
            bad_bin_numbers = [
                self.superbin_indices.index(bin_info[4]) for bin_info in bad_bins
            ]
            print(
                f"Sorted bad bins by decreasing stat uncertainty. Total bad bins: \
                    {len(bad_bins)}"
            )
        if not bad_bin_indices:
            print("No bad bins found exceeding the max_stat_uncert threshold.")
            return [], []
        return (
            bad_bin_indices,
            bad_bin_numbers,
        )  # return list of bad bins, each represented as a list of indices

    def _get_neighbor_bins(self, bad_bin_num) -> list:
        """
        Given an index into self.superbin_indices representing a "bad bin",
        return a list of neighboring superbins sorted by descending stat uncertainty,
        and secondarily by number of neighbors (also descending).
        """
        if not self.superbin_indices:
            print("No superbins.")
            return []

        superbins = self.superbin_indices
        bad_superbin = self.superbin_indices[bad_bin_num]


        neighbors = []
        neighbor_counts = []

        for superbin in superbins:
            count = 0
            if bad_superbin[0] in superbin:
                continue  # don't want to double count

            # loop over coordinates of each bin within superbin
            for coord in superbin:
                for bad_coord in bad_superbin:
                    if (abs(coord[0]-bad_coord[0])<=1 and abs(coord[1]-bad_coord[1])<=1):
                        if [coord] not in neighbors:
                            neighbors.append([coord])
                            count +=1

            if count > 0:
                neighbor_counts.append(count)

        # helper function for uncertainty to sort
        def stat_uncertainty(superbin):
            values = []
            for x, y in superbin:
                if self.debug:
                    values.append(self.hist[x, y])
                else:
                    values.append(self.hist.GetBinContent(x, y))
            total = sum(values)
            return np.sqrt(total) / total if total > 0 else 0

        # Pair each superbin with its sorting criteria:
        # (uncertainty, number of neighbors)
        criteria_with_bins = list(
            zip(
                [stat_uncertainty(sb) for sb in neighbors],
                [len(sb) for sb in neighbors],
                neighbors,
            )
        )

        # Sort based on uncertainty (descending), then number of neighbors (descending)
        criteria_with_bins.sort(key=lambda x: (-x[0], -x[1]))

        # Unpack the sorted bins
        neighbors = [x[2] for x in criteria_with_bins]

        if self.debug:
            print(f"Found {len(neighbors)} neighboring superbins.")

        return neighbors

    def print_bin_map(self):  # helper for testing and checking merging process
        print("Printing out formatted bin map")
        if not self.superbin_indices:
            print("No superbins to print. Please run _init_superbin_indices() first.")
            return
        for i, superbin in enumerate(self.superbin_indices):
            print(f"Superbin {i}:")
            for bin_index in superbin:
                if self.debug:
                    # In debug mode, assume hist is a numpy array
                    print(f"  Bin index: {bin_index} (value: {self.hist[bin_index[0], bin_index[1]]})")
                else:
                    # For ROOT histograms, use GetBinContent
                    value = self.hist.GetBinContent(bin_index[0], bin_index[1])
                    print(f"  Bin index: {bin_index} (value: {value})")
        print("End of bin map.") 
        return

    ## Getters and setters

    def get_bin_map(self) -> list: # external getter for superbin indices list 
        # Looks at the output file and returns the bin map as a list.
        if os.path.exists(self.output_file):
            with open(self.output_file, 'r') as f:
                bin_map = eval(f.read())
            return bin_map
        else:
            raise FileNotFoundError(f"Output file {self.output_file} does not exist. Please run _write_output() first.")

    def get_superbin_centroids(self) -> list:  # return list of superbin centroids
        centroids = []
        for superbin in self.superbin_indices:
            if self.debug:
                # In debug mode, assume hist is a numpy array
                x_coords = [bin_index[0] for bin_index in superbin]
                y_coords = [bin_index[1] for bin_index in superbin]
            else:
                # For ROOT histograms, use GetBinCenter
                x_coords = [
                    self.hist.GetXaxis().GetBinCenter(bin_index[0])
                    for bin_index in superbin
                ]
                y_coords = [
                    self.hist.GetYaxis().GetBinCenter(bin_index[1])
                    for bin_index in superbin
                ]
            centroid_x = np.mean(x_coords)
            centroid_y = np.mean(y_coords)
            centroids.append((centroid_x, centroid_y))
        return centroids

    ## testing and meta functions
    def _validate_inputs(self):
        if self.debug:
            if not isinstance(self.hist, np.ndarray):
                raise TypeError("In debug mode, `hist` must be a numpy array.")
            if not (
                (self.hist.ndim == 2) or
                (self.hist.ndim == 1 and self.hist.dtype == object and all(isinstance(row, np.ndarray) for row in self.hist))
            ):
                raise TypeError("In debug mode, `hist` must be a numpy array.")
        else:
            if not isinstance(self.hist, ROOT.TH2):
                raise TypeError("`hist` must be a ROOT.TH2 object.")
        if not isinstance(self.max_stat_uncert, float):
            raise TypeError("`max_stat_uncert` must be a float.")
        if not isinstance(self.output_dir, str):
            raise TypeError("`output_dir` must be a string.")

    def _get_merged_hist(self):  # testing: return np.array version of post-merge hist
        # from self.superbin_indices, get a new histogram in np.array format
        if not isinstance(self.hist, (np.ndarray, ROOT.TH2)):
            raise TypeError("`hist` must be a numpy array in debug mode or a ROOT.TH2 object.")
        if self.debug:
            if self.hist.ndim == 2:
                pass  # regular 2D array
            elif self.hist.ndim == 1 and self.hist.dtype == object and all(isinstance(row, np.ndarray) for row in self.hist):
                pass  # jagged array
            else:
                raise ValueError("In debug mode, hist must be a 2D numpy array or a jagged array (dtype=object, 1D array of arrays).")
                raise ValueError("In debug mode, `hist` cannot be an empty array.")
        else:
            if not isinstance(self.hist, ROOT.TH2):
                raise TypeError("`hist` must be a ROOT.TH2 object.")
            if self.hist.GetNbinsX() <= 0 or self.hist.GetNbinsY() <= 0:
                raise ValueError("`hist` must have positive number of bins in both dimensions.")
            if self.hist.GetEntries() == 0:
                raise ValueError("`hist` cannot be an empty ROOT histogram.")
            
        if not self.superbin_indices:
            print("No superbins initialized. Please run _init_superbin_indices() first.")
            return []
        merged_hist = np.zeros_like(self.final_hist)
        for superbin in self.superbin_indices:
            for bin_index in superbin:
                if self.debug:
                    # In debug mode, assume hist is a numpy array
                    if self.hist.ndim == 2:
                # Regular 2D array
                        if 0 <= bin_index[0] < merged_hist.shape[0] and 0 <= bin_index[1] < merged_hist.shape[1]:
                            merged_hist[bin_index[0], bin_index[1]] += self.hist[bin_index[0], bin_index[1]]
                    elif self.hist.ndim == 1 and self.hist.dtype == object:
                        # Jagged array: merged_hist and self.hist are 1D arrays of arrays
                        row, col = bin_index
                        if (
                            0 <= row < len(merged_hist)
                            and isinstance(merged_hist[row], np.ndarray)
                            and 0 <= col < len(merged_hist[row])
                        ):
                            merged_hist[row][col] += self.hist[row][col]
                    else:
                    # Skip invalid bin
                        continue
                else:
                    # For ROOT histograms, use GetBinContent
                    if 1 <= bin_index[0] <= self.hist.GetNbinsX() and 1 <= bin_index[1] <= self.hist.GetNbinsY():
                        merged_hist[bin_index[0]-1, bin_index[1]-1] += self.hist.GetBinContent(bin_index[0], bin_index[1])
                    else:
                        # Skip invalid bin
                        continue
        return merged_hist
    
    def _merged_hist_to_image(self): # testing: plot merged histogram as an image
        """
        Take the identified superbins, plot a 2D histogram plot
        with each superbin represented as a single bin and unique color
        """
        if not self.superbin_indices:
            print("No superbins initialized. Please run _init_superbin_indices() first.")
            return
        
        if self.debug:
            shape = self.hist.shape
        else:
            shape = (self.hist.GetNbinsX(), self.hist.GetNbinsY())
        label_map = np.full(shape, -1, dtype=int)  # -1 for unassigned bins
        #plot the merged histogram
        for superbin_idx, superbin in enumerate(self.superbin_indices):
            for bin_index in superbin:
                if self.debug:
                    x, y = bin_index
                    if 0 <= x < shape[0] and 0 <= y < shape[1]:
                        label_map[x, y] = superbin_idx
                else:
                    x, y = bin_index
                    if 1 <= x <= shape[0] and 1 <= y <= shape[1]:
                        label_map[x-1, y-1] = superbin_idx
        if label_map.size == 0:
            print("Merged histogram is empty. Cannot plot.")
            return
        
        all_bins = set((x, y) for x in range(shape[0]) for y in range(shape[1]))
        assigned_bins = set()
        for superbin in self.superbin_indices:
            for bin_index in superbin:
                if self.debug:
                    assigned_bins.add(bin_index)
                else:
                    assigned_bins.add((bin_index[0]-1, bin_index[1]-1))
        missing_bins = all_bins - assigned_bins
        if missing_bins:
            print(f"Warning: {len(missing_bins)} bins are not assigned to any superbin: {missing_bins}")

        from matplotlib.colors import ListedColormap
        import matplotlib.patches as patches

        n_superbins = len(self.superbin_indices)
        base_cmap = plt.get_cmap('tab20')
        colors = base_cmap.colors if n_superbins <= 20 else plt.cm.get_cmap('tab20b', n_superbins).colors
        cmap = ListedColormap(colors[:n_superbins])

        plt.figure(figsize=(10, 8))
        plt.imshow(label_map, cmap=cmap, interpolation='nearest', origin='lower',)
        plt.title('Merged Histogram with Superbins')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.xticks(ticks=np.arange(shape[1]), labels=np.arange(shape[1]))
        plt.yticks(ticks=np.arange(shape[0]), labels=np.arange(shape[0]))
        plt.grid(False)

        ax = plt.gca()
        for superbin_idx, superbin in enumerate(self.superbin_indices):
            # Get all bin coordinates for this superbin
            coords = [(b[0], b[1]) if self.debug else (b[0]-1, b[1]-1) for b in superbin]
            xs = [c[1] for c in coords]
            ys = [c[0] for c in coords]
            # Draw a rectangle around the bounding box of the superbin
            #min_x, max_x = min(xs), max(xs)
            #min_y, max_y = min(ys), max(ys)
            #rect = patches.Rectangle((min_x-0.5, min_y-0.5), max_x-min_x+1, max_y-min_y+1,
            #                        linewidth=2, edgecolor='black', facecolor='none')
            #ax.add_patch(rect)
            # Annotate centroid with stat uncertainty
            centroid_x = int(np.mean(xs))
            centroid_y = int(np.mean(ys))
            if self.debug:
                total_value = sum(self.hist[y, x] for y, x in coords)
            else:
                total_value = sum(self.hist.GetBinContent(y+1, x+1) for y, x in coords)
            stat_uncert = (np.sqrt(total_value) / total_value) if total_value > 0 else 0
            plt.text(centroid_x, centroid_y, f"{stat_uncert:.2f}", color='black',
                    ha='center', va='center', fontsize=10, fontweight='bold')


        plt.savefig(os.path.join(self.output_dir, 'merged_histogram.png'))
        plt.show()
        plt.close()
        print("Merged histogram plotted successfully.")
        # Save the plot
        print("Merged histogram saved as 'merged_histogram.png' in the output directory.")
        return       


    