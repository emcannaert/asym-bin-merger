from asym_bin_merger import AsymBinMerger
import ROOT

# Open ROOT file and retrieve histogram
r_file = ROOT.TFile.Open("./run_examples/random_test_hist.root", "READ")
hist = r_file.Get("random")  # Ensure it's a TH2F

# Merge with max stat uncertainty of 10%
bin_merger = AsymBinMerger(hist, 0.12, "randomized_hist")

# Retrieve outputs
bin_merger.run()

# Open ROOT file and retrieve histogram
r_file = ROOT.TFile.Open("./run_examples/radial_test_hist.root", "READ")
hist = r_file.Get("radial")  # Ensure it's a TH2F

# Merge with max stat uncertainty of 10%
bin_merger = AsymBinMerger(hist, 0.14, "radial_hist")

# Retrieve outputs
bin_merger.run()

# Open ROOT file and retrieve histogram
r_file = ROOT.TFile.Open("./run_examples/quarter_radial_test_hist.root", "READ")
hist = r_file.Get("quarter_radial")  # Ensure it's a TH2F

# Merge with max stat uncertainty of 10%
bin_merger = AsymBinMerger(hist, 0.15, "quarter_radial_hist")

# Retrieve outputs
bin_merger.run()

