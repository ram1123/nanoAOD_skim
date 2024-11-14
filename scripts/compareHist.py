import uproot
import matplotlib.pyplot as plt
import numpy as np

# Open the ROOT file and access the TTree
file = uproot.open("skimmed_nano.root")
tree = file["Events"]

# Define the histogram bins
pt_bins = np.linspace(0, 150, 51)  # 50 bins from 0 to 150

# Load nMuon and Muon_pt branches
nMuon = tree["nMuon"].array(library="np")
muon_pt = tree["Muon_pt"].array(library="np")
muon_uncorrected_pt = tree["Muon_uncorrected_pt"].array(library="np")

# Apply the cut nMuon > 0 and select only Muon_pt[0] for events with at least one muon
muon_pt = muon_pt[nMuon > 0]
muon_uncorrected_pt = muon_uncorrected_pt[nMuon > 0]

# Extract the first muon pt (Muon_pt[0]) for each event
first_muon_pt = np.array([event[0] for event in muon_pt if len(event) > 0])
first_muon_uncorrected_pt = np.array([event[0] for event in muon_uncorrected_pt if len(event) > 0])

# Create histograms
counts1, bin_edges1 = np.histogram(first_muon_pt, bins=pt_bins)
counts2, bin_edges2 = np.histogram(first_muon_uncorrected_pt, bins=pt_bins)

# Plotting without displaying
fig, ax = plt.subplots()

# Plot histogram1 in blue
ax.step(bin_edges1[:-1], counts1, where="mid", label="Muon_pt[0]", color="blue")

# Plot histogram2 in red on the same canvas
ax.step(bin_edges2[:-1], counts2, where="mid", label="Muon_uncorrected_pt[0]", color="red")

# Customize the plot
ax.legend(title="First Muon Transverse Momentum")
ax.set_xlabel("Transverse Momentum (GeV)")
ax.set_ylabel("Counts")

# Save the plot to a PNG file
fig.savefig("histogram_comparison.png")
plt.close(fig)
