import ROOT
from collections import OrderedDict
import json

# Load the JSON data
cutFlowPath = 'CutFlowData.json'  # Make sure this path is correct
with open(cutFlowPath, 'r') as file:
    cutFlowData = json.load(file, object_pairs_hook=OrderedDict)

# Print the cut flow data
for cut, count in cutFlowData.items():
    print('{:<40} {:>10}'.format(cut, count))

# Create a histogram
nBins = len(cutFlowData)
histCutFlow = ROOT.TH1F("cutFlow", "Cut Flow", nBins, 0, nBins)

# Fill the histogram
for i, (cut, count) in enumerate(cutFlowData.items(), start=1):
    binNum = histCutFlow.FindBin(i)
    histCutFlow.SetBinContent(binNum, count)
    histCutFlow.GetXaxis().SetBinLabel(binNum, cut)

# Optional: Make labels vertical
histCutFlow.LabelsOption("v")
histCutFlow.GetXaxis().LabelsOption("di")
# histCutFlow.GetXaxis().SetLabelOffset(0.02)

# Draw the histogram
canvas = ROOT.TCanvas("canvas", "Cut Flow", 1200, 600)
histCutFlow.Draw("HIST TEXT0")
canvas.Draw()

# Save the canvas as a PNG file
canvas.SaveAs("cutFlowHistogram.png")
