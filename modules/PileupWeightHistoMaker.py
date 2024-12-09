import ROOT
import yaml
from data.Run3_2022_LHC_Simulation_10h_2h_cfi import *

def getPUweight(pileupdatafile, inputhistoname, outputhistoname, outdir):

    input_file = ROOT.TFile.Open(pileupdatafile)
    input_hist = input_file.Get(inputhistoname)
    output_hist = ROOT.TH1D(outputhistoname, outputhistoname, input_hist.GetNbinsX(), input_hist.GetXaxis().GetXmin(), input_hist.GetXaxis().GetXmax())
    sum_data_entries = 0
    sum_MC_entries = 0
    for j in range(1, input_hist.GetNbinsX() + 1):
        sum_data_entries += input_hist.GetBinContent(j)
        sum_MC_entries += probValue[j-1]
    for i in range(1, input_hist.GetNbinsX() + 1):
        bin_content = input_hist.GetBinContent(i)/(probValue[i-1]*sum_data_entries/sum_MC_entries)
        output_hist.SetBinContent(i, bin_content)
    output_file = ROOT.TFile.Open(outdir, "RECREATE")
    output_hist.Write()
    input_file.Close()
    output_file.Close()

def MakeHisto(year):

    if year == 2022:
        cfgFile = "Input_2022.yml"    
    with open(cfgFile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        inputdataNPV = cfg['inputdataNPV']
        inputdataNPVup = cfg['inputdataNPVup']
        inputdataNPVdown = cfg['inputdataNPVdown']
        outputdataNPV = cfg['outputdataNPV']
        outputdataNPVup = cfg['outputdataNPVup']
        outputdataNPVdown = cfg['outputdataNPVdown']
        outputhistoname = cfg['PUweightHistoName']

    getPUweight(inputdataNPV, "pileup", outputhistoname, outputdataNPV)
    getPUweight(inputdataNPVup, "pileup", outputhistoname, outputdataNPVup)
    getPUweight(inputdataNPVdown, "pileup", outputhistoname, outputdataNPVdown)

MakeHisto(2022)

