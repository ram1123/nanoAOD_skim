#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>

void stitchDYsamples() {
    // Create a TChain to hold all the trees
    TChain chain("Events");

    // Add the ROOT files to the chain
    chain.Add("/eos/user/a/avijay/DYsamples/dy_0to50/dy_0to50.root");
    chain.Add("/eos/user/a/avijay/DYsamples/dy_50to100/dy_50to100.root");
    chain.Add("/eos/user/a/avijay/DYsamples/dy_100to250/dy_100to250.root");
    chain.Add("/eos/user/a/avijay/DYsamples/dy_250to400/dy_250to400.root");
    chain.Add("/eos/user/a/avijay/DYsamples/dy_400to650/dy_400to650.root");
    chain.Add("/eos/user/a/avijay/DYsamples/dy_650toInf/dy_650toInf.root");

    // Variable to hold the pTZ1 value
    Float_t pTZ1;
    // Set the branch address
    chain.SetBranchAddress("pTZ1", &pTZ1);

    // Define histograms for different pT bins
    TH1F *h0to50 = new TH1F("h0to50", "DY pT: 0 to 50 GeV", 100, 0, 50);
    TH1F *h50to100 = new TH1F("h50to100", "DY pT: 50 to 100 GeV", 50, 50, 100);
    TH1F *h100to250 = new TH1F("h100to250", "DY pT: 100 to 250 GeV", 150, 100, 250);
    TH1F *h250to400 = new TH1F("h250to400", "DY pT: 250 to 400 GeV", 150, 250, 400);
    TH1F *h400to650 = new TH1F("h400to650", "DY pT: 400 to 650 GeV", 250, 400, 650);
    TH1F *h650toInf = new TH1F("h650toInf", "DY pT: 650 to Inf GeV", 350, 650, 1000); // Adjust max bin as needed

    // Assign different colors for clarity
    h0to50->SetLineColor(kBlue);
    h50to100->SetLineColor(kGreen);
    h100to250->SetLineColor(kRed);
    h250to400->SetLineColor(kCyan);
    h400to650->SetLineColor(kMagenta);
    h650toInf->SetLineColor(kYellow);

    // Loop over all entries
    Long64_t nentries = chain.GetEntries();
    for (Long64_t i=0; i<nentries; i++) {
        chain.GetEntry(i);

        // Fill the appropriate histogram based on pTZ1 value
        if (pTZ1 <= 50) h0to50->Fill(pTZ1);
        else if (pTZ1 <= 100) h50to100->Fill(pTZ1);
        else if (pTZ1 <= 250) h100to250->Fill(pTZ1);
        else if (pTZ1 <= 400) h250to400->Fill(pTZ1);
        else if (pTZ1 <= 650) h400to650->Fill(pTZ1);
        else h650toInf->Fill(pTZ1);
    }

    // Draw histograms on a canvas
    TCanvas *c = new TCanvas("c", "DY Stitched Samples", 800, 600);
    h0to50->Draw();
    h50to100->Draw("same");
    h100to250->Draw("same");
    h250to400->Draw("same");
    h400to650->Draw("same");
    h650toInf->Draw("same");

    // Optionally, add legend, axis labels, etc. here

    c->BuildLegend();
    c->Update();
    c->SaveAs("DY_Stitched_Samples.png"); // Save the canvas to a file
}

int main() {
    stitchDYsamples();
    return 0;
}

