#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <iostream>

void normalizeAndStitchDYsamples() {
    double luminosity = 60000; // pb^-1

    // Cross-sections in pb
    double crossSections[] = {1485.0, 397.4, 97.2, 3.701, 0.5086, 0.04728};

    // File paths
    const char* filePaths[] = {
        "/eos/user/a/avijay/DYsamples/dy_0to50/dy_0to50.root",
        // "/eos/user/a/avijay/DYsamples/dy_50to100/dy_50to100.root",
        // "/eos/user/a/avijay/DYsamples/dy_100to250/dy_100to250.root",
        // "/eos/user/a/avijay/DYsamples/dy_250to400/dy_250to400.root",
        // "/eos/user/a/avijay/DYsamples/dy_400to650/dy_400to650.root",
        // "/eos/user/a/avijay/DYsamples/dy_650toInf/dy_650toInf.root"
    };

    // Open output file
    TFile outputFile("normalizedDYsamples.root", "RECREATE");

    for (int i = 0; i < 6; i++) {
        // Open each file
        TFile* file = TFile::Open(filePaths[i]);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filePaths[i] << std::endl;
            continue;
        }

        // Access the "Runs" tree and the "genEventCount" branch
        TTree* tree;
        file->GetObject("Runs", tree);
        if (!tree) {
            std::cerr << "Runs tree not found in file: " << filePaths[i] << std::endl;
            file->Close();
            continue;
        }

        Long64_t genEventCount, totalEvents = 0;
        tree->SetBranchAddress("genEventCount", &genEventCount);

        for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
            tree->GetEntry(j);
            totalEvents += genEventCount;
        }

        // Calculate normalization factor based on total events, cross-section, and luminosity
        double normalizationFactor = (crossSections[i] * luminosity) / totalEvents;

        std::cout << "Total events in file " << filePaths[i] << ": " << totalEvents << "\tNormalization factor: " << normalizationFactor << std::endl;

        // Now, access the "Events" tree to fill histograms for pTZ1
        TTree *eventsTree;
        file->GetObject("Events", eventsTree);
        if (!eventsTree) {
            std::cerr << "Events tree not found in file: " << filePaths[i] << std::endl;
            file->Close();
            continue;
        }

        Float_t pTZ1;
        Bool_t passZZ2l2qSelection;
        eventsTree->SetBranchAddress("pTZ1", &pTZ1);
        eventsTree->SetBranchAddress("passZZ2l2qSelection", &passZZ2l2qSelection);

        TCanvas *c = new TCanvas("c", "c", 800, 600);
        // Create histogram for pTZ1
        TH1F *hist = new TH1F(Form("hist_pTZ1_%d", i), "pTZ1 distribution; pTZ1 (GeV); Events", 100, 0, 2000);

        Long64_t nEntries = eventsTree->GetEntries();
        for (Long64_t k = 0; k < 100; ++k)
        {
            eventsTree->GetEntry(k);
            std::cout << "passZZ2l2qSelection: " << passZZ2l2qSelection << std::endl;
            if (!passZZ2l2qSelection)
                continue;
            // std::cout << "pTZ1: " << pTZ1 << std::endl;
            hist->Fill(pTZ1);
        }
        // Normalize the histogram
        hist->Scale(normalizationFactor);
        c->SaveAs(Form("pTZ1_%d.png", i));

        // Write the histogram to the output file

        hist->Write();

        // Close the input file
        file->Close();
    }
    // Close the output file
    outputFile.Close();
}
