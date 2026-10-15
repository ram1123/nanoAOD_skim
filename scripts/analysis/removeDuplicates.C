// Remove duplicate data events after skimming, when combining several primary datasets
// (e.g. DoubleMuon, SingleMuon, EGamma, MuonEG): the same collision can be stored in more
// than one of them. An event is identified by (run, luminosityBlock, event), and the FIRST
// occurrence is kept, so list the input files in dataset-priority order.
//
// Usage:
//   root -l -b -q 'scripts/analysis/removeDuplicates.C+("inputs.txt", "Data_noDuplicates.root")'
//
// inputs.txt: one ROOT file (or a pattern accepted by TChain::Add, e.g. dir/*.root) per line;
// empty lines and lines starting with '#' are ignored.
// Only the Events tree is written; MC needs no duplicate removal.

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

void removeDuplicates(const char *inputList, const char *outputFile, const char *treeName = "Events")
{
    std::vector<TString> inputs;
    std::ifstream list(inputList);
    if (!list) {
        std::cerr << "ERROR: cannot open input list " << inputList << std::endl;
        return;
    }
    std::string line;
    while (std::getline(list, line)) {
        TString entry(line);
        entry = entry.Strip(TString::kBoth);
        if (entry.IsNull() || entry.BeginsWith("#"))
            continue;
        inputs.push_back(entry);
    }
    if (inputs.empty()) {
        std::cerr << "ERROR: no input files in " << inputList << std::endl;
        return;
    }

    TChain chain(treeName);
    std::vector<int> firstTreeOfInput;  // chain tree number where each input line starts
    for (const auto &in : inputs) {
        firstTreeOfInput.push_back(chain.GetNtrees());
        if (chain.Add(in) == 0)
            std::cerr << "WARNING: no files added for " << in << std::endl;
    }

    // NanoAOD types: run and luminosityBlock are UInt_t, event is ULong64_t
    UInt_t run = 0, lumi = 0;
    ULong64_t event = 0;
    chain.SetBranchAddress("run", &run);
    chain.SetBranchAddress("luminosityBlock", &lumi);
    chain.SetBranchAddress("event", &event);

    TFile out(outputFile, "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "ERROR: cannot create " << outputFile << std::endl;
        return;
    }
    TTree *outTree = chain.CloneTree(0);
    if (!outTree) {
        std::cerr << "ERROR: no clonable tree named " << treeName << " in the input files" << std::endl;
        return;
    }

    std::set<std::tuple<UInt_t, UInt_t, ULong64_t>> seen;
    std::vector<Long64_t> nRead(inputs.size(), 0), nDup(inputs.size(), 0);
    const Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << " entries in " << chain.GetNtrees() << " files" << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        if (chain.GetEntry(i) <= 0)
            continue;
        if (i % 1000000 == 0)
            std::cout << "  " << i << " / " << nEntries << std::endl;

        // which input line this entry came from (for the per-input summary)
        size_t idx = 0;
        while (idx + 1 < firstTreeOfInput.size() && chain.GetTreeNumber() >= firstTreeOfInput[idx + 1])
            idx++;
        nRead[idx]++;

        if (seen.insert(std::make_tuple(run, lumi, event)).second)
            outTree->Fill();
        else
            nDup[idx]++;
    }

    out.cd();
    outTree->Write();
    const Long64_t nKept = outTree->GetEntries();
    out.Close();

    std::cout << "\nPer input (first occurrence kept; order = priority):" << std::endl;
    Long64_t totDup = 0;
    for (size_t k = 0; k < inputs.size(); k++) {
        std::cout << "  " << inputs[k] << ": read " << nRead[k] << ", duplicates removed " << nDup[k] << std::endl;
        totDup += nDup[k];
    }
    std::cout << "Kept " << nKept << " events, removed " << totDup << " duplicates -> " << outputFile << std::endl;
}
