#include <iostream>
#include <set>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

void removeDuplicates()
{

    TString prefix = "Muon";
    TString outputFile = "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/Data_noDuplicates.root";

        TString infiles[] = {"/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/DoubleMuon.root",
                             "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/SingleMuon.root",
                             "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/MuonEG.root",
                             "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/Muon.root",
                             "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/EGamma.root"};

    TChain *oldtree = new TChain("Events");
    for (TString infile : infiles)
    {
        oldtree->Add(infile);
    }

    // Run from Runs tree

    Long64_t nentries = oldtree->GetEntries();
    std::cout << nentries << " total entries." << std::endl;
    ULong64_t Run, LumiSect, Event;
    bool passedZ4lSelection;
    oldtree->SetBranchAddress("run", &Run);
    oldtree->SetBranchAddress("luminosityBlock", &LumiSect);
    oldtree->SetBranchAddress("event", &Event);

    // Create a new file + a clone of old tree in new file
    // TFile *newfile = new TFile(
    //     prefix + "_noDuplicates.root", "recreate");
    TFile *newfile = new TFile(outputFile, "recreate");
    TTree *newtree = oldtree->CloneTree(0);

    std::set<TString> runlumieventSet;
    int nremoved = 0;
    for (Long64_t i = 0; i < nentries; i++)
    {
        if (i % 100000 == 0)
            std::cout << i << "/" << nentries << std::endl;
        oldtree->GetEntry(i);

        TString s_Run = std::to_string(Run);
        TString s_Lumi = std::to_string(LumiSect);
        TString s_Event = std::to_string(Event);
        TString runlumievent = s_Run + ":" + s_Lumi + ":" + s_Event;

        if (runlumieventSet.find(runlumievent) == runlumieventSet.end())
        {
            runlumieventSet.insert(runlumievent);
            newtree->Fill();
        }
        else
        {
            nremoved++;
        }
        // if (passedZ4lSelection) newtree->Fill();
    }

    std::cout << nremoved << " duplicates." << std::endl;
    newtree->Print();
    newtree->AutoSave();
    // delete oldfile;
    delete newfile;
}
