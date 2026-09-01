/////-----------------/////
// script to make the input root file for the limit extraction using Higgs combine tool
//////----------------/////

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <iostream>
#include <vector>

using namespace std;

struct Sample {
    TString name;
    TString filepath;
    double xsec;
    double nMC;
    double nNeg;
};

//========================================================

int make_shapes() {

    // Output file
    TFile *fout = TFile::Open("input_root_file_combine_v4.root", "UPDATE");
    if (!fout || fout->IsZombie()) {
        cout << "Cannot open output file" << endl;
        return 1;
    }

    // Luminosity from data
    double lumi = 59547.0; // pb^-1

    //---------------------------------
    //define backgrounds
    //---------------------------------
    vector<Sample> samples = {
        //--dibosons---
        {"GluGluToContinToZZTo2e2nu",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToContinToZZTo2e2nu.root",
         0.01719, 500000.0, 1.0},

        {"GluGluToContinToZZTo2mu2nu",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToContinToZZTo2mu2nu.root",
         0.01719, 500000.0, 1.0},

        {"GluGluToWWToENEN",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToWWToENEN.root",
         0.0368, 4904000.0, 1.0},

        {"GluGluToWWToENMN",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToWWToENMN.root",
         0.0368, 4928000.0, 1.0},

        {"GluGluToWWToMNMN",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToWWToMNMN.root",
         0.0368, 4986000.0, 1.0},

        {"WWTo2L2Nu",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WWTo2L2Nu.root",
            11.09, 9994000.0, 1.0},

        {"WZTo2Q2L",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WZTo2Q2L.root",
            6.419, 28576996.0, 1.0},

        {"WZTo3LNu",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WZTo3LNu.root",
            5.257, 9821283.0, 1.0},

        {"ZZTo2L2Nu",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ZZTo2L2Nu.root",
            0.9738, 56886000.0, 1.0},

        {"ZZTo2Q2L",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ZZTo2Q2L.root",
            3.676, 29357938.0, 1.0},
        
        
        // -- top --
        {"TTTo2L2Nu",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/TTTo2L2Nu.root",
         88.29, 145020000.0, 1.0},

        {"ST_s-channel_4f_leptonDecays",
        "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_s-channel_4f_leptonDecays.root",
         3.549, 19365999.0, 1.0},

        {"ST_t-channel_top_4f_InclusiveDecays",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_t-channel_top_4f_InclusiveDecays.root",
         119.7, 178336000.0, 1.0},

         {"ST_t-channel_antitop_4f_InclusiveDecays",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_t-channel_antitop_4f_InclusiveDecays.root",
            67.93, 95627000.0, 1.0},

        {"ST_tW_top_5f_NoFullyHadronicDecays",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_tW_top_5f_NoFullyHadronicDecays.root",
            34.91, 11270430.0, 1.0},

        {"ST_tW_antitop_5f_NoFullyHadronicDecays",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_tW_antitop_5f_NoFullyHadronicDecays.root",
            34.91, 10949620.0, 1.0},

        {"tZq_ll_4f_ckm_NLO",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/tZq_ll_4f_ckm_NLO.root",
            0.07561, 11916000.0, 1.0},

        {"TTWJetsToLNu",
         "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/TTWJetsToLNu.root",
         0.2163, 10450000.0, 1.0},

        {"TTZToLLNuNu",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/TTZToLLNuNu.root",
            0.2439, 19608000.0, 1.0}, 

        // -- triboson ----
        {"WWZ_4F",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WWZ_4F.root",
            0.1707, 10209999.0, 1.0},
            
        {"ZZZ",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ZZZ.root",
            0.01476, 10139000.0, 1.0}, 

        {"WZZ",
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WZZ.root",
            0.05709, 10294000.0, 1.0},

    };

    // Common weight
    TString weight = "puWeight";

    // Variable
    TString var = "HZZ2l2nu_ZZmT";

    // Histogram settings
    int nbins = 100;
    float xmin = 0;
    float xmax = 1000;

    //========================================================
    // LOOP OVER SAMPLES
    //========================================================
    for (auto &s : samples) {

        cout << "Processing: " << s.name << endl;

        TFile *fin = TFile::Open(s.filepath);
        if (!fin || fin->IsZombie()) {
            cout << "Cannot open file: " << s.filepath << endl;
            continue;
        }

        TTree *tree = (TTree*)fin->Get("Events");
        if (!tree) {
            cout << "Missing tree in: " << s.filepath << endl;
            fin->Close();
            continue;
        }

        // ---- Scale factor ----
        double scale = (s.xsec * lumi) / (s.nMC - 2.0 * s.nNeg);

        // ---- Histogram names ----
        //TString hname     = s.name + "_alpha";
        TString hname     = s.name;
        TString hnameUp   = s.name + "_alphaUp";
        TString hnameDown = s.name + "_alphaDown";

        // ---- Create histogram ----
        TH1F *h = new TH1F(hname, ";m_{T};Events", nbins, xmin, xmax);
        h->Sumw2();

        // ---- Fill with weight ----
        tree->Draw(var + " >> " + hname, weight, "goff");

        // ---- Handle overflow ----
        int N = h->GetNbinsX();
        h->SetBinContent(N, h->GetBinContent(N) + h->GetBinContent(N+1));

        // ---- Apply normalization ----
        h->Scale(scale);

        // ---- Create systematics (for now identical) ----
        TH1* hUp   = (TH1*)h->Clone(hnameUp);
        TH1* hDown = (TH1*)h->Clone(hnameDown);

        // ---- Write to file ----
        fout->cd();
        h->Write("", TObject::kOverwrite);
        hUp->Write("", TObject::kOverwrite);
        hDown->Write("", TObject::kOverwrite);

        // Cleanup
        fin->Close();

        cout << "Done: " << s.name << endl;
    }

    //------------------
    //add data
    //-------------------
    TFile *fdata = TFile::Open(
        "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/data/data_2018_noDuplicates.root"
    );

    if (!fdata || fdata->IsZombie()) {
        cout << "Cannot open data file" << endl;
        return 1;
    }

    TTree *tdata = (TTree*)fdata->Get("Events");

    if (!tdata) {
        cout << "Missing tree" << endl;
        return 1;
    }

    TH1F *h_data = new TH1F("data_obs",";m_{T};Events",100,0,1000);
    tdata->Draw("HZZ2l2nu_ZZmT >> data_obs", "", "goff");

    // write data hist
    fout->cd();
    h_data->Write("data_obs", TObject::kOverwrite);

    //----------------
    //add signal
    //----------------
    cout << "Processing signal: ggH125" << endl;

    TFile *fsig = TFile::Open("/eos/user/a/avijay/Mergedrootfiles_hzz_new/v15/signal_v15/signal_m125.root");

    if (!fsig || fsig->IsZombie()) {
        cout << "Cannot open signal file" << endl;
        return 1;
    }

    TTree *tsig = (TTree*)fsig->Get("Events");

    if (!tsig) {
        cout << "Missing signal tree" << endl;
        return 1;
    }

    // ---- signal inputs ----
    double sig_xsec = 39.71;    
    double sig_ngen = 486250.0;  

    double sig_scale = (sig_xsec * lumi) / sig_ngen;

    // ---- histograms ----
    TString hsig     = "signal_alpha";
    TH1F *h_sig = new TH1F(hsig, ";m_{T};Events", nbins, xmin, xmax);
    h_sig->Sumw2();

    tsig->Draw(var + " >> signal_alpha", weight, "goff");
    // ---- overflow ----
    int N = h_sig->GetNbinsX();
    h_sig->SetBinContent(N, h_sig->GetBinContent(N) + h_sig->GetBinContent(N+1));

    // ---- scale ----
    h_sig->Scale(sig_scale);
    fout->cd();
    h_sig->Write("signal_alpha", TObject::kOverwrite);


    // ===================== drellyan=====================

cout << "Processing DY inclusive..." << endl;

// DY inputs
vector<tuple<TString, double, double>> dy_files = {
    {"/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/drellyan_ptz_0to50.root",   1485,   196207761.0},
    {"/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/drellyan_ptz_50to100.root", 397.4,  122967660.0},
    {"/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/drellyan_ptz_100to250.root",97.2,   79527324.0},
    {"/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/drellyan_ptz_250to400.root",3.701,  24195330.0},
    {"/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/drellyan_ptz_400to650.root",0.5086, 3936102.0},
    {"/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/drellyan_ptz_650toInf.root",0.04728,3994997.0}
};

// Create inclusive histograms
TH1F *h_dy = new TH1F("dy_alpha", ";m_{T};Events", nbins, xmin, xmax);
h_dy->Sumw2();

for (auto &finfo : dy_files) {

    TString fname;
    double xsec, ngen;
    tie(fname, xsec, ngen) = finfo;

    TString fullpath = "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/" + fname;

    cout << "Processing DY file: " << fullpath << endl;

    TFile *f = TFile::Open(fullpath);
    if (!f || f->IsZombie()) continue;

    TTree *tree = (TTree*)f->Get("Events");
    if (!tree) continue;

    double scale = (xsec * lumi) / ngen;

    TH1F *h_tmp = new TH1F("h_tmp","",nbins,xmin,xmax);
    h_tmp->Sumw2();

    tree->Draw(var + ">>h_tmp", weight, "goff");

    h_tmp->Scale(scale);
    h_dy->Add(h_tmp);

    delete h_tmp;
    f->Close();
}

// overflow
int Ndy = h_dy->GetNbinsX();
h_dy->SetBinContent(Ndy, h_dy->GetBinContent(Ndy) + h_dy->GetBinContent(Ndy+1));

// clones
TH1* h_dyUp   = (TH1*)h_dy->Clone("dy_alphaUp");
TH1* h_dyDown = (TH1*)h_dy->Clone("dy_alphaDown");

// rename nominal properly
h_dy->SetName("dy_alpha");

// write
fout->cd();
h_dy->Write("", TObject::kOverwrite);
h_dyUp->Write("", TObject::kOverwrite);
h_dyDown->Write("", TObject::kOverwrite);

cout << "DY done!" << endl;

    fout->Close();

    cout << "\n All backgrounds processed successfully!\n";

    return 0;
}