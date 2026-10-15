#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1.h>
#include <TString.h>
#include <fstream>
#include <iomanip>

#include <iostream>
#include <vector>

using namespace std;

//====================================================
// Structures
//====================================================

struct Sample {

    TString name;
    vector<TString> files;

    double xsec;
    double ngen;

    TString category;
};

struct WJetBin {

    vector<TString> files;

    double xsec;

    double ngen_total;
};

//====================================================
// Helper function
//====================================================

void add_sample_to_histogram(
    TH1F* h_target,
    Sample &s,
    TString var,
    TString totalWeight,
    double lumi,
    int nbins,
    float xmin,
    float xmax
) {

    cout << "   -> " << s.name << endl;

    //-----------------------------------------
    // Temporary histogram
    //-----------------------------------------
    TH1F *h_tmp_total = new TH1F(
        "h_tmp_total",
        "",
        nbins,
        xmin,
        xmax
    );

    h_tmp_total->Sumw2();

    //-----------------------------------------
    // Loop over files
    //-----------------------------------------
    for (auto &fname : s.files) {

        TFile *f = TFile::Open(fname);

        if (!f || f->IsZombie()) {

            cout << "Cannot open: "
                 << fname << endl;

            continue;
        }

        TTree *tree =
            (TTree*)f->Get("Events");

        if (!tree) {

            cout << "Missing tree in "
                 << fname << endl;

            f->Close();
            continue;
        }

        TH1F *h_tmp = new TH1F(
            "h_tmp",
            "",
            nbins,
            xmin,
            xmax
        );

        h_tmp->Sumw2();

        tree->Draw(
            var + ">>h_tmp",
            totalWeight,
            "goff"
        );

        h_tmp_total->Add(h_tmp);

        delete h_tmp;

        f->Close();
    }

    //-----------------------------------------
    // Overflow
    //-----------------------------------------
    int N = h_tmp_total->GetNbinsX();

    h_tmp_total->SetBinContent(
        N,
        h_tmp_total->GetBinContent(N)
        + h_tmp_total->GetBinContent(N+1)
    );

    //-----------------------------------------
    // Normalize
    //-----------------------------------------
    double scale =
        (s.xsec * lumi) / s.ngen;

    h_tmp_total->Scale(scale);

    //-----------------------------------------
    // Add to category histogram
    //-----------------------------------------
    h_target->Add(h_tmp_total);

    delete h_tmp_total;
}

//====================================================
// Main
//====================================================

int make_shapes_final() {

    //------------------------------------------------
    // Output file
    //------------------------------------------------
    TFile *fout = TFile::Open(
        "input_root_file_combine_grouped_with_signal_selection.root",
        "RECREATE"
    );

    if (!fout || fout->IsZombie()) {

        cout << "Cannot open output file"
             << endl;

        return 1;
    }

    ofstream yieldFile("yields_with_signal_selection.txt");

    yieldFile << fixed << setprecision(3);

    yieldFile << "=============================\n";
    yieldFile << "Process           Yield\n";
    yieldFile << "=============================\n";

    //------------------------------------------------
    // Settings
    //------------------------------------------------
    double lumi = 59547.0;

    TString weight = "puWeight";
    TString var    = "HZZ2l2nu_ZZmT";

//------------------------------------------------
// Signal Region Selection
//------------------------------------------------

TString SRcut =

    "(HZZ2l2qNu_isELE==1)"                    // electron channel

    " && ((pTL1>25)&&(pTL2>25))"              // lepton pt

    " && (HZZ2l2qNu_cutOppositeChargeFlag==0)" // opposite charge

    " && (passZZ2l2nuSelection==1)"           // 2l2nu selection

    " && (HZZ2l2qNu_nMediumBtagJets==0)"      // no b-tag

    " && (pT_MET>100)"                        // MET

    " && (massZ1>76&&massZ1<106)"             // Z mass window

    " && (pTZ1>55)"                           // Z boson pt

    " && (HZZ2l2nu_ifVBF==0)"                 // non-VBF

    " && (HZZ2l2qNu_nJets==0)"                // zero jet

    " && (Triggers_HZZ2l2nu_SingleLep==1)";  // trigger

    int nbins  = 100;
    float xmin = 0;
    float xmax = 1000;

    TString totalWeight = weight + " * (" + SRcut + ")";

    //------------------------------------------------
    // Define ALL samples
    //------------------------------------------------

    vector<Sample> samples = {

        //================================================
        // DIBOSON
        //================================================

        {
            "GluGluToContinToZZTo2e2nu",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToContinToZZTo2e2nu.root"
            },
            0.01719,
            500000.0,
            "diboson"
        },

        {
            "GluGluToContinToZZTo2mu2nu",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToContinToZZTo2mu2nu.root"
            },
            0.01719,
            500000.0,
            "diboson"
        },

        {
            "GluGluToWWToENEN",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToWWToENEN.root"
            },
            0.0368,
            4904000.0,
            "diboson"
        },

        {
            "GluGluToWWToENMN",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToWWToENMN.root"
            },
            0.0368,
            4928000.0,
            "diboson"
        },

        {
            "GluGluToWWToMNMN",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/GluGluToWWToMNMN.root"
            },
            0.0368,
            4986000.0,
            "diboson"
        },

        {
            "WWTo2L2Nu",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WWTo2L2Nu.root"
            },
            11.09,
            9994000.0,
            "diboson"
        },

        {
            "WZTo2Q2L",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WZTo2Q2L.root"
            },
            6.419,
            28576996.0,
            "diboson"
        },

        {
            "WZTo3LNu",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WZTo3LNu.root"
            },
            5.257,
            9821283.0,
            "diboson"
        },

        {
            "ZZTo2L2Nu",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ZZTo2L2Nu.root"
            },
            0.9738,
            56886000.0,
            "diboson"
        },

        {
            "ZZTo2Q2L",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ZZTo2Q2L.root"
            },
            3.676,
            29357938.0,
            "diboson"
        },

        //================================================
        // TOP
        //================================================

        {
            "TTTo2L2Nu",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/TTTo2L2Nu.root"
            },
            88.29,
            145020000.0,
            "top"
        },

        {
            "ST_s-channel",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_s-channel_4f_leptonDecays.root"
            },
            3.549,
            19365999.0,
            "top"
        },

        {
            "ST_t-channel_top",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_t-channel_top_4f_InclusiveDecays.root"
            },
            119.7,
            178336000.0,
            "top"
        },

        {
            "ST_t-channel_antitop",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_t-channel_antitop_4f_InclusiveDecays.root"
            },
            67.93,
            95627000.0,
            "top"
        },

        {
            "ST_tW_top",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_tW_top_5f_NoFullyHadronicDecays.root"
            },
            34.91,
            11270430.0,
            "top"
        },

        {
            "ST_tW_antitop",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ST_tW_antitop_5f_NoFullyHadronicDecays.root"
            },
            34.91,
            10949620.0,
            "top"
        },

        {
            "tZq",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/tZq_ll_4f_ckm_NLO.root"
            },
            0.07561,
            11916000.0,
            "top"
        },

        {
            "TTWJets",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/TTWJetsToLNu.root"
            },
            0.2163,
            10450000.0,
            "top"
        },

        {
            "TTZ",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/TTZToLLNuNu.root"
            },
            0.2439,
            19608000.0,
            "top"
        },

        //================================================
        // TRIBOSON
        //================================================

        {
            "WWZ",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WWZ_4F.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WWZ_4F_ext.root"
            },
            0.1707,
            10209999.0,
            "triboson"
        },

        {
            "ZZZ",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ZZZ.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/ZZZ_ext.root"
            },
            0.01476,
            10139000.0,
            "triboson"
        },

        {
            "WZZ",
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WZZ.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WZZ_ext.root"
            },
            0.05709,
            10294000.0,
            "triboson"
        }
    };

    //================================================
    // CATEGORY HISTOGRAMS
    //================================================

    TH1F *h_top = new TH1F(
        "top",
        ";m_{T};Events",
        nbins,
        xmin,
        xmax
    );

    TH1F *h_diboson = new TH1F(
        "diboson",
        ";m_{T};Events",
        nbins,
        xmin,
        xmax
    );

    TH1F *h_triboson = new TH1F(
        "triboson",
        ";m_{T};Events",
        nbins,
        xmin,
        xmax
    );

    h_top->Sumw2();
    h_diboson->Sumw2();
    h_triboson->Sumw2();

    //================================================
    // LOOP OVER SAMPLES
    //================================================

    for (auto &s : samples) {

        if (s.category == "top") {

            add_sample_to_histogram(
                h_top,
                s,
                var,
                totalWeight,
                lumi,
                nbins,
                xmin,
                xmax
            );
        }

        else if (s.category == "diboson") {

            add_sample_to_histogram(
                h_diboson,
                s,
                var,
                totalWeight,
                lumi,
                nbins,
                xmin,
                xmax
            );
        }

        else if (s.category == "triboson") {

            add_sample_to_histogram(
                h_triboson,
                s,
                var,
                totalWeight,
                lumi,
                nbins,
                xmin,
                xmax
            );
        }
    }

    //================================================
    // WRITE GROUPED HISTOGRAMS
    //================================================

    fout->cd();

    h_top->Write("", TObject::kOverwrite);
    h_diboson->Write("", TObject::kOverwrite);
    h_triboson->Write("", TObject::kOverwrite);

    //------------------------------------------------
    // Systematics clones
    //------------------------------------------------

    ((TH1*)h_top->Clone("top_alphaUp"))
        ->Write("", TObject::kOverwrite);

    ((TH1*)h_top->Clone("top_alphaDown"))
        ->Write("", TObject::kOverwrite);

    ((TH1*)h_diboson->Clone("diboson_alphaUp"))
        ->Write("", TObject::kOverwrite);

    ((TH1*)h_diboson->Clone("diboson_alphaDown"))
        ->Write("", TObject::kOverwrite);

    ((TH1*)h_triboson->Clone("triboson_alphaUp"))
        ->Write("", TObject::kOverwrite);

    ((TH1*)h_triboson->Clone("triboson_alphaDown"))
        ->Write("", TObject::kOverwrite);
        
    yieldFile << left
        << setw(15) << "top"
        << h_top->Integral() << endl;

    yieldFile << left
        << setw(15) << "diboson"
        << h_diboson->Integral() << endl;

    yieldFile << left
        << setw(15) << "triboson"
        << h_triboson->Integral() << endl;

    ///////##########################################
    //====================================================
    // DATA
    //====================================================

    TFile *fdata = TFile::Open(
        "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/data/data_2018_noDuplicates.root"
    );

    TTree *tdata =
        (TTree*)fdata->Get("Events");

    TH1F *h_data = new TH1F(
        "data_obs",
        ";m_{T};Events",
        100,
        0,
        1000
    );

    tdata->Draw(
        "HZZ2l2nu_ZZmT >> data_obs",
        SRcut,
        "goff"
    );

    fout->cd();
    h_data->Write(
        "data_obs",
        TObject::kOverwrite
    );

    yieldFile << left
          << setw(15) << "data_obs"
          << h_data->Integral() << endl;
    //====================================================
    // SIGNAL
    //====================================================

    cout << "Processing signal..." << endl;

    TFile *fsig = TFile::Open(
        "/eos/user/a/avijay/Mergedrootfiles_hzz_new/v15/signal_v09/mass1000.root"
    );

    TTree *tsig =
        (TTree*)fsig->Get("Events");

    double sig_xsec = 39.71;
    double sig_ngen = 486250.0;

    double sig_scale =
        (sig_xsec * lumi) / sig_ngen;

    TH1F *h_sig = new TH1F(
        "signal",
        ";m_{T};Events",
        nbins,
        xmin,
        xmax
    );

    h_sig->Sumw2();

    tsig->Draw(
        var + " >> signal",
        totalWeight,
        "goff"
    );

    int Nsig = h_sig->GetNbinsX();

    h_sig->SetBinContent(
        Nsig,
        h_sig->GetBinContent(Nsig)
        + h_sig->GetBinContent(Nsig+1)
    );

    h_sig->Scale(sig_scale);

    fout->cd();
    h_sig->Write(
        "signal",
        TObject::kOverwrite
    );

    yieldFile << left
          << setw(15) << "signal"
          << h_sig->Integral() << endl;
    //====================================================
    // DY
    //====================================================

    cout << "Processing DY..." << endl;

    vector<tuple<TString,double,double>> dy_files = {

        {"drellyan_ptz_0to50.root",
         1485,
         196207761.0},

        {"drellyan_ptz_50to100.root",
         397.4,
         122967660.0},

        {"drellyan_ptz_100to250.root",
         97.2,
         79527324.0},

        {"drellyan_ptz_250to400.root",
         3.701,
         24195330.0},

        {"drellyan_ptz_400to650.root",
         0.5086,
         3936102.0},

        {"drellyan_ptz_650toInf.root",
         0.04728,
         3994997.0}
    };

    TH1F *h_dy = new TH1F(
        "dy",
        ";m_{T};Events",
        nbins,
        xmin,
        xmax
    );

    h_dy->Sumw2();

    for (auto &finfo : dy_files) {

        TString fname;
        double xsec, ngen;

        tie(fname, xsec, ngen) = finfo;

        TString fullpath =
            "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/"
            + fname;

        TFile *f = TFile::Open(fullpath);

        if (!f || f->IsZombie()) continue;

        TTree *tree =
            (TTree*)f->Get("Events");

        if (!tree) continue;

        double scale =
            (xsec * lumi) / ngen;

        TH1F *h_tmp = new TH1F(
            "h_tmp",
            "",
            nbins,
            xmin,
            xmax
        );

        h_tmp->Sumw2();

        tree->Draw(
            var + ">>h_tmp",
            totalWeight,
            "goff"
        );

        h_tmp->Scale(scale);

        h_dy->Add(h_tmp);

        delete h_tmp;

        f->Close();
    }

    int Ndy = h_dy->GetNbinsX();

    h_dy->SetBinContent(
        Ndy,
        h_dy->GetBinContent(Ndy)
        + h_dy->GetBinContent(Ndy+1)
    );

    TH1* h_dyUp =
        (TH1*)h_dy->Clone("dy_alphaUp");

    TH1* h_dyDown =
        (TH1*)h_dy->Clone("dy_alphaDown");

    fout->cd();

    h_dy->Write("", TObject::kOverwrite);
    h_dyUp->Write("", TObject::kOverwrite);
    h_dyDown->Write("", TObject::kOverwrite);


    yieldFile << left
          << setw(15) << "dy"
          << h_dy->Integral() << endl;
    //====================================================
    // WJETS
    //====================================================

    cout << "Processing WJets..." << endl;

    vector<WJetBin> wjets_bins = {

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-70To100.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-70To100_ext1.root"
            },
            1283.0,
            137713215.0
        },

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-100To200.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-100To200_ext1.root"
            },
            1244.0,
            120246316.0
        },

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-200To400.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-200To400_ext1.root"
            },
            337.8,
            114981255.0
        },

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-400To600.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-400To600_ext1.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-400To600_ext2.root"
            },
            44.93,
            19330678.0
        },

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-600To800.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-600To800_ext1.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-600To800_ext2.root"
            },
            11.19,
            29105780.0
        },

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-800To1200.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-800To1200_ext1.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-800To1200_ext2.root"
            },
            4.926,
            14179696.0
        },

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-1200To2500.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-1200To2500_ext1.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-1200To2500_ext2.root"
            },
            1.152,
            27123534.0
        },

        {
            {
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-2500ToInf_ext1.root",
                "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/WJetsToLNu_HT-2500ToInf_ext2.root"
            },
            0.02646,
            17298556.0
        }
    };

    TH1F *h_wjets = new TH1F(
        "wjets",
        ";m_{T};Events",
        nbins,
        xmin,
        xmax
    );

    h_wjets->Sumw2();

    for (auto &bin : wjets_bins) {

        TH1F *h_ht = new TH1F(
            "h_ht",
            "",
            nbins,
            xmin,
            xmax
        );

        h_ht->Sumw2();

        for (auto &fname : bin.files) {

            TFile *f = TFile::Open(fname);

            if (!f || f->IsZombie()) continue;

            TTree *tree =
                (TTree*)f->Get("Events");

            if (!tree) continue;

            TH1F *h_tmp = new TH1F(
                "h_tmp",
                "",
                nbins,
                xmin,
                xmax
            );

            h_tmp->Sumw2();

            tree->Draw(
                var + ">>h_tmp",
                totalWeight,
                "goff"
            );

            h_ht->Add(h_tmp);

            delete h_tmp;

            f->Close();
        }

        double scale =
            (bin.xsec * lumi)
            / bin.ngen_total;

        h_ht->Scale(scale);

        h_wjets->Add(h_ht);

        delete h_ht;
    }

    int Nw = h_wjets->GetNbinsX();

    h_wjets->SetBinContent(
        Nw,
        h_wjets->GetBinContent(Nw)
        + h_wjets->GetBinContent(Nw+1)
    );

    TH1* h_wjetsUp =
        (TH1*)h_wjets->Clone("wjets_alphaUp");

    TH1* h_wjetsDown =
        (TH1*)h_wjets->Clone("wjets_alphaDown");

    fout->cd();

    h_wjets->Write("", TObject::kOverwrite);
    h_wjetsUp->Write("", TObject::kOverwrite);
    h_wjetsDown->Write("", TObject::kOverwrite);

    yieldFile << left
          << setw(15) << "wjets"
          << h_wjets->Integral() << endl;

    //====================================================
    //////###########################################

    fout->Close();

    cout << "\nAll grouped backgrounds done!\n";

    return 0;
}