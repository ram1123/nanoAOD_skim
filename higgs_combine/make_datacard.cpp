#include <TFile.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

//----------------------------------------------------
double getYield(TFile *f, const string &hname) {

    TH1 *h = (TH1*)f->Get(hname.c_str());
    if (!h) {
        cout << "WARNING missing: " << hname << endl;
        return 0.0;
    }

    return h->Integral();
}

//----------------------------------------------------
string histName(string proc) {

    // signal stays signal_alpha
    if (proc == "signal_alpha") return "signal_alpha";

    // everything else uses _alpha
    return proc + "_alpha";
}

//----------------------------------------------------
int make_datacard() {

    TFile *f = TFile::Open("input_root_file_combine_v3.root");

    if (!f || f->IsZombie()) {
        cout << "Cannot open file" << endl;
        return 1;
    }

    //------------------------------------------------
    // EXACT PROCESS ORDER YOU WANT IN DATACARD
    //------------------------------------------------
    vector<string> proc = {

        "signal_alpha",
        "dy",
        "TTTo2L2Nu",
        "GluGluToContinToZZTo2mu2nu",
        "GluGluToContinToZZTo2e2nu",
        "GluGluToWWToENEN",
        "GluGluToWWToENMN",
        "GluGluToWWToMNMN",
        "ST_s-channel_4f_leptonDecays",
        "ST_t-channel_antitop_4f_InclusiveDecays",
        "ST_t-channel_top_4f_InclusiveDecays",
        "ST_tW_antitop_5f_NoFullyHadronicDecays",
        "ST_tW_top_5f_NoFullyHadronicDecays",
        "TTWJetsToLNu",
        "TTZToLLNuNu",
        "tZq_ll_4f_ckm_NLO",
        "WWTo2L2Nu",
        "WWZ_4F",
        "WZTo2Q2L",
        "WZTo3LNu",
        "WZZ",
        "ZZTo2L2Nu",
        "ZZTo2Q2L",
        "ZZZ"
    };

    //------------------------------------------------
    // GET YIELDS
    //------------------------------------------------
    map<string,double> yield;

    for (auto &p : proc) {

        string hname = histName(p);
        yield[p] = getYield(f, hname);
    }

    //------------------------------------------------
    // DATA
    //------------------------------------------------
    double data_obs = getYield(f, "data_obs");

    //------------------------------------------------
    // WRITE DATACARD
    //------------------------------------------------
    ofstream out("datacard.txt");

    int nproc = proc.size();

    out << "imax 1\n";
    out << "jmax " << (nproc - 1) << "\n";
    out << "kmax *\n";
    out << "----------------\n";

    out << "shapes * * input_root_file_combine_v3.root $PROCESS $PROCESS_$SYSTEMATIC\n";
    out << "----------------\n";

    out << "bin bin1\n";
    out << "observation " << data_obs << "\n";
    out << "------------------------------\n";

    //------------------------------------------------
    // BIN LINE
    //------------------------------------------------
    out << "bin ";
    for (int i = 0; i < nproc; i++) out << "bin1 ";
    out << "\n";

    //------------------------------------------------
    // PROCESS NAMES (EXACT ORDER YOU WANT)
    //------------------------------------------------
    out << "process ";
    for (auto &p : proc) out << p << " ";
    out << "\n";

    //------------------------------------------------
    // PROCESS IDS
    //------------------------------------------------
    out << "process ";
    for (int i = 0; i < nproc; i++) out << i << " ";
    out << "\n";

    //------------------------------------------------
    // RATE
    //------------------------------------------------
    out << "rate ";
    for (auto &p : proc) {
        out << yield[p] << " ";
    }
    out << "\n";

    //------------------------------------------------
    // SYSTEMATIC LINE
    //------------------------------------------------
    out << "--------------------------------\n";
    out << "alpha shape ";

    for (auto &p : proc) {

        if (p == "signal_alpha") {
            out << "- ";
        } else {
            out << "1 ";
        }
    }

    out << "\n";

    out.close();

    cout << "\n Datacard created successfully: datacard.txt\n";

    return 0;
}