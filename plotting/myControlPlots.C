// #include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
// #include <string>
#include <vector>
// #include <assert.h>
// #include "TROOT.h"
// #include "TLatex.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
// #include "TBox.h"
// #include "TH1F.h"
#include "THStack.h"
// #include "TString.h"
#include "TLegend.h"
// #include "TLine.h"
// #include "TPad.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TSystem.h"
// #include "TF1.h"
// #include "TFitResult.h"
// #include "TFitResultPtr.h"
#include "TError.h"

#include "tdrstyle.C"
#include "utils.C" // Tokenize
#include "CMS_lumi.C"
#include <cmath>

// #include <Python.h>

typedef struct
{
    int index;
    TString samplename;
    TString treefilename;
    double xsecpblumi;
    double otherscale;
    int nMCevents;
    int MCnegEvent;
    int colorcode;
    int stackit;
} SampleInfo_t;

#include "plotvar_t.h"

double intLUMIinvpb;

// Per-event MC weight multiplied into the cut string. Default keeps the old behaviour;
// a RunMacro can override it, e.g. gMCWeight = "overallEventWeight" (genWeight*puWeight*L1PreFiringWeight_Nom).
TString gMCWeight = "puWeight";
// Optional per-event data weight (e.g. "gamma_prescaleWeight" for prescaled photon triggers); empty = unweighted.
TString gDataWeight = "";

//=====================================================================================
// SYNOPSIS:
//   1. Prepare "InData" and "OutDir" directories; e.g., "ln -s . OutDir" to go to current dir
//   2. Prepare "cuttable.txt" of cut names and cut strings
//   3. root [0] .L myControlPlots.C+
//      root [1] myControlPlots("cuttable.txt","sampleFileName.txt")
//
// ====================================================================================
// Self Function
// ====================================================================================

void cmspre(double intlumifbinv)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);

    latex.SetTextAlign(31); // align right
    latex.DrawLatex(0.95, 0.92, "#sqrt{s} = 13 TeV");
    latex.SetTextAlign(31); // align right
    latex.DrawLatex(0.73, 0.93, Form("#int #font[12]{L} dt = %.1f fb^{-1}", (float)intlumifbinv));

    latex.SetTextAlign(11); // align left
    latex.SetTextSize(0.06);
    latex.DrawLatex(0.145, 0.92, "CMS #bf{#it{preliminary}}");
}

//======================================================================

void loadCutString(const char *filename, TString &cutstring)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        std::cout << "Error, file " << TString(filename) << " not found." << std::endl;;;
        exit(-1);
    }

    char line[512];

    for (int i = 0; !feof(fp) && fgets(line, 512, fp); i++)
    {
        if (!strlen(line) || line[0] == '#')
            continue; // comments are welcome

        if (cutstring.Length())
            cutstring += " && ";

        std::string strline(line);
        strline.pop_back(); // shed the \n
        std::vector<std::string> fields;

        // expect columns with fields cutname, cutvalue, possible embedded spaces both
        // within and between, so " " or "\t" cannot be used as delimiters. Require quotes
        // instead.
        //
        Tokenize(strline, fields, "\"");

        // for (size_t j=0; j<fields.size(); j++)
        // std::cout << j << ": \"" << fields[j] << "\"" << std::endl;;

        assert(fields.size() == 3);
        cutstring += TString(fields.at(2));
    }
} // loadCutString

//======================================================================

class Sample
{
public:
    Sample(const SampleInfo_t &sinfo)
    {
        info_ = sinfo;
        tree_ = 0;
        // std::cout << "sample = " << name_ << std::endl;;
        TFile *f = new TFile(sinfo.treefilename, "READ");
        if (!f)
        {
            std::cerr << "Couldn't find file " << sinfo.treefilename << std::endl;;
            return;
        }
        // TFile *f =  TFile::Open("root://cmsxrootd.fnal.gov/"+sinfo.treefilename, "READ"); if (!f) { std::cerr << "Couldn't find file " << sinfo.treefilename << std::endl;; return; }
        // tree_ =  (TTree *)f->Get("otree"); if (!tree_) { std::cerr << "Couldn't find tree otree in file " << sinfo.treefilename << std::endl;; return; }
        tree_ = (TTree *)f->Get("Events");
        if (!tree_)
        {
            std::cerr << "Couldn't find tree Events in file " << sinfo.treefilename << std::endl;;
            return;
        }
    }
    ~Sample()
    {
        if (tree_)
            delete tree_;
    }
    TTree *Tree() const { return tree_; }
    TString name() const { return info_.samplename; }
    TString filename() const { return info_.treefilename; }
    bool stackit() const { return (info_.stackit != 0); }
    int colorcode() const { return info_.colorcode; }
    double otherscale() const { return info_.otherscale; }
    TH1D *Draw(const plotVar_t &pv, const TCut &cut, const TCut &cutSQ)
    {
        if (!tree_)
        {
            std::cerr << "No tree to draw from." << std::endl;;
            return 0;
        }
        double tmp = 0.;
        std::cout << "\tDrawing " << pv.plotvar << " for sample = " << info_.samplename << " ... ";
        TString hname = TString("th1") + pv.outfile + Form("%d", info_.index);
        TH1D *histo = new TH1D(hname, hname, pv.NBINS, pv.MINRange, pv.MAXRange);
        assert(histo);
        histo->Sumw2();
        std::cout << tree_->Draw(pv.plotvar + TString(">>") + hname, cut, "goff") << " entries, ";
        histo->SetBinContent(pv.NBINS, histo->GetBinContent(pv.NBINS) + histo->GetBinContent(pv.NBINS + 1));
        // std::cout << histo->Integral() << " " << tmp << " " << "weighted entries, (";
        std::cout << histo->IntegralAndError(0, histo->GetNbinsX() + 1, tmp) << " " << tmp << " "
             << "weighted entries";

        if (info_.nMCevents)
        {
            // std::cout<<"\n===> Evetns = "<<info_.xsecpblumi<<"\t"<<info_.nMCevents<<"\t"<<info_.MCnegEvent<<"\t"<<info_.colorcode<<std::endl;;
            histo->Scale((info_.xsecpblumi * info_.otherscale) / (info_.nMCevents - 2 * info_.MCnegEvent));
            std::cout << ", " << histo->IntegralAndError(1, histo->GetNbinsX(), tmp) << " " << tmp << " " << (histo->Integral(1, histo->GetNbinsX() + 1) * info_.xsecpblumi * info_.otherscale) / (info_.nMCevents - 2 * info_.MCnegEvent) << " "
                 << " scaled events in window";
        }
        std::cout << std::endl;;

        return histo;
    }

private:
    SampleInfo_t info_;
    TTree *tree_;
};

//======================================================================

void loadSamples(const char *filename, std::vector<Sample *> &samples)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        std::cout << "Error, file " << TString(filename) << " not found." << std::endl;;
        exit(-1);
    }

    char line[512];

    intLUMIinvpb = -1; // obvious error condition

    for (int i = 0; !feof(fp) && fgets(line, 512, fp); i++)
    {
        if (!strlen(line) || line[0] == '#')
            continue; // comments are welcome

        std::string strline(line);
        strline.pop_back(); // shed the \n
        std::vector<std::string> fields;

        // expect columns with fields cutname, cutvalue, possible embedded spaces both
        // within and between, so " " or "\t" cannot be used as delimiters. Require quotes
        // instead.
        //
        Tokenize(strline, fields, " \t");

        // for (size_t j=0; j<fields.size(); j++)
        // std::cout << j << ": \"" << fields[j] << "\"" << std::endl;;

        assert(fields.size() == 8);

        SampleInfo_t s;
        s.index = i;
        s.samplename = fields[0];
        s.treefilename = fields[1];
        s.xsecpblumi = str2dbl(fields[2]);
        s.otherscale = str2dbl(fields[3]);
        s.nMCevents = str2int(fields[4]);
        s.MCnegEvent = str2int(fields[5]);
        s.colorcode = str2int(fields[6]);
        s.stackit = str2int(fields[7]);

        std::cout << "Loading sample " << s.samplename << " -> " << s.treefilename << std::endl;;

        if (!samples.size())
        {
            if (s.samplename.EqualTo("data"))
            {
                intLUMIinvpb = s.xsecpblumi;
                s.xsecpblumi = 1;
                std::cout << "intLUMI = " << intLUMIinvpb << " pb^-1" << std::endl;;
            }
            else
            {
                std::cerr << "First sample in the table must be 'data'" << std::endl;;
                exit(-1);
            }
        }
        else
        {
            s.otherscale *= intLUMIinvpb;
        }

        samples.push_back(new Sample(s));
    }
} // loadSamples

//======================================================================

void myControlPlots(const char *cuttablefilename,
                    const char *samplefilename,
                    const plotVar_t plotvars[],
                    const std::string OutRootFile = "testrk.root",
                    const std::string OutDir = "OutDir",
                    const int ScaleSignal = 0,
                    const std::string RecreateAppend = "RECREATE",
                    const int isData = 1)
{
    gErrorIgnoreLevel = kError;

    if (int(ScaleSignal) != 1 && int(ScaleSignal) != 0)
    {
        std::cout << "=========================================\n\n"
                  << std::endl;;
        std::cerr << "Error: Scale Signal values should be  1 or 0.\n"
                  << std::endl;;
        exit(-1);
    }

    std::ofstream Logfile;

    TH1D::SetDefaultSumw2(1);

    // TFile *outDC = new TFile("htt_mt_inputs.root","RECREATE");

    TString unwtcutstring;

    loadCutString(cuttablefilename, unwtcutstring);

    TString blinddatacutstring = unwtcutstring + "&&(1)";
    //	For blinding data uncomment this
    // TString blinddatacutstring = unwtcutstring + "&& (PuppiAK8_jet_mass_so<65 || PuppiAK8_jet_mass_so>105)&&(1)";

    //  const char* the_cut = "1";
    //  double BINWIDTH = ((MAXRange-MINRange)/NBINS);

    // Get the input trees:

    std::vector<Sample *> samples;

    loadSamples(samplefilename, samples);

    // Data
    Sample *sdata = samples[0];

    if (sdata->Tree())
        std::cout << "ndata =" << sdata->Tree()->GetEntries() << std::endl;;
    TFile *f;

    if (strcmp(RecreateAppend.c_str(), "RECREATE") == 0)
    {
        f = new TFile(OutRootFile.c_str(), "RECREATE");
        std::cout << "RECREATE file " << OutRootFile.c_str() << std::endl;;
    }
    else if (strcmp(RecreateAppend.c_str(), "UPDATE") == 0)
    {
        f = new TFile(OutRootFile.c_str(), "UPDATE");
        std::cout << "UPDATE file " << OutRootFile.c_str() << std::endl;;
    }
    else
    {
        std::cout << "=========================================\n\n"
                  << std::endl;;
        std::cerr << "Error: Root file option can either be RECREATE or UPDATE...\n"
                  << std::endl;;
        exit(-1);
    }
    //============================================================
    //  VARIABLE LOOP
    //============================================================

    for (int ivar = 0;; ivar++)
    {
        plotVar_t pv = plotvars[ivar];
        //  Create OutDir
        gSystem->mkdir(OutDir.c_str(), kTRUE);
        // copy file index.php from current directory  to OutDir
        gSystem->CopyFile("index.php", TString(OutDir) + TString("/index.php"));

        TString outfile = TString(OutDir) + TString("/") + TString(gSystem->BaseName(cuttablefilename)).ReplaceAll(".txt", "") + TString("_") + pv.outfile;
        Logfile.open(outfile + ".log");
        TString temp = TString("Yield of Variable : ") + pv.outfile;
        const std::string spaces(temp.Length(), ' ');
        const std::string second = "* " + spaces + " *";
        const std::string first(second.size(), '*');
        Logfile << first << std::endl;;
        Logfile << second << std::endl;;
        Logfile << "* " << temp << " *" << std::endl;;
        Logfile << second << std::endl;;
        if (int(ScaleSignal) == 1)
            Logfile << "* \t signal is scaled by 50..." << std::endl;;
        Logfile << second << std::endl;;
        Logfile << first << std::endl;;

        if (!pv.plotvar.Length())
            break;

        std::cout << pv.plotvar << "\t" << pv.MINRange << "\t" << pv.MAXRange << "\t" << pv.NBINS << "\tTHE CUT " << std::endl;;

        if (!pv.plotvar.Contains("_smeared"))
            if (sdata->Tree() &&
                sdata->Tree()->Draw(pv.plotvar, "", "goff", 1) == -1)
            { // check if the variable exists in the tree
                std::cout << "\t...can't be plotted!" << std::endl;;
                continue;
            }

        //
        //	ONE LEPTON CUTS
        //
        // TCut the_cut(TString("L1_PrefweightUp*btag0Wgt*genWeight*trig_eff_Weight*id_eff_Weight*pu_Weight*(")+unwtcutstring+TString(")"));
        // TCut the_cut(TString("L1_PrefweightDown*btag0Wgt*genWeight*trig_eff_Weight*id_eff_Weight*pu_Weight*(")+unwtcutstring+TString(")"));
        // TCut the_cut(TString("*(") + unwtcutstring + TString(")"));
        //TCut the_cut(TString("1*(") + unwtcutstring + TString(")")); // PV weight not added

        TCut the_cut(gMCWeight + TString("*(") + unwtcutstring + TString(")"));  //####AV### (was hard-coded puWeight)
        // TCut the_cut(TString("trig_eff_Weight*btag0Wgt*genWeight*id_eff_Weight*pu_Weight*(")+unwtcutstring+TString(")"));
        // TCut the_cut(TString("btag1Wgt*genWeight*trig_eff_Weight*id_eff_Weight*pu_Weight*(")+unwtcutstring+TString(")"));	// For Top control region

        TCut nullcut("");

        std::cout << "the_cut = " << the_cut << std::endl;;

        const double BINWIDTH = ((pv.MAXRange - pv.MINRange) / pv.NBINS);

        map<TString, TH1D *> m_histos;
        map<TString, bool> m_stacked;

        double totevents = 0.;

        //============================================================
        // DRAW THE VARIABLE FOR ALL SAMPLES, CREATE HISTOS
        //============================================================

        for (size_t isamp = 0; isamp < samples.size(); isamp++)
        {
            Sample *s = samples[isamp];

            m_stacked[s->name()] = false;

            TH1D *h;

            if (s->name().EqualTo("data"))
            {
                plotVar_t pvnosmear = pv;
                // pvnosmear.plotvar.ReplaceAll("_smeared",""); // no smearing for data!
                TString datacut = gDataWeight.Length() ? gDataWeight + "*(" + blinddatacutstring + ")" : blinddatacutstring;
                h = s->Draw(pvnosmear, TCut(datacut), nullcut); // effwt*puwt==1 for data! -- NO IT DOESN'T NECESSARILY!
            }
            // else if (s->name().EqualTo("aQGCX100")){
            else if (s->name().EqualTo("aQGC"))
            {
                if (int(ScaleSignal) == 1)
                {
                    h = s->Draw(pv, the_cut * "50" * "(LHEWeight[1121]/LHEWeight[0])", the_cut * "50" * "(LHEWeight[1121]/LHEWeight[0])");
                    std::cout << "====> Scale Signal = " << ScaleSignal << std::endl;;
                }
                else
                    h = s->Draw(pv, the_cut * "(LHEWeight[1121]/LHEWeight[0])", the_cut * "(LHEWeight[1121]/LHEWeight[0])");
                if (s->stackit())
                {
                    totevents += h->Integral(1, h->GetNbinsX() + 1);
                }
            }
            // else if (s->name().EqualTo("WV_EWKX100")){
            else if (s->name().EqualTo("WV_EWK"))
            {
                if (int(ScaleSignal) == 1)
                    h = s->Draw(pv, the_cut * "50", the_cut * "50");
                else
                    h = s->Draw(pv, the_cut, the_cut);
                if (s->stackit())
                {
                    totevents += h->Integral(1, h->GetNbinsX() + 1);
                }
            }
            //else if (s->name().EqualTo("Z+jets"))
            else if (s->name().EqualTo("DYJets"))
            {
                // add additional cut to the cut string: Weight_nJets_FromDataMC
                //h = s->Draw(pv, the_cut* "Weight_nJets_FromDataMC", the_cut* "Weight_nJets_FromDataMC");  //##AV##
		//h = s->Draw(pv, the_cut* "Weight_nJets_FromDataMC_ge1jet", the_cut* "Weight_nJets_FromDataMC_ge1jet");  //##AV##
		//std::cout << "\t\t==> Updated cut: " << the_cut* "Weight_nJets_FromDataMC" << std::endl;;
               h = s->Draw(pv,  the_cut,  the_cut);    //uncommented ##AV##
                if (s->stackit())
                {
                    totevents += h->Integral(1, h->GetNbinsX() + 1);
                }
            }
            else
            {
                h = s->Draw(pv, the_cut, the_cut);
                // scaling for Non resonant backgrounds.
                const double NRB_SF = 1; //0.0786

                TString fname = s->filename();

                bool isNRB = false;

                // Top samples
                if (s->name().EqualTo("Top"))
                isNRB = false;

                // WJets samples
                if (s->name().EqualTo("WJets"))
                isNRB = false;

                // WW samples
                if (fname.Contains("WWTo2L2Nu") ||
                fname.Contains("GluGluToWW"))
                isNRB = true;

                // WWZ only
                if (fname.Contains("WWZ"))
                isNRB = false;

                //// WWW only
                if (fname.Contains("WWW"))
                isNRB = true;

                if (isNRB)
                {
                    h->Scale(NRB_SF);

                std::cout << " --> Applied NRB SF = " << NRB_SF << " to " << fname << std::endl;
    }

                if (s->stackit())
                {
                    totevents += h->Integral(1, h->GetNbinsX() + 1);
                }
            }

            if (!h)
                continue;

            map<TString, TH1D *>::iterator mit = m_histos.find(s->name());
            if (mit == m_histos.end())
            {
                h->SetName(s->name());
                h->SetTitle(s->name());
                if (s->stackit())
                {
                    if (s->name().EqualTo("Triboson"))
                    {
                        h->SetLineColor(TColor::GetColor(222, 90, 106));
                        h->SetFillColor(TColor::GetColor(222, 90, 106));
                        h->SetLineWidth(0);
                    }
                    else if (s->name().EqualTo("W+jets"))
                    {
                        h->SetLineColor(TColor::GetColor(222, 90, 106));
                        h->SetFillColor(TColor::GetColor(222, 90, 106));
                        h->SetLineWidth(0);
                    }
                    else if (s->name().EqualTo("top"))
                    {
                        h->SetLineColor(TColor::GetColor(155, 152, 204));
                        h->SetFillColor(TColor::GetColor(155, 152, 204));
                        h->SetLineWidth(0);
                    }
                    else if (s->name().EqualTo("TTbar"))
                    {
                        h->SetLineColor(TColor::GetColor(155, 152, 204));
                        h->SetFillColor(TColor::GetColor(155, 152, 204));
                        h->SetLineWidth(0);
                    }
                    else if (s->name().EqualTo("SingleTop"))
                    {
                        h->SetLineColor(TColor::GetColor(0, 84, 159));
                        h->SetFillColor(TColor::GetColor(0, 84, 159));
                        h->SetLineWidth(0);
                    }
                    else if (s->name().EqualTo("Z+jets"))
                    {
                        h->SetLineColor(TColor::GetColor(248, 206, 104));
                        h->SetFillColor(TColor::GetColor(248, 206, 104));
                        h->SetLineWidth(0);
                    }
                    else if (s->name().EqualTo("Diboson"))
                    {
                        h->SetLineColor(TColor::GetColor(250, 202, 255));
                        h->SetFillColor(TColor::GetColor(250, 202, 255));
                        h->SetLineWidth(0);
                    }
                    else
                    {
                        h->SetLineColor(s->colorcode());
                        h->SetFillColor(s->colorcode());
                        h->SetLineWidth(0);
                    }
                }
                m_histos[s->name()] = h;
            }
            else
            {
                mit->second->Add(h);
            }
        }

        //============================================================
        // COUNT EVENTS, RENORM TO DATA, CONSTRUCT THE TSTACK & LEGEND
        //============================================================

        TH1D *th1data = m_histos["data"];
        // TH1D *aqgc = m_histos["aQGC"];
        double maxval = -9e99;

        double ndata = 1., renorm = 1.;
        if (th1data)
        {
            ndata = th1data->Integral(1, th1data->GetNbinsX() + 1);
            renorm = ndata / totevents;
            int maxbin = th1data->GetMaximumBin();
            maxval = std::max(maxval, th1data->GetBinContent(maxbin));

            std::cout << "Total MC	= " << totevents << std::endl;;
            std::cout << "data	= " << ndata << std::endl;;
            std::cout << "data/MC	= " << renorm << std::endl;;
            Logfile << "Total MC	= " << totevents << std::endl;;
            Logfile << "data	= " << ndata << std::endl;;
            Logfile << "data/MC	= " << renorm << std::endl;;
        }

        // Setup the stack and total
        THStack *hs = new THStack("hs", "MC contribution");
        TH1D *th1tot = new TH1D("th1tot", "th1tot", pv.NBINS, pv.MINRange, pv.MAXRange);

        // Set up the legend
        float legX0 = 0.65, legX1 = 0.99, legY0 = 0.54, legY1 = 0.88;
        TLegend *Leg = new TLegend(legX0, legY0, legX1, legY1);
        Leg->SetFillColor(0);
        Leg->SetFillStyle(0);
        Leg->SetBorderSize(0);

        if (th1data)
        {
            if (TString(cuttablefilename).Contains("Mu"))
                Leg->AddEntry(th1data, "Observed", "LP");
            else
                Leg->AddEntry(th1data, "Observed", "LP");
        }

        std::vector<double> binErrSQ(pv.NBINS, 0.);

        std::vector<std::pair<TString, TH1D *>> v_legentries;

        TString oldsamplename;

        for (size_t isamp = 1; isamp < samples.size(); isamp++)
        {
            Sample *s = samples[isamp];
            if (s->name() == oldsamplename)
                continue;

            map<TString, TH1D *>::iterator mit = m_histos.find(s->name());
            TH1D *h = mit->second;
            h->SetLineWidth(3.);
            h->SetTitleSize(0.055, "Y");
            h->SetTitleOffset(1.600, "Y");
            h->SetLabelOffset(0.014, "Y");
            h->SetLabelSize(0.050, "Y");
            h->SetLabelFont(42, "Y");
            h->SetTitleSize(0.055, "X");
            h->SetTitleOffset(1.300, "X");
            h->SetLabelOffset(0.014, "X");
            h->SetLabelSize(0.050, "X");
            h->SetLabelFont(42, "X");
            h->SetMarkerSize(0.6);
            h->GetYaxis()->SetTitleFont(42);
            h->GetXaxis()->SetTitleFont(42);

            std::cout << s->name() << "	= " << Form("%7g", h->Integral(1, h->GetNbinsX() + 1)) << std::endl;;
            Logfile << s->name() << "	= " << Form("%7g", h->Integral(1, h->GetNbinsX() + 1)) << std::endl;;

            if (s->stackit())
            {
                // Compose the stack
                hs->Add(h);
                th1tot->Add(h);
                m_stacked[s->name()] = true;
                for (int ibin = 1; ibin <= pv.NBINS; ibin++)
                    binErrSQ[ibin - 1] += h->GetBinError(ibin) * h->GetBinError(ibin);
            }
            v_legentries.push_back(*mit);

            int maxbin = h->GetMaximumBin();
            maxval = std::max(maxval, h->GetBinContent(maxbin));

            oldsamplename = s->name();
            h->Write();
        }

        std::cout << "maxval " << maxval << std::endl;;

        // Reverse the order for the legend
        for (std::vector<std::pair<TString, TH1D *>>::reverse_iterator
                 rit = v_legentries.rbegin();
             rit != v_legentries.rend();
             rit++)
        {
            if (rit->first == "aQGC")
            {
                if (int(ScaleSignal) == 1)
                    Leg->AddEntry(rit->second, rit->first + TString("X50"), "L"); // "F");
                else
                    Leg->AddEntry(rit->second, rit->first + TString("(FT2=-0.5)"), "L"); // "F");
            }
            else if (rit->first == "WV_EWK")
            {
                if (int(ScaleSignal) == 1)
                    Leg->AddEntry(rit->second, rit->first + TString("X50"), "L"); // "F");
                else
                    Leg->AddEntry(rit->second, rit->first, "L"); // "F");
            }
            else
                Leg->AddEntry(rit->second, rit->first, "F"); // "F");
        }
        TH1D *th1totClone = (TH1D *)th1tot->Clone("th1totClone");
        th1totClone->SetMarkerStyle(0);
        th1totClone->SetFillStyle(3003);
        th1totClone->SetFillColor(11);
        th1totClone->SetLineColor(0);

        for (int ibin = 1; ibin <= th1totClone->GetNbinsX(); ++ibin)
        {
            th1totClone->SetBinError(ibin, std::sqrt(binErrSQ[ibin - 1]));
        }

        //============================================================
        // SETUP THE CANVAS
        //============================================================
        setTDRStyle();

        TCanvas *c1 = new TCanvas(pv.plotvar, pv.plotvar, 10, 10, 800, 600);
        TPad *d1;
        TPad *d2;

        c1->SetFillColor(0);
        c1->SetBorderMode(0);
        c1->SetBorderSize(10);
        // Set margins to reasonable defaults
        c1->SetLeftMargin(0.18);
        c1->SetRightMargin(0.05);
        c1->SetTopMargin(0.08);
        c1->SetBottomMargin(0.15);
        // Setup a frame which makes sense
        c1->SetFrameFillStyle(0);
        c1->SetFrameLineStyle(0);
        c1->SetFrameBorderMode(0);
        c1->SetFrameBorderSize(10);
        c1->SetFrameFillStyle(0);
        c1->SetFrameLineStyle(0);
        c1->SetFrameBorderMode(0);
        c1->SetFrameBorderSize(10);

        if (th1data)
        {
            c1->Divide(1, 2);
            d1 = (TPad *)c1->GetPad(1);
            d1->SetPad(0.01, 0.30, 0.95, 0.99);
            // c1->cd();
            d2 = (TPad *)c1->GetPad(2);
            d2->SetPad(0.01, 0.02, 0.95, 0.30);
            d1->cd();
            gPad->SetBottomMargin(0.0);
        }

        gPad->SetTopMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetLeftMargin(0.14);

        gPad->SetLogy(0);

        if (th1tot->GetEntries())
            Leg->AddEntry(th1tot, "Stat. Uncertainty", "f");

        Leg->SetFillColor(0);

        TH1D *th1totempty = new TH1D("th1totempty", "th1totempty", pv.ANBINS, pv.AMINRange, pv.AMAXRange);
        if (th1data)
        {
            th1data->SetMarkerStyle(20);
            th1data->SetLineColor(kBlack);
            th1data->SetMarkerSize(1.3);
            th1data->SetLineWidth(3);
        }

        th1tot->SetFillStyle(3013);
        th1tot->SetFillColor(13);
        th1tot->SetLineColor(kBlack);
        th1tot->SetLineWidth(1);
        th1tot->SetMarkerStyle(0);

        char tmpc[100];
        sprintf(tmpc, "Events");
        if (pv.slog == 1)
            sprintf(tmpc, "Events/ %.1f", BINWIDTH);
        if (pv.slog == 2)
            sprintf(tmpc, "Events/ %.3f", BINWIDTH);
        if (pv.slog == 3)
            sprintf(tmpc, "Events/ ( %.0f GeV)", BINWIDTH);
        if (pv.slog == 6)
            sprintf(tmpc, "Events/ ( %.1f rad)", BINWIDTH);
        th1totempty->SetYTitle(tmpc);
        th1totempty->GetYaxis()->SetTitleOffset(1.2);
        th1totempty->GetYaxis()->SetLabelOffset(0.01);
        th1totempty->GetYaxis()->SetLabelSize(0.04);

        if (th1tot->GetEntries())
        {
            th1tot->SetMinimum(0.1);
            th1totempty->SetMinimum(0.1);
        }

        //    th1totempty->SetMaximum(2.5*maxval);
        th1totempty->SetMaximum(1.9 * maxval);
        if (pv.slog == 9)
            th1totempty->SetMaximum(1.9 * maxval);
        // th1totempty->SetMaximum(350.);

        // Draw it all

        th1totempty->Draw();

        if (th1data)
        {
            th1data->SetMinimum(0.1);
            th1data->SetBinErrorOption(TH1D::kPoisson);
            th1data->Draw("e0same");
        }

        // th1tot->Draw("e2same");

        hs->SetMaximum(400.0);
        hs->Draw("samehist");
        if (pv.drawleg == 1)
            Leg->Draw();

        th1tot->Draw("e2same");

        if (th1data)
            th1data->Draw("e1same");
        // if(aqgc)
        //   aqgc->Draw("same");

        // draw unstacked histos separately.
        oldsamplename.Clear();
        TH1D *h1 = new TH1D("h1", "", 100, 0, 100);
        TH1D *h2 = new TH1D("h2", "", 100, 0, 100);
        for (size_t isamp = 0; isamp < samples.size(); isamp++)
        {
            Sample *s = samples[isamp];
            assert(s);
            if (s->name() == oldsamplename)
                continue;
            if (!s->stackit())
            {
                map<TString, TH1D *>::iterator mit = m_histos.find(s->name());
                if (mit != m_histos.end())
                {
                    TH1D *h = mit->second;
                    if (h && s->name() == "WV_EWK")
                    // if (h && s->name()=="WV_EWKX100")
                    {
                        h->SetFillStyle(0.);
                        // aqgc->SetLineStyle(11);
                        h->SetLineWidth(3.);
                        h->SetLineColor(kBlue + 3);
                        if (int(ScaleSignal) == 1)
                        {
                            std::cout << "Significance (SM EWK) = " << (h->Integral(1, h->GetNbinsX() + 1) / 50) / std::sqrt((h->Integral(1, h->GetNbinsX() + 1) / 50) + totevents) << std::endl;;
                            Logfile << "Significance (SM EWK) = " << (h->Integral(1, h->GetNbinsX() + 1) / 50) / std::sqrt((h->Integral(1, h->GetNbinsX() + 1) / 50) + totevents) << std::endl;;
                        }
                        else
                        {
                            std::cout << "Significance (SM EWK) = " << (h->Integral(1, h->GetNbinsX() + 1)) / std::sqrt((h->Integral(1, h->GetNbinsX() + 1)) + totevents) << std::endl;;
                            Logfile << "Significance (SM EWK) = " << (h->Integral(1, h->GetNbinsX() + 1)) / std::sqrt((h->Integral(1, h->GetNbinsX() + 1)) + totevents) << std::endl;;
                        }
                        h->Draw("histsame");
                        h1 = (TH1D *)h->Clone();
                    }
                    if (h && s->name() == "aQGC")
                    // if (h && s->name()=="aQGCX100")
                    {
                        h->SetFillStyle(0.);
                        // aqgc->SetLineStyle(11);
                        h->SetLineWidth(3.);
                        h->SetLineColor(kRed + 3);
                        // std::cout << "Significance (aQGC)   = " << (h->Integral(1,h->GetNbinsX()+1))/std::sqrt((h->Integral(1,h->GetNbinsX()+1))+totevents) << std::endl;;
                        std::cout << "Significance (aQGC)   = " << (h->Integral(1, h->GetNbinsX() + 1) / 50) / std::sqrt((h->Integral(1, h->GetNbinsX() + 1) / 50) + totevents) << std::endl;;
                        // Logfile << "Significance (aQGC)   = " << (h->Integral(1,h->GetNbinsX()+1))/std::sqrt((h->Integral(1,h->GetNbinsX()+1))+totevents) << std::endl;;
                        Logfile << "Significance (aQGC)   = " << (h->Integral(1, h->GetNbinsX() + 1) / 50) / std::sqrt((h->Integral(1, h->GetNbinsX() + 1) / 50) + totevents) << std::endl;;
                        h->Draw("histsame");
                        h2 = (TH1D *)h->Clone();
                        // h->Draw("e1same");
                    }
                    if (isData)
                        h->Write();
                }
            }
            oldsamplename = s->name();
        }

        // cmspre(intLUMIinvpb/1000.0);
        //CMS_lumi(d1, 4, 10);
        // th1data->Draw("Axissame");
        gPad->RedrawAxis();
        TH1D *hhratio;
        TH1D *th1emptyclone;
        TBox *errbox;

        if (th1data)
        {
            d2->cd();
            // TH1D    * hhratio    = (TH1D*) th1data->Clone("hhratio")  ;
            hhratio = (TH1D *)th1data->Clone("hhratio");
            hhratio->Sumw2();
            hhratio->SetStats(0);

            gPad->SetLeftMargin(0.14);
            gPad->SetTopMargin(0);
            gPad->SetRightMargin(0.05);
            gPad->SetFrameBorderSize(0);
            gPad->SetBottomMargin(0.3);
            gPad->SetTickx();

            hhratio->SetMarkerSize(1.25);
            //  hhratio->GetYaxis()->SetRangeUser(0.48,1.52);
            hhratio->GetYaxis()->SetRangeUser(0.5, 1.5);
            hhratio->GetXaxis()->SetTitle(pv.xlabel);
            hhratio->GetXaxis()->SetTitleOffset(0.9);
            hhratio->GetXaxis()->SetTitleSize(0.15);
            hhratio->GetXaxis()->SetLabelSize(0.15);
            hhratio->GetYaxis()->SetTitleSize(0.1);
            hhratio->GetYaxis()->SetTitleOffset(0.5);
            hhratio->GetYaxis()->CenterTitle(true);
            hhratio->GetYaxis()->SetLabelSize(0.1);
            std::cout << hhratio->GetNbinsX() << std::endl;;
            std::cout << th1tot->GetNbinsX() << std::endl;;
            hhratio->Divide(th1tot);
            double binError(0.0), mcbinentry(0.0), mcerror(0.0);
            for (int i = 0; i < hhratio->GetNbinsX(); ++i)
            {
                binError = hhratio->GetBinError(i);
                mcerror = th1tot->GetBinError(i);
                mcbinentry = th1tot->GetBinContent(i);
                if (mcbinentry > 0.)
                    mcerror /= mcbinentry;
                else
                    mcerror = 0.0;
                binError = std::sqrt(binError * binError + mcerror * mcerror);
                hhratio->SetBinError(i, binError);
            }
            th1emptyclone = new TH1D("th1emptyclone", "th1emptyclone", pv.ANBINS, pv.AMINRange, pv.AMAXRange);
            //th1emptyclone->GetYaxis()->SetRangeUser(0.1, 2.8000);
	    th1emptyclone->GetYaxis()->SetRangeUser(0.5, 1.5);
            th1emptyclone->GetXaxis()->SetTitle(pv.xlabel);
            th1emptyclone->GetXaxis()->SetTitleOffset(0.9);
            th1emptyclone->GetXaxis()->SetTitleSize(0.15);
            th1emptyclone->GetXaxis()->SetLabelSize(0.15);
            th1emptyclone->SetYTitle("Data/MC");
            th1emptyclone->GetYaxis()->SetTitleSize(0.13);
            th1emptyclone->GetXaxis()->SetNdivisions(505);
            th1emptyclone->GetYaxis()->SetNdivisions(505);
            th1emptyclone->GetYaxis()->SetTitleOffset(0.5);
            th1emptyclone->GetYaxis()->CenterTitle(true);
            th1emptyclone->GetYaxis()->SetLabelSize(0.1);
            th1emptyclone->Draw();

            // TBox *errbox = new TBox(pv.AMINRange,0.974,pv.AMAXRange,1.026); // lumi systematic uncertainty
            errbox = new TBox(pv.AMINRange, 0.974, pv.AMAXRange, 1.026); // lumi systematic uncertainty
            errbox->SetFillColor(kGray);
            errbox->Draw();

            hhratio->Draw("esame");
            TLine *line1;
            line1 = new TLine(pv.AMINRange, 0.5, pv.AMAXRange, 0.5);
            line1->SetLineStyle(2);
            line1->SetLineWidth(1);
            line1->SetLineColor(1);
            line1->Draw();
            TLine *line2;
            line2 = new TLine(pv.AMINRange, 1.0, pv.AMAXRange, 1.0);
            line2->SetLineStyle(2);
            line2->SetLineWidth(1);
            line2->SetLineColor(1);
            line2->Draw();
            TLine *line3;
            line3 = new TLine(pv.AMINRange, 1.5, pv.AMAXRange, 1.5);
            line3->SetLineStyle(2);
            line3->SetLineWidth(1);
            line3->SetLineColor(1);
            line3->Draw();
        }

        c1->Print(outfile + ".pdf");
        //c1->Print(outfile + ".png");
        c1->Print(outfile + ".root");
        c1->Print(outfile + ".C");
        // c1->Print(outfile+".tex");
        //-----------------------------------------------------------------
        //	Save log plot
        //-----------------------------------------------------------------
        TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
        TPad *d3, *d4;

        c2->SetFillColor(0);
        c2->SetBorderMode(0);
        c2->SetBorderSize(10);
        // Set margins to reasonable defaults
        c2->SetLeftMargin(0.18);
        c2->SetRightMargin(0.05);
        c2->SetTopMargin(0.08);
        c2->SetBottomMargin(0.15);
        // Setup a frame which makes sense
        c2->SetFrameFillStyle(0);
        c2->SetFrameLineStyle(0);
        c2->SetFrameBorderMode(0);
        c2->SetFrameBorderSize(10);
        c2->SetFrameFillStyle(0);
        c2->SetFrameLineStyle(0);
        c2->SetFrameBorderMode(0);
        c2->SetFrameBorderSize(10);

        c2->Divide(1, 2);
        d3 = (TPad *)c2->GetPad(1);
        d3->SetPad(0.01, 0.30, 0.95, 0.99);
        // c1->cd();
        d4 = (TPad *)c2->GetPad(2);
        d4->SetPad(0.01, 0.02, 0.95, 0.30);
        d3->cd();
        gPad->SetBottomMargin(0.0);

        gPad->SetTopMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetLeftMargin(0.14);

        gPad->SetLogy(1); 

        th1totempty->SetMaximum(800 * maxval);
        th1totempty->SetMinimum(0.001);
        th1totempty->Draw();

        if (th1data)
        {
            th1data->SetMinimum(0.001);
            th1data->SetBinErrorOption(TH1D::kPoisson);
            th1data->Draw("e0same");
        }

        th1tot->Draw("e2same");

        hs->SetMinimum(0.001);
        hs->Draw("samehist");
        if (pv.drawleg == 1)
            Leg->Draw();

        th1tot->Draw("e2same");
        if (h1->GetEntries() > 0)
            h1->Draw("histsame");
        if (h2->GetEntries() > 0)
            h2->Draw("histsame");
        if (th1data)
            th1data->Draw("e1same");

        CMS_lumi(d3, 4, 10);
        gPad->RedrawAxis();
        d4->cd();
        gPad->SetLeftMargin(0.14);
        gPad->SetTopMargin(0);
        gPad->SetRightMargin(0.05);
        gPad->SetFrameBorderSize(0);
        gPad->SetBottomMargin(0.3);
        gPad->SetTickx();
        th1emptyclone->Draw();
        errbox->Draw();

        hhratio->Draw("esame");
        TLine *line;
        line = new TLine(pv.AMINRange, 1.0, pv.AMAXRange, 1.0);
        line->SetLineStyle(1);
        line->SetLineWidth(1);
        line->SetLineColor(1);
        line->Draw();

        c2->Print(outfile + "_log.png");
        //-----------------------------------------------------------------
        //	END:	Save log plot
        //-----------------------------------------------------------------

        Logfile.close();
        c1->Write();
    } // var loop

    // f.Write();
    f->Close();

} // myControlPlots
