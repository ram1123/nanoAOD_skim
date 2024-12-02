#include <TLatex.h>
// #include "RooRealVar.h"
// #include "RooDataSet.h"
#include <TString.h>
// #include "RooFormulaVar.h"
// #include "RooAddPdf.h"
// #include "RooExponential.h"
// #include "RooGaussModel.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
// #include "RooHist.h"
// #include "RooGenericPdf.h"
// #include "RooTruthModel.h"
// #include "RooGaussModel.h"
// #include "RooDecay.h"
// #include "RooProdPdf.h"
// #include "RooEffProd.h"
#include "TF1.h"
#include "TH1F.h"

using namespace std;
void dylep1pt()
{

    TFile *f1 = new TFile("/eos/user/a/avijay/DYsamples/dy_0to50/dy_0to50.root");
    TFile *f2 = new TFile("/eos/user/a/avijay/DYsamples/dy_50to100/dy_50to100.root");
    TFile *f3 = new TFile("/eos/user/a/avijay/DYsamples/dy_100to250/dy_100to250.root");
    TFile *f4 = new TFile("/eos/user/a/avijay/DYsamples/dy_250to400/dy_250to400.root");
    TFile *f5 = new TFile("/eos/user/a/avijay/DYsamples/dy_400to650/dy_400to650.root");
    TFile *f6 = new TFile("/eos/user/a/avijay/DYsamples/dy_650toInf/dy_650toInf.root");

    TTree *tree1 = (TTree *)f1->Get("Events");
    TTree *tree2 = (TTree *)f2->Get("Events");
    TTree *tree3 = (TTree *)f3->Get("Events");
    TTree *tree4 = (TTree *)f4->Get("Events");
    TTree *tree5 = (TTree *)f5->Get("Events");
    TTree *tree6 = (TTree *)f6->Get("Events");

    // Loop for pTL1
    Float_t pTL11;
    Float_t pTL12;
    Float_t pTL13;
    Float_t pTL14;
    Float_t pTL15;
    Float_t pTL16;

    Int_t m1_nentries_ = tree1->GetEntries();
    Int_t m2_nentries_ = tree2->GetEntries();
    Int_t m3_nentries_ = tree3->GetEntries();
    Int_t m4_nentries_ = tree4->GetEntries();
    Int_t m5_nentries_ = tree5->GetEntries();
    Int_t m6_nentries_ = tree6->GetEntries();

    tree1->SetBranchAddress("pTL1", &pTL11);
    tree2->SetBranchAddress("pTL1", &pTL12);
    tree3->SetBranchAddress("pTL1", &pTL13);
    tree4->SetBranchAddress("pTL1", &pTL14);
    tree5->SetBranchAddress("pTL1", &pTL15);
    tree6->SetBranchAddress("pTL1", &pTL16);

    TH1F *h1 = new TH1F("h1", "", 40, 0, 400);
    TH1F *h2 = new TH1F("h2", "", 40, 0, 400);
    TH1F *h3 = new TH1F("h3", "", 40, 0, 400);
    TH1F *h4 = new TH1F("h4", "", 40, 0, 400);
    TH1F *h5 = new TH1F("h5", "", 40, 0, 400);
    TH1F *h6 = new TH1F("h6", "", 40, 0, 400);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas *c1 = new TCanvas("c1", "hists with different scales", 600, 400);
    TPad *gpad1 = new TPad("gpad1", "", 0, 0, 1, 1);
    gpad1->Draw();
    gpad1->cd();
    for (int j = 0; j < m1_nentries_; j++)
    {
        tree1->GetEntry(j);
        h1->Fill(pTL11);
        //      double content = h1->GetBinContent(j);
        //      h1->SetBinContent(j, content * 0.81);
    }
    for (int j = 0; j < m2_nentries_; j++)
    {
        tree2->GetEntry(j);
        h2->Fill(pTL12);
        //      double content1 = h2->GetBinContent(j);
        //      h2->SetBinContent(j, content1 * 0.37);
    }
    for (int j = 0; j < m3_nentries_; j++)
    {
        tree3->GetEntry(j);
        h3->Fill(pTL13);
        //      double content2 = h3->GetBinContent(j);
        //      h3->SetBinContent(j, content2 * 0.39);
    }
    for (int j = 0; j < m4_nentries_; j++)
    {
        tree4->GetEntry(j);
        h4->Fill(pTL14);
        //      double content3 = h4->GetBinContent(j);
        //      h4->SetBinContent(j, content3 * 0.04);
    }
    for (int j = 0; j < m5_nentries_; j++)
    {
        tree5->GetEntry(j);
        h5->Fill(pTL15);
        //      double content4 = h5->GetBinContent(j);
        //       h5->SetBinContent(j, content4 * 0.015);
    }
    for (int j = 0; j < m6_nentries_; j++)
    {
        tree6->GetEntry(j);
        h6->Fill(pTL16);
        //      double content5 = h6->GetBinContent(j);
        //      h6->SetBinContent(j, content5 * 0.0011);
    }

    h1->Scale(0.81);
    h2->Scale(0.37);
    h3->Scale(0.39);
    h4->Scale(0.04);
    h5->Scale(0.015);
    h6->Scale(0.0011);

    h1->Add(h2);
    h1->Add(h3);
    h1->Add(h4);
    h1->Add(h5);
    h1->Add(h6);
    h1->SetLineColor(kRed);
    h1->SetLineWidth(3);
    h2->SetLineColor(kViolet);
    h2->SetLineWidth(3);
    h3->SetLineColor(kCyan);
    h3->SetLineWidth(3);
    h4->SetLineColor(kBlue);
    h4->SetLineWidth(3);
    h5->SetLineColor(kGreen);
    h5->SetLineWidth(3);
    h6->SetLineColor(kBlack);
    h6->SetLineWidth(3);

    // hist6->GetXaxis()->SetTitle("P_{T} of leading lepton GeV/c");
    // hist6->GetXaxis()->SetLabelSize(0.045);
    // hist6->GetYaxis()->SetTitle("Events / 10 GeV/c");
    // hist6->GetYaxis()->SetLabelSize(0.045);
    // h1->Scale(0.81);

    h1->Integral();
    h1->Draw();
    // h2->Draw("same");
    // h3->Draw("same");
    // h4->Draw("same");
    // h5->Draw("same");
    // h6->Draw("same");

    /*TLegend *pl = new TLegend(0.3,0.4,0.6,0.6);
    pl->SetTextSize(0.03);
    pl->SetFillColor(0);

    TLegendEntry *plr = pl->AddEntry(hist1, "DrellYan_Pt_0to50",  "l");
    TLegendEntry *pls = pl->AddEntry(hist2, "DrellYan_Pt_50to100",  "l");
    TLegendEntry *plv = pl->AddEntry(hist3, "DrellYan_Pt_100to250",  "l");
    TLegendEntry *plw = pl->AddEntry(hist4, "DrellYan_Pt_250to400",  "l");
    TLegendEntry *plx = pl->AddEntry(hist5, "DrellYan_Pt_400to650",  "l");
    TLegendEntry *ply = pl->AddEntry(hist6, "DrellYan_Pt_650toInf",  "l");


    plr->SetMarkerSize(3.0);
    plr->SetMarkerColor(kRed);
    plr->SetTextSize(0.045);
    pls->SetMarkerSize(3.0);
    pls->SetMarkerColor(kViolet);
    pls->SetTextSize(0.045);
    plv->SetMarkerSize(3.0);
    plv->SetMarkerColor(kCyan);
    plv->SetTextSize(0.045);
    plw->SetMarkerSize(3.0);
    plw->SetMarkerColor(kBlue);
    plw->SetTextSize(0.045);
    plx->SetMarkerSize(3.0);
    plx->SetMarkerColor(kGreen);
    plx->SetTextSize(0.045);
    ply->SetMarkerSize(3.0);
    ply->SetMarkerColor(kBlack);
    ply->SetTextSize(0.045);

    pl->Draw();
*/
}
