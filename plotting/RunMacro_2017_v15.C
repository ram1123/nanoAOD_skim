// 2017 (UL, NanoAOD v15) control plots with the AN-2016/325 event selection (see make_config.py):
// low-MET CR, m_ll-sideband (DY) CR, emu CR (= N_emu^in of the alpha method), alpha-method sidebands,
// and the single-photon CR (separate skim / inputs).
// Run from nanoAOD_skim/plotting:
//   root -l -b -q 'RunMacro_2017_v15.C("LowMET")'      // one region
//   root -l -b -q 'RunMacro_2017_v15.C("all")'         // all dilepton regions
//   root -l -b -q 'RunMacro_2017_v15.C("gammaCR")'     // single-photon CR
#include "controlplotvars_2017.h"

void RunMacro_2017_v15(TString which = "all",
                       TString outBase = "/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325/plots")
{
    const TString dir = ".";
    gROOT->ProcessLine(".L myControlPlots.C+");
    gROOT->ProcessLine("gMCWeight = \"overallEventWeight\";");   // genWeight * puWeight * L1PreFiringWeight_Nom
    gROOT->ProcessLine("lumi_13TeV = \"41.5 fb^{-1}\";");        // CMS_lumi label (default is the 2018 value)

    const std::vector<TString> regions = {"LowMET", "mllSideband", "emu", "alphaSB", "gammaCR"};

    for (const auto &region : regions) {
        if (which != region && !(which == "all" && region != "gammaCR"))
            continue;
        const bool isGamma = (region == "gammaCR");
        // gamma CR: own skim/inputs; data weighted by the prescale of the photon-pT window's trigger
        const TString samples = dir + (isGamma ? "/HZZ2l2nu_Samples13TeV_2017_v15_gamma.txt"
                                               : "/HZZ2l2nu_Samples13TeV_2017_v15.txt");
        gROOT->ProcessLine(isGamma ? "gDataWeight = \"gamma_prescaleWeight\";" : "gDataWeight = \"\";");
        std::vector<TString> flavs = {"ee", "mumu"};
        std::vector<TString> jetcats = {"0jet", "ge1jet", "vbf"};
        if (region == "emu" || isGamma)
            flavs = {""};
        if (region == "alphaSB") {
            flavs = {"ee", "mumu", "emu"};
            jetcats = {"incl"};
        }
        for (const auto &flav : flavs)
            for (const auto &jc : jetcats) {
                TString tag = region + "_" + jc + (flav.Length() ? "_" + flav : TString(""));
                TString cuts = dir + "/cuts/HZZ2l2nu_2017_Cuts_" + tag + ".txt";
                TString outdir = outBase + "/" + region + "/" + tag;
                gSystem->mkdir(outdir, kTRUE);
                std::cout << "\n===== " << tag << " =====" << std::endl;
                // per-region ranges/binning; no jet variables in the 0-jet category, VBF variables in VBF
                TString kind = (jc == "0jet") ? "_0jet" : ((jc == "vbf" && !isGamma) ? "_vbf" : "_jets");
                TString vars = "controlplotvars_2017_" + region + kind;
                gROOT->ProcessLine(Form("myControlPlots(\"%s\", \"%s\", %s, \"%s/%s.root\", \"%s\", 0)",
                                        cuts.Data(), samples.Data(), vars.Data(), outdir.Data(), tag.Data(), outdir.Data()));
            }
    }
}
