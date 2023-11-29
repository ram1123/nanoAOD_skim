#ifndef H4LTools_h
#define H4LTools_h

#include <utility>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TSpline.h>
#include <vector>
#include "yaml-cpp/yaml.h"
#include "RoccoR.h"
#include "../JHUGenMELA/MELA/interface/Mela.h"

class H4LTools {
    public:
      H4LTools(int year, std::string DATAPATH);

      float elePtcut, MuPtcut, eleEtacut, MuEtacut, elesip3dCut, Musip3dCut,Zmass,MZ1cut,MZcutup,MZcutdown,MZZcut,HiggscutUp,HiggscutDown;
      float eleLoosedxycut,eleLoosedzcut,MuLoosedxycut,MuLoosedzcut,MuTightdxycut,MuTightdzcut,MuTightTrackerLayercut,MuTightpTErrorcut,MuHighPtBound,eleIsocut,MuIsocut;
      float fsrphotonPtcut,fsrphotonEtacut,fsrphotonIsocut,fsrphotondRlcut,fsrphotondRlOverPtcut, JetPtcut,JetEtacut;
      float eleBDTWPLELP,eleBDTWPMELP,eleBDTWPHELP,eleBDTWPLEHP,eleBDTWPMEHP,eleBDTWPHEHP;
      float HZZ2l2q_Leading_Lep_pT, HZZ2l2q_SubLeading_Lep_pT, HZZ2l2q_Lep_eta, HZZ2l2q_MZLepcutdown, HZZ2l2q_MZLepcutup;

      void InitializeElecut(float elePtcut_,float eleEtacut_,float elesip3dCut_,float eleLoosedxycut_,float eleLoosedzcut_,float eleIsocut_,float eleBDTWPLELP_,float eleBDTWPMELP_, float eleBDTWPHELP_,float eleBDTWPLEHP_,float eleBDTWPMEHP_,float eleBDTWPHEHP_){
        elePtcut = elePtcut_;
        eleEtacut = eleEtacut_;
        elesip3dCut = elesip3dCut_;
        eleLoosedxycut = eleLoosedxycut_;
        eleLoosedzcut = eleLoosedzcut_;
        eleIsocut = eleIsocut_;
        eleBDTWPLELP = eleBDTWPLELP_;
        eleBDTWPMELP = eleBDTWPMELP_;
        eleBDTWPHELP = eleBDTWPHELP_;
        eleBDTWPLEHP = eleBDTWPLEHP_;
        eleBDTWPMEHP = eleBDTWPMEHP_;
        eleBDTWPHEHP = eleBDTWPHEHP_;
      }

      void Initialize2l2qEvtCut(float HZZ2l2q_Leading_Lep_pT_, float HZZ2l2q_SubLeading_Lep_pT_, float HZZ2l2q_Lep_eta_, float HZZ2l2q_MZLepcutdown_, float HZZ2l2q_MZLepcutup_) {
        HZZ2l2q_Leading_Lep_pT = HZZ2l2q_Leading_Lep_pT_;
        HZZ2l2q_SubLeading_Lep_pT = HZZ2l2q_SubLeading_Lep_pT_;
        HZZ2l2q_Lep_eta = HZZ2l2q_Lep_eta_;
        HZZ2l2q_MZLepcutdown = HZZ2l2q_MZLepcutdown_;
        HZZ2l2q_MZLepcutup = HZZ2l2q_MZLepcutup_;
      }

      void InitializeMucut(float MuPtcut_, float MuEtacut_, float Musip3dCut_, float MuLoosedxycut_, float MuLoosedzcut_, float MuIsocut_, float MuTightdxycut_, float MuTightdzcut_, float MuTightTrackerLayercut_, float MuTightpTErrorcut_, float MuHighPtBound_)
      {
        MuPtcut = MuPtcut_;
        MuEtacut = MuEtacut_;
        Musip3dCut = Musip3dCut_;
        MuLoosedxycut = MuLoosedxycut_;
        MuLoosedzcut = MuLoosedzcut_;
        MuIsocut = MuIsocut_;
        MuTightdxycut = MuTightdxycut_;
        MuTightdzcut = MuTightdzcut_;
        MuTightTrackerLayercut = MuTightTrackerLayercut_;
        MuTightpTErrorcut = MuTightpTErrorcut_;
        MuHighPtBound = MuHighPtBound_;
      }
      void InitializeFsrPhotonCut(float fsrphotonPtcut_, float fsrphotonEtacut_, float fsrphotonIsocut_, float fsrphotondRlcut_, float fsrphotondRlOverPtcut_){
        fsrphotonPtcut = fsrphotonPtcut_;
        fsrphotonEtacut = fsrphotonEtacut_;
        fsrphotonIsocut = fsrphotonIsocut_;
        fsrphotondRlcut = fsrphotondRlcut_;
        fsrphotondRlOverPtcut = fsrphotondRlOverPtcut_;
      }
      void InitializeJetcut(float JetPtcut_, float JetEtacut_){
        JetPtcut = JetPtcut_;
        JetEtacut = JetEtacut_;
      }
      void InitializeEvtCut(float MZ1cut_,float MZZcut_,float HiggscutDown_,float HiggscutUp_,float Zmass_,float MZcutdown_, float MZcutup_){
        MZ1cut = MZ1cut_;
        MZZcut = MZZcut_;
        HiggscutDown = HiggscutDown_;
        HiggscutUp = HiggscutUp_;
        Zmass = Zmass_;
        MZcutdown = MZcutdown_;
        MZcutup = MZcutup_;
      }

      void SetElectrons(float Electron_pt_, float Electron_eta_, float Electron_phi_, float Electron_mass_, float Electron_dxy_,float Electron_dz_,
                        float Electron_sip3d_, float Electron_mvaFall17V2Iso_, int Electron_pdgId_, float Electron_pfRelIso03_all_){
        Electron_pt.push_back(Electron_pt_);
        Electron_phi.push_back(Electron_phi_);
        Electron_eta.push_back(Electron_eta_);
        Electron_mass.push_back(Electron_mass_);
        Electron_dxy.push_back(Electron_dxy_);
        Electron_dz.push_back(Electron_dz_);
        Electron_sip3d.push_back(Electron_sip3d_);
        Electron_mvaFall17V2Iso.push_back(Electron_mvaFall17V2Iso_);
        Electron_pdgId.push_back(Electron_pdgId_);
        Electron_pfRelIso03_all.push_back(Electron_pfRelIso03_all_);
      }

      void SetJets(float Jet_pt_, float Jet_eta_, float Jet_phi_, float Jet_mass_, int Jet_jetId_, float Jet_btagDeepC_,
                         int Jet_puId_){
        Jet_pt.push_back(Jet_pt_);
        Jet_phi.push_back(Jet_phi_);
        Jet_eta.push_back(Jet_eta_);
        Jet_mass.push_back(Jet_mass_);
        Jet_btagDeepC.push_back(Jet_btagDeepC_);
        Jet_jetId.push_back(Jet_jetId_);
        Jet_puId.push_back(Jet_puId_); //1 or 0?
      }

      void SetFatJets(float Jet_pt_, float Jet_eta_, float Jet_phi_, float Jet_mass_, int Jet_jetId_, float Jet_btagDeepB_,
                      float Jet_PNZvsQCD_)
      {
        FatJet_pt.push_back(Jet_pt_);
        FatJet_eta.push_back(Jet_eta_);
        FatJet_phi.push_back(Jet_phi_);
        FatJet_SDmass.push_back(Jet_mass_);
        FatJet_jetId.push_back(Jet_jetId_);
        FatJet_btagDeepB.push_back(Jet_btagDeepB_);
        FatJet_PNZvsQCD.push_back(Jet_PNZvsQCD_); // 1 or 0?
      }

      void SetMET(float MET_pt_, float MET_phi_, float MET_sumEt_)
      {
        MET_pt = MET_pt_;
        MET_phi = MET_phi_;
        MET_sumEt = MET_sumEt_;
	std::cout<<"Inside header file: MET_sumEt = " << MET_sumEt_ << "\t" << MET_sumEt << std::endl;
      }

      
      void SetMuons(float Muon_pt_, float Muon_eta_, float Muon_phi_, float Muon_mass_, bool Muon_isGlobal_, bool Muon_isTracker_,
                        float Muon_dxy_, float Muon_dz_,float Muon_sip3d_, float Muon_ptErr_,
                        int Muon_nTrackerLayers_, bool Muon_isPFcand_, int Muon_pdgId_,int Muon_charge_, float Muon_pfRelIso03_all_
                        ){
        Muon_pt.push_back(Muon_pt_); 
        Muon_phi.push_back(Muon_phi_);
        Muon_eta.push_back(Muon_eta_);
        Muon_mass.push_back(Muon_mass_);
        Muon_isGlobal.push_back(Muon_isGlobal_);
        Muon_isTracker.push_back(Muon_isTracker_);
        Muon_dxy.push_back(Muon_dxy_);
        Muon_dz.push_back(Muon_dz_);
        Muon_sip3d.push_back(Muon_sip3d_);
        Muon_ptErr.push_back(Muon_ptErr_);
        Muon_nTrackerLayers.push_back(Muon_nTrackerLayers_);
        Muon_isPFcand.push_back(Muon_isPFcand_);
        Muon_pdgId.push_back(Muon_pdgId_);
        Muon_charge.push_back(Muon_charge_);
        Muon_pfRelIso03_all.push_back(Muon_pfRelIso03_all_);

      }
      void SetMuonsGen(int Muon_genPartIdx_){
        Muon_genPartIdx.push_back(Muon_genPartIdx_);
      }
      /*void SetMuons(TTreeReaderArray<float> *Muon_pt_, TTreeReaderArray<float> *Muon_eta_,
                        TTreeReaderArray<float> *Muon_phi_, TTreeReaderArray<float> *Muon_mass_, TTreeReaderArray<bool> *Muon_isGlobal_, TTreeReaderArray<bool> *Muon_isTracker_,
                        TTreeReaderArray<float> *Muon_dxy_, TTreeReaderArray<float> *Muon_dz_,TTreeReaderArray<float> *Muon_sip3d_, TTreeReaderArray<float> *Muon_ptErr_,
                        TTreeReaderArray<int> *Muon_nTrackerLayers_, TTreeReaderArray<bool> *Muon_isPFcand_, TTreeReaderArray<int> *Muon_pdgId_,TTreeReaderArray<int> *Muon_charge_, TTreeReaderArray<float> *Muon_pfRelIso03_all_,
                        TTreeReaderArray<int> *Muon_genPartIdx_){
        Muon_pt = Muon_pt_;
        Muon_phi = Muon_phi_;
        Muon_eta = Muon_eta_;
        Muon_mass = Muon_mass_;
        Muon_isGlobal = Muon_isGlobal_;
        Muon_isTracker = Muon_isTracker_;
        Muon_dxy = Muon_dxy_;
        Muon_dz = Muon_dz_;
        Muon_sip3d = Muon_sip3d_;
        Muon_ptErr = Muon_ptErr_;
        Muon_nTrackerLayers = Muon_nTrackerLayers_;
        Muon_isPFcand = Muon_isPFcand_;
        Muon_pdgId = Muon_pdgId_;
        Muon_charge = Muon_charge_;
        Muon_pfRelIso03_all = Muon_pfRelIso03_all_;
        Muon_genPartIdx = Muon_genPartIdx_;
      }*/
      void SetFsrPhotons(float FsrPhoton_dROverEt2_, float FsrPhoton_eta_,
                        float FsrPhoton_phi_, float FsrPhoton_pt_, float FsrPhoton_relIso03_){
        FsrPhoton_dROverEt2.push_back(FsrPhoton_dROverEt2_);
        FsrPhoton_phi.push_back(FsrPhoton_phi_);
        FsrPhoton_eta.push_back(FsrPhoton_eta_);
        FsrPhoton_pt.push_back(FsrPhoton_pt_);
        FsrPhoton_relIso03.push_back(FsrPhoton_relIso03_);
      }
      /*void SetFsrPhotons(TTreeReaderArray<float> *FsrPhoton_dROverEt2_, TTreeReaderArray<float> *FsrPhoton_eta_,
                        TTreeReaderArray<float> *FsrPhoton_phi_, TTreeReaderArray<float> *FsrPhoton_pt_,
                        TTreeReaderArray<float> *FsrPhoton_relIso03_){
        FsrPhoton_dROverEt2 = FsrPhoton_dROverEt2_;
        FsrPhoton_phi = FsrPhoton_phi_;
        FsrPhoton_eta = FsrPhoton_eta_;
        FsrPhoton_pt = FsrPhoton_pt_;
        FsrPhoton_relIso03 = FsrPhoton_relIso03_;
      }*/
      void SetGenParts(float GenPart_pt_){
        GenPart_pt.push_back(GenPart_pt_);
      }
      /*void SetGenParts(TTreeReaderArray<float> *GenPart_pt_){
        GenPart_pt = GenPart_pt_;
      }*/
      void SetObjectNum(unsigned nElectron_,unsigned nMuon_,unsigned nJet_,unsigned nFsrPhoton_){
        nElectron = nElectron_;
        nMuon = nMuon_;
        nJet = nJet_;
        nFsrPhoton = nFsrPhoton_;
      }
      void SetObjectNumGen(unsigned nGenPart_){
        nGenPart = nGenPart_;
      }

      std::vector<unsigned int> goodLooseElectrons2012();
      std::vector<unsigned int> goodLooseMuons2012();
      std::vector<unsigned int> goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex);
      std::vector<unsigned int> goodElectrons2015_noIso_noBdt(std::vector<unsigned int> Electronindex);
      std::vector<bool> passTight_BDT_Id();
      std::vector<bool> passTight_Id();
      std::vector<unsigned int> goodFsrPhotons();
      unsigned doFsrRecovery(TLorentzVector Lep);
      std::vector<TLorentzVector> BatchFsrRecovery(std::vector<TLorentzVector> LepList);
      std::vector<TLorentzVector> ElectronFsr();
      std::vector<TLorentzVector> MuonFsr();
      std::vector<float> ElectronFsrPt();
      std::vector<float> ElectronFsrEta();
      std::vector<float> ElectronFsrPhi();
      std::vector<float> MuonFsrPt();
      std::vector<float> MuonFsrEta();
      std::vector<float> MuonFsrPhi();
      std::vector<unsigned int> SelectedJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu);
      std::vector<unsigned int> SelectedFatJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu);

      std::vector<TLorentzVector> Zlist;
      std::vector<TLorentzVector> Zlistnofsr;
      std::vector<int> Zflavor; //mu->13, e->11
      std::vector<int> Zlep1index;
      std::vector<int> Zlep2index;
      std::vector<float> Zlep1pt;
      std::vector<float> Zlep1eta;
      std::vector<float> Zlep1phi;
      std::vector<float> Zlep1mass;
      std::vector<float> Zlep1chg;
      std::vector<float> Zlep2pt;
      std::vector<float> Zlep2eta;
      std::vector<float> Zlep2phi;
      std::vector<float> Zlep2mass;
      std::vector<float> Zlep2chg;
      std::vector<float> Zlep1ptNoFsr;
      std::vector<float> Zlep1etaNoFsr;
      std::vector<float> Zlep1phiNoFsr;
      std::vector<float> Zlep1massNoFsr;
      std::vector<float> Zlep2ptNoFsr;
      std::vector<float> Zlep2etaNoFsr;
      std::vector<float> Zlep2phiNoFsr;
      std::vector<float> Zlep2massNoFsr;
      std::vector<unsigned int> jetidx;
      std::vector<unsigned int> FatJetidx;

      int nTightEle;
      int nTightMu;
      int nTightEleChgSum;
      int nTightMuChgSum;
      bool flag4e;
      bool flag4mu;
      bool flag2e2mu;

      bool flag2e;
      bool flag2mu;
      bool flag2l;
      bool flag2e_met;
      bool flag2mu_met;
      bool flag2l_met; 

      void LeptonSelection();
      std::vector<unsigned int> looseEle,looseMu,bestEle,bestMu, tighteleforjetidx, tightmuforjetidx;
      std::vector<unsigned int> Electronindex;
      std::vector<unsigned int> Muonindex;
      std::vector<bool> AllEid;
      std::vector<bool> AllMuid;
      std::vector<TLorentzVector> Elelist;
      std::vector<TLorentzVector> Mulist;
      std::vector<TLorentzVector> ElelistFsr;
      std::vector<TLorentzVector> MulistFsr;
      std::vector<int> Elechg;
      std::vector<int> Muchg;
      std::vector<float> Muiso,Eiso;
      std::vector<bool> Eid;
      std::vector<bool> muid;

      std::vector<int> TightEleindex;
      std::vector<int> TightMuindex;
      void Initialize(){
        Electron_pt.clear();Electron_phi.clear();Electron_eta.clear();Electron_mass.clear();Electron_dxy.clear();Electron_dz.clear();Electron_sip3d.clear();
        Electron_mvaFall17V2Iso.clear();Electron_pdgId.clear();Electron_pfRelIso03_all.clear();
        Muon_pt.clear();Muon_phi.clear();Muon_eta.clear();Muon_mass.clear();Muon_dxy.clear();Muon_dz.clear();Muon_sip3d.clear();Muon_ptErr.clear();Muon_pfRelIso03_all.clear();
        Muon_nTrackerLayers.clear();Muon_genPartIdx.clear();Muon_pdgId.clear();Muon_charge.clear();
        Muon_isTracker.clear();Muon_isGlobal.clear();Muon_isPFcand.clear();
        Jet_pt.clear();Jet_phi.clear();Jet_eta.clear();Jet_mass.clear();Jet_btagDeepC.clear();
        MET_pt = 0.0; MET_phi = 0.0;  ////new
	MET_sumEt = 0.0;

        Jet_jetId.clear();Jet_puId.clear();
        FatJet_pt.clear();FatJet_phi.clear();FatJet_eta.clear();FatJet_SDmass.clear();FatJet_btagDeepB.clear(); FatJet_PNZvsQCD.clear();

        FatJet_jetId.clear();
        FsrPhoton_dROverEt2.clear();FsrPhoton_phi.clear();FsrPhoton_eta.clear();FsrPhoton_pt.clear();FsrPhoton_relIso03.clear();
        GenPart_pt.clear();
        Zlist.clear();
        Zlistnofsr.clear();
        Zflavor.clear();
        Zlep1index.clear();
        Zlep2index.clear();
        Zlep1pt.clear(); Zlep1eta.clear(); Zlep1phi.clear(); Zlep1mass.clear();
        Zlep2pt.clear(); Zlep2eta.clear(); Zlep2phi.clear(); Zlep2mass.clear();
        Zlep1chg.clear(); Zlep2chg.clear();
        Zlep1ptNoFsr.clear(); Zlep1etaNoFsr.clear(); Zlep1phiNoFsr.clear(); Zlep1massNoFsr.clear();
        Zlep2ptNoFsr.clear(); Zlep2etaNoFsr.clear(); Zlep2phiNoFsr.clear(); Zlep2massNoFsr.clear();
        jetidx.clear();
        FatJetidx.clear();

        Muon_Pt_Corrected.clear();
        looseEle.clear(); looseMu.clear(); bestEle.clear(); bestMu.clear();  tighteleforjetidx.clear();  tightmuforjetidx.clear();
        Electronindex.clear();  Muonindex.clear(); AllEid.clear(); AllMuid.clear(); Elelist.clear(); Mulist.clear(); ElelistFsr.clear(); Mulist.clear();
        Elechg.clear(); Muchg.clear(); Muiso.clear();Eiso.clear(); Eid.clear(); muid.clear(); TightEleindex.clear(); TightMuindex.clear();
        nElectron = 0; nMuon = 0; nJet = 0; nFsrPhoton = 0; nGenPart = 0;
        nTightEle = 0; nTightMu = 0; nTightEleChgSum = 0; nTightMuChgSum = 0; 
        pTL1 = -999; etaL1 = -999; phiL1 = -999; massL1 = -999;
        pTL2 = -999; etaL2 = -999; phiL2 = -999; massL2 = -999;
        pTL3 = -999; etaL3 = -999; phiL3 = -999; massL3 = -999;
        pTL4 = -999; etaL4 = -999; phiL4 = -999; massL4 = -999;

        pTj1 = -99;  etaj1 = -99;  phij1 = -99;  mj1 = -99;
        pTj2 = -99;  etaj2 = -99;  phij2 = -99;  mj2 = -99;

        flag4e=false; flag4mu=false; flag2e2mu=false;
        flag2e=false; flag2mu=false; flag2l=false; flag2e_met=false; flag2l_met=false; flag2mu_met=false;
      }
      bool isFSR=true;
      unsigned int Zsize=0;
      TSpline *spline_g4;
      TSpline *spline_g2;
      TSpline *spline_L1;
      TSpline *spline_L1Zgs;
      bool findZCandidate();
      bool ZZSelection_4l();
      bool ZZSelection_2l2q();
      bool ZZSelection_2l2nu();
      TLorentzVector Z1;
      TLorentzVector Z1nofsr;
      TLorentzVector Z2;
      TLorentzVector Z2_2j;
      TLorentzVector Z2_met;
      TLorentzVector Z2nofsr;
      TLorentzVector ZZsystem;
      TLorentzVector ZZsystemnofsr;
      TLorentzVector ZZ_2jsystem;
      TLorentzVector ZZ_metsystem;
      TLorentzVector ZZ_2jsystemnofsr;
      TLorentzVector ZZ_metsystemnofsr;

      RoccoR  *calibrator;
      Mela* mela;
      float ApplyRoccoR(bool isMC, int charge, float pt, float eta, float phi, float genPt, float nLayers);
      std::vector<float> Muon_Pt_Corrected;
      void MuonPtCorrection(bool isMC);
      float me_0plus_JHU, me_qqZZ_MCFM, p0plus_m4l, bkg_m4l;
      float D_bkg_kin, D_bkg, D_g4, D_g1g4, D_0m, D_CP, D_0hp, D_int, D_L1, D_L1_int, D_L1Zg, D_L1Zgint;
      float D_bkg_kin_vtx_BS;
      float p0minus_VAJHU, Dgg10_VAMCFM, pg1g4_VAJHU;
      float p0plus_VAJHU, p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, pDL1_VAJHU, pD_L1Zgint; //, p0plus_VAJHU;
      float getDg4Constant(float ZZMass);
      float getDg2Constant(float ZZMass);
      float getDL1Constant(float ZZMass);
      float getDL1ZgsConstant(float ZZMass);

      int cut2e_m40_180, cut2mu_m40_180, cut2l_m40_180;
      int cutMETlt150;
      int cutMETgt150;
      int cut2l_met_m40_180, cut2e_met_m40_180, cut2mu_met_m40_180;
      int cut2e, cut2mu, cut2l, cut2l1J, cut2l2j, cut2l1Jor2j, cut2l1met;
      int cut2e_met, cut2mu_met, cut2l_met;
      int cut4e, cut4mu, cut2e2mu, cutZZ4e, cutZZ4mu, cutZZ2e2mu, cutm4l4e, cutm4l4mu, cutm4l2e2mu, cutghost2e2mu, cutQCD2e2mu, cutLepPt2e2mu, cutghost4e, cutQCD4e, cutLepPt4e, cutghost4mu, cutQCD4mu, cutLepPt4mu;
      float pTL1, etaL1, phiL1, massL1, pTL2, etaL2, phiL2, massL2, pTL3, etaL3, phiL3, massL3, pTL4, etaL4, phiL4, massL4;
      float pTj1, etaj1, phij1, mj1, pTj2, etaj2, phij2, mj2;



    private:
      std::vector<float> Electron_pt,Electron_phi,Electron_eta,Electron_mass,Electron_dxy,Electron_dz,Electron_sip3d;
      std::vector<float> Electron_mvaFall17V2Iso,Electron_pfRelIso03_all;
      std::vector<int> Electron_pdgId;

      std::vector<float> Jet_pt,Jet_phi,Jet_eta,Jet_mass,Jet_btagDeepC;
      std::vector<int> Jet_jetId,Jet_puId;
      float MET_pt, MET_phi;
      float MET_sumEt; 

      std::vector<float> FatJet_pt, FatJet_phi, FatJet_eta, FatJet_SDmass, FatJet_btagDeepB, FatJet_PNZvsQCD;
      std::vector<int> FatJet_jetId;

      std::vector<float> Muon_pt,Muon_phi,Muon_eta,Muon_mass,Muon_dxy,Muon_dz,Muon_sip3d,Muon_ptErr,Muon_pfRelIso03_all;
      std::vector<int> Muon_nTrackerLayers,Muon_genPartIdx,Muon_pdgId,Muon_charge;
      std::vector<bool> Muon_isTracker,Muon_isGlobal,Muon_isPFcand;

      std::vector<float> FsrPhoton_dROverEt2,FsrPhoton_phi,FsrPhoton_pt,FsrPhoton_relIso03,FsrPhoton_eta;

      std::vector<float> GenPart_pt;


      unsigned nElectron,nMuon,nJet,nGenPart,nFsrPhoton;



};

H4LTools::H4LTools(int year, std::string DATAPATH){
  std::cout<<"year"<<" "<<year<<std::endl;
  calibrator = new RoccoR(DATAPATH);
  mela = new Mela(13.0, 125.0, TVar::SILENT);
  mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  TFile *gConstant_g4 = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_g4.root");
  spline_g4 = (TSpline*) gConstant_g4->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4");
  gConstant_g4->Close();
  delete gConstant_g4;
  TFile *gConstant_g2 = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_g2.root");
  spline_g2 = (TSpline*) gConstant_g2->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2");
  gConstant_g2->Close();
  delete gConstant_g2;
  TFile *gConstant_L1 = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_L1.root");
  spline_L1 = (TSpline*) gConstant_L1->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1");
  gConstant_L1->Close();
  delete gConstant_L1;
  TFile *gConstant_L1Zgs = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_L1Zgs.root");
  spline_L1Zgs = (TSpline*) gConstant_L1Zgs->Get("sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs");
  gConstant_L1Zgs->Close();
  delete gConstant_L1Zgs;

  cut2e2mu = 0;
  cut4e = 0;
  cut4mu = 0;
  cutghost2e2mu = 0;
  cutghost4e = 0;
  cutghost4mu = 0;
  cutLepPt2e2mu = 0;
  cutLepPt4e = 0;
  cutLepPt4mu = 0;
  cutQCD2e2mu = 0;
  cutQCD4e = 0;
  cutQCD4mu = 0;
  cutZZ2e2mu = 0;
  cutZZ4e = 0;
  cutZZ4mu = 0;
  cutm4l2e2mu = 0;
  cutm4l4e = 0;
  cutm4l4mu = 0;
  cutMETlt150 = 0;
  cutMETgt150 = 0;

  cut2e = 0;
  cut2mu = 0;
  cut2l = 0;
  cut2l1J = 0;
  cut2l2j = 0;
  cut2l1Jor2j = 0;
  cut2e_m40_180 = 0;
  cut2mu_m40_180 = 0;
  cut2l_m40_180 = 0;

  cut2e_met = 0;
  cut2mu_met = 0;
  cut2l_met = 0;
  cut2l_met_m40_180 = 0;
  cut2e_met_m40_180 = 0;
  cut2mu_met_m40_180 = 0;
  cut2l1met = 0;    
}
#endif

