#include "../include/H4LTools.h"
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <vector>


std::vector<unsigned int> H4LTools::goodLooseElectrons2012(){
    std::vector<unsigned int> LooseElectronindex;
    for (unsigned int i=0; i<Electron_pt.size(); i++){
        if (DEBUG) std::cout << "Inside goodLooseElectrons2012:: Electron_pt[" << i << "] = " << Electron_pt[i] << std::endl;
        if (Electron_pt[i] < elePtcut) continue;
        if (fabs(Electron_eta[i]) > eleEtacut) continue;
        if ((fabs(Electron_eta[i])<1.4442)||(fabs(Electron_eta[i])>1.5660))
        {
            if (fabs(Electron_eta[i])<eleEtacut){
                LooseElectronindex.push_back(i);
            }
        }
    }
    return LooseElectronindex;
}

std::vector<unsigned int> H4LTools::goodLooseMuons2012(){
    std::vector<unsigned int> LooseMuonindex;
    for (unsigned int i=0; i<Muon_eta.size(); i++){
        if (DEBUG) std::cout << "Inside goodLooseMuons2012:: Muon_pt[" << i << "] = " << Muon_pt[i] << std::endl;
        if (Muon_pt[i] < MuPtcut) continue;
        if (fabs(Muon_eta[i]) > MuEtacut) continue;
        if ((Muon_isGlobal[i]||Muon_isTracker[i]||Muon_isPFcand[i])&&(Muon_mediumId[i])){
            LooseMuonindex.push_back(i);
        }

    }
    return LooseMuonindex;
}

std::vector<unsigned int> H4LTools::goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex)
{
    std::vector<unsigned int> bestMuonindex;
    for (unsigned int i=0; i<Muonindex.size(); i++){
        if (Muon_pt[Muonindex[i]] < MuPtcut) continue;
        if (fabs(Muon_eta[Muonindex[i]]) > MuEtacut) continue;
        if ((Muon_isGlobal[Muonindex[i]]||Muon_isTracker[Muonindex[i]])&&(Muon_mediumId[Muonindex[i]]))
        {
            if ((fabs(Muon_dxy[Muonindex[i]]) < MuLoosedxycut) && (fabs(Muon_dz[Muonindex[i]]) < MuLoosedzcut))
            {
                if (Muon_sip3d[Muonindex[i]] < Musip3dCut && year == 2022)
                {
                    bestMuonindex.push_back(Muonindex[i]);
                } else
                {
                    bestMuonindex.push_back(Muonindex[i]);
                }
            }
        }
    }
    return bestMuonindex;
}
std::vector<unsigned int> H4LTools::goodElectrons2015_noIso_noBdt(std::vector<unsigned int> Electronindex)
{
    std::vector<unsigned int> bestElectronindex;
    for (unsigned int i=0; i<Electronindex.size(); i++)
    {
        if (Electronindex[i] >= Electron_pt.size())
        {
            std::cerr << "ERROR: Electronindex out of bounds! i = " << i
                      << ", Electronindex[i] = " << Electronindex[i]
                      << ", Electron_pt.size() = " << Electron_pt.size() << std::endl;
            continue;
        }
        if ((Electron_pt[Electronindex[i]])>elePtcut)
        {
                if((fabs(Electron_dxy[Electronindex[i]])<eleLoosedxycut)&&(fabs(Electron_dz[Electronindex[i]])<eleLoosedzcut))
                {
                    if (Electron_sip3d[Electronindex[i]] < elesip3dCut && year == 2022)
                    {
                        bestElectronindex.push_back(Electronindex[i]);
                    } else
                    {
                        bestElectronindex.push_back(Electronindex[i]);
                    }
                }
        }
    }

    return bestElectronindex;
}

std::vector<bool> H4LTools::passTight_BDT_Id()
{
    std::vector<bool> tightid;
    float cutVal, mvaVal;
    cutVal = 1000;
    mvaVal = -1;
    for (unsigned int i = 0; i < Electron_pt.size(); i++)
    {
        mvaVal = Electron_mvaFall17V2Iso_WP90[i];
        tightid.push_back(mvaVal);
    }
    return tightid;
}

std::vector<bool> H4LTools::passTight_BDT_Id_ZZ4l()
{
    std::vector<bool> tightid;
    float cutVal, mvaVal;
    cutVal = 1000;
    mvaVal = -1;
    // unsigned nE = (*nElectron).Get()[0];
    for (unsigned int i = 0; i < Electron_pt.size(); i++)
    {
        if (Electron_pt[i] < 10)
        {
            if (fabs(Electron_eta[i]) < 0.8)
                cutVal = eleBDTWPLELP;
            if ((fabs(Electron_eta[i]) >= 0.8) && (fabs(Electron_eta[i]) < 1.479))
                cutVal = eleBDTWPMELP;
            if (fabs(Electron_eta[i]) >= 1.479)
                cutVal = eleBDTWPHELP;
        }
        else
        {
            if (fabs(Electron_eta[i]) < 0.8)
                cutVal = eleBDTWPLEHP;
            if ((fabs(Electron_eta[i]) >= 0.8) && (fabs(Electron_eta[i]) < 1.479))
                cutVal = eleBDTWPMEHP;
            if (fabs(Electron_eta[i]) >= 1.479)
                cutVal = eleBDTWPHEHP;
        }

        mvaVal = Electron_mvaFall17V2Iso[i];
        if (mvaVal > cutVal)
        {
            tightid.push_back(true);
        }
        else
        {
            tightid.push_back(false);
        }
    }

    return tightid;
}

std::vector<bool> H4LTools::passTight_Id(){
    std::vector<bool> tightid;
    //unsigned nMu = (*nMuon).Get()[0];
    for (unsigned int i=0; i<Muon_pt.size(); i++){
        if (Muon_pt[i]<MuHighPtBound){
            tightid.push_back(Muon_isPFcand[i]);
        }
        else{
            tightid.push_back(Muon_isPFcand[i]||(((Muon_ptErr[i]/Muon_pt[i])<MuTightpTErrorcut)&&(fabs(Muon_dxy[i])<MuTightdxycut)&&(fabs(Muon_dz[i])<MuTightdzcut)&&(Muon_nTrackerLayers[i]>MuTightTrackerLayercut)));
        }

    }

    return tightid;
}

std::vector<unsigned int> H4LTools::goodFsrPhotons(){
    std::vector<unsigned int> goodFsrPhoton;
    //unsigned nFsr = (*nFsrPhoton).Get()[0];
    for (unsigned int i=0; i<FsrPhoton_pt.size(); i++){
        if((FsrPhoton_pt[i]>fsrphotonPtcut)&&(fabs(FsrPhoton_eta[i])<fsrphotonEtacut)&&(FsrPhoton_relIso03[i]<fsrphotonIsocut)){
            goodFsrPhoton.push_back(i);
        }
    }
    return goodFsrPhoton;
}

std::vector<unsigned int> H4LTools::SelectedJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu)
{
    std::vector<unsigned int> goodJets;
    for (unsigned int i = 0; i < Jet_pt.size(); i++)
    {
        if ((Jet_pt[i] <= JetPtcut)) continue;
        if (fabs(Jet_eta[i]) >= JetEtacut) continue;
        if (Jet_jetId[i] <= 0) continue;
        if ((Jet_pt[i] < 50) && (Jet_puId[i] != 7)) continue;
        if (DEBUG)
            std::cout << "DEBUG: Jet_pt.size() = " << Jet_pt.size() << ";" << " JetID = " << Jet_jetId[i] << ";" << "Jet_pt = " << Jet_pt[i] << ";" << " puID = " << Jet_puId[i] << std::endl;
        int overlaptag = 0;
        TLorentzVector jettest;
        jettest.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        for (unsigned int ie = 0; ie < ele.size(); ie++)
        {
            TLorentzVector eletest;
            eletest.SetPtEtaPhiM(Electron_pt[ele[ie]], Electron_eta[ele[ie]], Electron_phi[ele[ie]], Electron_mass[ele[ie]]);
            if (eletest.DeltaR(jettest) < 0.5)
                overlaptag++;
        }
        for (unsigned int im = 0; im < mu.size(); im++)
        {
            TLorentzVector mutest;
            mutest.SetPtEtaPhiM(Muon_pt[mu[im]], Muon_eta[mu[im]], Muon_phi[mu[im]], Muon_mass[mu[im]]);
            if (mutest.DeltaR(jettest) < 0.5)
                overlaptag++;
        }
        if (overlaptag == 0)
            goodJets.push_back(i);
    }
    njets_pt30_eta4p7 = goodJets.size();

    // count the number of tight, medium and loose b-tagged jets
    // Reference: https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/#ak4-b-tagging
    for (unsigned int i = 0; i < goodJets.size(); i++) // FIXME: These variables seems to be wrong.
    {
        if (DEBUG)
            std::cout << "Jet_btagDeepFlavB[" << goodJets[i] << "]: " << Jet_btagDeepFlavB[goodJets[i]] << std::endl;

        if (Jet_btagDeepFlavB[goodJets[i]] > btag_deepJet_Tight)
            HZZ2l2qNu_nTightBtagJets++;
        if (Jet_btagDeepFlavB[goodJets[i]] > btag_deepJet_Medium)
            HZZ2l2qNu_nMediumBtagJets++;
        if (Jet_btagDeepFlavB[goodJets[i]] > btag_deepJet_Loose)
            HZZ2l2qNu_nLooseBtagJets++;
    }

    return goodJets;
}

// Pre-selection for Fat-jets
std::vector<unsigned int> H4LTools::SelectedFatJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu)
{
    std::vector<unsigned int> goodJets;
    // unsigned nJ = (*nJet).Get()[0];
    // std::cout<<"DEBUG: FatJet_pt.size() = " << FatJet_pt.size() << std::endl;
    for (unsigned int i = 0; i < FatJet_pt.size(); i++)
    {
      int overlaptag = 0;
      TLorentzVector jettest;
      if (FatJet_pt[i] < 200)
        continue;
      if (fabs(FatJet_eta[i]) > 2.4)
        continue;

      jettest.SetPtEtaPhiM(FatJet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
      for (unsigned int ie = 0; ie < ele.size(); ie++)
      {
        TLorentzVector eletest;
        eletest.SetPtEtaPhiM(Electron_pt[ele[ie]], Electron_eta[ele[ie]], Electron_phi[ele[ie]], Electron_mass[ele[ie]]);
        if (eletest.DeltaR(jettest) < 0.8)
            overlaptag++;
      }
      for (unsigned int im = 0; im < mu.size(); im++)
      {
        TLorentzVector mutest;
        mutest.SetPtEtaPhiM(Muon_pt[mu[im]], Muon_eta[mu[im]], Muon_phi[mu[im]], Muon_mass[mu[im]]);
        if (mutest.DeltaR(jettest) < 0.8)
            overlaptag++;
      }

      if (overlaptag == 0)
      {
        // std::cout<<"jetindex: "<<i<<"jetID "<<(*Jet_jetId)[i]<<" puID "<<(*Jet_puId)[i]<<std::endl;
        if ((Jet_jetId[i] > 0) && (Jet_puId[i] == 7))
        {
            goodJets.push_back(i);
        }
      }
    }
    // std::cout<<"DEBUG: goodJets.size() = " << goodJets.size() << std::endl;

    return goodJets;
}
//  END: Pre-selection for Fat-jets

unsigned H4LTools::doFsrRecovery(TLorentzVector Lep){
    // This Function returns the index for the possible FsrPhoton
    unsigned int FsrIdx = 999; //only Idx>0 works, pay attention!
    std::vector<unsigned int> BestFsrPhotons;
    BestFsrPhotons = goodFsrPhotons();
    float dRl,dRlOverPt;
    dRl = 999;
    dRlOverPt = 999;
    for(unsigned int i=0;i<BestFsrPhotons.size();i++){
        TLorentzVector fsrcand;
        fsrcand.SetPtEtaPhiM(FsrPhoton_pt[BestFsrPhotons[i]],FsrPhoton_eta[BestFsrPhotons[i]],FsrPhoton_phi[BestFsrPhotons[i]],0);
        float dRlC,dRlOverPtC;
        dRlC = fsrcand.DeltaR(Lep);
        if ((dRlC<fsrphotondRlcut)&&(FsrPhoton_dROverEt2[BestFsrPhotons[i]]<fsrphotondRlOverPtcut)){
            if(FsrPhoton_dROverEt2[BestFsrPhotons[i]]<dRlOverPt){
                dRl = dRlC;
                dRlOverPt = FsrPhoton_dROverEt2[BestFsrPhotons[i]];
                FsrIdx = BestFsrPhotons[i];
            }
        }
    }

    return FsrIdx;

}

unsigned H4LTools::doFsrRecovery_Run3(std::vector<unsigned int> goodfsridx, unsigned lepidx, int lepflavor){//lepflavor 11 or 13

    unsigned matchedfsridx = 999;
    if(lepflavor == 11){
        for(unsigned fsridx=0; fsridx<goodfsridx.size(); fsridx++){
            if(FsrPhoton_electronIdx[goodfsridx[fsridx]] == lepidx){
                matchedfsridx = fsridx;
                break;
            }
        }
    }
    if(lepflavor == 13){
        for(unsigned fsridx=0; fsridx<goodfsridx.size(); fsridx++){
            if(FsrPhoton_muonIdx[goodfsridx[fsridx]] == lepidx){
                matchedfsridx = fsridx;
                break;
            }
        }
    }
    return matchedfsridx;
}
void H4LTools::BatchFsrRecovery_Run3(){
    unsigned fsridx;
    std::vector<unsigned> fsrlist;
    fsrlist = goodFsrPhotons();
    for(unsigned int i=0; i<Electron_pt.size(); i++){
        TLorentzVector fsr,lep;
        lep.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
        fsridx = doFsrRecovery_Run3(fsrlist,i,11);
        if(fsridx<900){
            fsr.SetPtEtaPhiM(FsrPhoton_pt[fsrlist[fsridx]], FsrPhoton_eta[fsrlist[fsridx]], FsrPhoton_phi[fsrlist[fsridx]], 0);
            lep = lep + fsr;
            Electrondressed_Run3.push_back(lep);
        }
        else{
            Electrondressed_Run3.push_back(lep);
        }
    }
    for(unsigned int j=0; j<Muon_pt.size(); j++){
        TLorentzVector fsr,lep;
        lep.SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
        fsridx = doFsrRecovery_Run3(fsrlist,j,13);
        if (DEBUG) std::cout << "fsridx = " << fsridx << std::endl;
        if(fsridx<900){
            fsr.SetPtEtaPhiM(FsrPhoton_pt[fsrlist[fsridx]], FsrPhoton_eta[fsrlist[fsridx]], FsrPhoton_phi[fsrlist[fsridx]], 0);
            lep = lep + fsr;
            Muondressed_Run3.push_back(lep);
        }
        else{
            Muondressed_Run3.push_back(lep);
        }
    }
}

std::vector<TLorentzVector> H4LTools::BatchFsrRecovery(std::vector<TLorentzVector> LepList){

    std::vector<TLorentzVector> LepFsrList;

    for(unsigned int i=0;i<LepList.size();i++){
        int Fsrtag;
        Fsrtag = doFsrRecovery(LepList[i]);
        if (Fsrtag>900){
            LepFsrList.push_back(LepList[i]);
            continue;
        }
        TLorentzVector FsrPhoton;
        FsrPhoton.SetPtEtaPhiM(FsrPhoton_pt[Fsrtag],FsrPhoton_eta[Fsrtag],FsrPhoton_phi[Fsrtag],0);
        TLorentzVector LepFsrRecovery;
        LepFsrRecovery = FsrPhoton + LepList[i];
        LepFsrList.push_back(LepFsrRecovery);

        if (DEBUG)
            std::cout << "Inside BatchFsrRecovery:: LepFsrRecovery.Pt() = " << LepFsrRecovery.Pt() << "\t before fsr recovery: " << LepList[i].Pt() << std::endl;
    }
    return LepFsrList;
}
std::vector<TLorentzVector> H4LTools::ElectronFsr(){
    std::vector<TLorentzVector> leplist;
    std::vector<TLorentzVector> leplistfsr;
    //unsigned nlep = (*nElectron).Get()[0];
    for(unsigned int i=0;i<Electron_pt.size();i++){
        TLorentzVector Lep;
        Lep.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
        leplist.push_back(Lep);
    }
    leplistfsr = BatchFsrRecovery(leplist);
    return leplistfsr;
}

std::vector<TLorentzVector> H4LTools::MuonFsr(){
    std::vector<TLorentzVector> leplist;
    std::vector<TLorentzVector> leplistfsr;
    //unsigned nlep = (*nMuon).Get()[0];
    for(unsigned int i=0;i<Muon_pt.size();i++){
        TLorentzVector Lep;
        Lep.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],Muon_mass[i]);
        leplist.push_back(Lep);
    }
    leplistfsr = BatchFsrRecovery(leplist);
    return leplistfsr;
}

std::vector<float> H4LTools::ElectronFsrPt(){
    std::vector<float> lepPt;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepPt.push_back(Electrondressed_Run3[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::ElectronFsrEta(){
    std::vector<float> lepEta;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepEta.push_back(Electrondressed_Run3[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::ElectronFsrPhi(){
    std::vector<float> lepPhi;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepPhi.push_back(Electrondressed_Run3[i].Phi());
    }
    return lepPhi;
}

std::vector<float> H4LTools::MuonFsrPt(){
    std::vector<float> lepPt;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepPt.push_back(Muondressed_Run3[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::MuonFsrEta(){
    std::vector<float> lepEta;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepEta.push_back(Muondressed_Run3[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::MuonFsrPhi(){
    std::vector<float> lepPhi;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepPhi.push_back(Muondressed_Run3[i].Phi());
    }
    return lepPhi;
}
/*std::vector<float> H4LTools::ElectronFsrPt(){
    std::vector<float> lepPt;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPt.push_back(leplistfsr[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::ElectronFsrEta(){
    std::vector<float> lepEta;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepEta.push_back(leplistfsr[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::ElectronFsrPhi(){
    std::vector<float> lepPhi;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPhi.push_back(leplistfsr[i].Phi());
    }
    return lepPhi;
}

std::vector<float> H4LTools::MuonFsrPt(){
    std::vector<float> lepPt;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPt.push_back(leplistfsr[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::MuonFsrEta(){
    std::vector<float> lepEta;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepEta.push_back(leplistfsr[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::MuonFsrPhi(){
    std::vector<float> lepPhi;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPhi.push_back(leplistfsr[i].Phi());
    }
    return lepPhi;
}*/

void H4LTools::LeptonSelection(){
    looseEle = goodLooseElectrons2012();
    looseMu = goodLooseMuons2012();
    bestEle = goodElectrons2015_noIso_noBdt(looseEle);
    bestMu = goodMuons2015_noIso_noPf(looseMu);
    Electronindex = bestEle;
    Muonindex = bestMu;
    if (year == 2022) AllEid = passTight_BDT_Id_ZZ4l();
    else AllEid = passTight_BDT_Id();
    AllMuid = passTight_Id();
    for (unsigned int iuj=0;iuj<looseEle.size();iuj++){
        if(AllEid[looseEle[iuj]]) tighteleforjetidx.push_back(looseEle[iuj]);
    }
    for (unsigned int juj=0;juj<looseMu.size();juj++){
        if(AllMuid[looseMu[juj]]) tightmuforjetidx.push_back(looseMu[juj]);
    }
    jetidx = SelectedJets(tighteleforjetidx,tightmuforjetidx);
    FatJetidx = SelectedFatJets(tighteleforjetidx, tightmuforjetidx);

    for(unsigned int ie=0; ie<Electronindex.size();ie++)
    {
        if (Electronindex[ie] >= Electron_pt.size())
        {
            std::cerr << "ERROR: Electronindex out of bounds! ie = " << ie
                      << ", Electronindex[ie] = " << Electronindex[ie]
                      << ", Electron_pt.size() = " << Electron_pt.size() << std::endl;
            continue;
        }
        if(Electron_pdgId[Electronindex[ie]]>0){
            Elechg.push_back(-1);
        }
        else{
            Elechg.push_back(1);
        }
        TLorentzVector Ele;
        if (DEBUG)
            std::cout << "L#524: Inside LeptonSelection:: Electron_pt[" << Electronindex[ie] << "] = " << Electron_pt[Electronindex[ie]] << std::endl;
        Ele.SetPtEtaPhiM(Electron_pt[Electronindex[ie]],Electron_eta[Electronindex[ie]],Electron_phi[Electronindex[ie]],Electron_mass[Electronindex[ie]]);
        Elelist.push_back(Ele);
        if (year == 2022)
        {
            ElelistFsr.push_back(Electrondressed_Run3[Electronindex[ie]]);
        }
        Eiso.push_back(Electron_pfRelIso03_all[Electronindex[ie]]);
        Eid.push_back(AllEid[Electronindex[ie]]);
    }

    for(unsigned int imu=0; imu<Muonindex.size();imu++){
        if(Muon_pdgId[Muonindex[imu]]>0){
            Muchg.push_back(-1);
        }
        else{
            Muchg.push_back(1);
        }
        TLorentzVector Mu(0.,0.,0.,0.);
        if (DEBUG)
        {
            std::cout << "L#545: Inside LeptonSelection:: Muon_pt[" << Muonindex[imu] << "] = " << Muon_pt[Muonindex[imu]] << std::endl;
            std::cout << "L#546: Inside LeptonSelection:: Muon_eta[" << Muonindex[imu] << "] = " << Muon_eta[Muonindex[imu]] << std::endl;
            std::cout << "L#547: Inside LeptonSelection:: Muon_phi[" << Muonindex[imu] << "] = " << Muon_phi[Muonindex[imu]] << std::endl;
            std::cout << "L#548: Inside LeptonSelection:: Muon_mass[" << Muonindex[imu] << "] = " << Muon_mass[Muonindex[imu]] << std::endl;
        }
        Mu.SetPtEtaPhiM(Muon_pt[Muonindex[imu]],Muon_eta[Muonindex[imu]],Muon_phi[Muonindex[imu]],Muon_mass[Muonindex[imu]]);

        if (DEBUG) std::cout << "L#552: Inside LeptonSelection:: Mulist size = " << Mulist.size() << std::endl;

        Mulist.push_back(Mu);

        if (DEBUG)
        {
            std::cout << "L#558: Inside LeptonSelection:: Mulist size = " << Mulist.size() << std::endl;
            std::cout << "L#559: Inside LeptonSelection:: Muondressed_Run3 size = " << Muondressed_Run3.size() << std::endl;
            std::cout << "L#560: Inside LeptonSelection:: Muonindex[" << imu << "] = " << Muonindex[imu] << std::endl;
            std::cout << "L#561: Muondressed_Run3[" << Muonindex[imu] << "].Pt() = " << Muondressed_Run3[Muonindex[imu]].Pt() << std::endl;
            std::cout << "L#561: Muondressed_Run3[" << Muonindex[imu] << "].Pt() = " << Muondressed_Run3[Muonindex[imu]].Pt() << std::endl;
            std::cout << "L#563: Muon_pfRelIso03_all[" << Muonindex[imu] << "] = " << Muon_pfRelIso03_all[Muonindex[imu]] << std::endl;
            std::cout << "L#564: Muon_pfRelIso03_all[" << Muonindex[imu] << "] = " << Muon_pfRelIso03_all[Muonindex[imu]] << std::endl;
        }
        if (year == 2022)
        {
            MulistFsr.push_back(Muondressed_Run3[Muonindex[imu]]);
        }
        muid.push_back(AllMuid[Muonindex[imu]]);
        Muiso.push_back(Muon_pfRelIso03_all[Muonindex[imu]]);
    }

    if (DEBUG)
    {
        std::cout << "Electronindex.size() = " << Electronindex.size() << std::endl;
        std::cout << "Muonindex.size() = " << Muonindex.size() << std::endl;
    }

    if (year != 2022)
    {
        ElelistFsr = BatchFsrRecovery(Elelist);
        MulistFsr = BatchFsrRecovery(Mulist);
    }

    if (DEBUG)
    {
        std::cout << "ElelistFsr.size() = " << ElelistFsr.size() << std::endl;
        std::cout << "MulistFsr.size() = " << MulistFsr.size() << std::endl;
        // print Eid, Muid, Eiso, Muiso
        for (unsigned int i = 0; i < Eid.size(); i++)
        {
            std::cout << "Eid[" << i << "] = " << Eid[i] << std::endl;
        }
    }

    if (DEBUG)
        std::cout << "Line#574: Size of Eid = " << Eid.size() << std::endl;

    for (unsigned int ae = 0; ae < Eid.size(); ae++)
    {
        float RelEleIsoNoFsr;
        RelEleIsoNoFsr = Eiso[ae];
        unsigned FsrEleidx;
        if (year == 2022)
            FsrEleidx = doFsrRecovery_Run3(goodFsrPhotons(), Electronindex[ae], 11);
        else
            FsrEleidx = doFsrRecovery(Elelist[ae]);

        if (DEBUG)
            std::cout << "FsrEleidx = " << FsrEleidx << std::endl;
        if (isFSR && (FsrEleidx < 900))
        {
            TLorentzVector fsrele;
            fsrele.SetPtEtaPhiM(FsrPhoton_pt[FsrEleidx], FsrPhoton_eta[FsrEleidx], FsrPhoton_phi[FsrEleidx], 0);
            if (Elelist[ae].DeltaR(fsrele) > 0.01)
            {
                RelEleIsoNoFsr = RelEleIsoNoFsr - FsrPhoton_pt[FsrEleidx] / Elelist[ae].Pt();
            }
        }

        if (DEBUG)
            std::cout << "Line#597: Eid[" << ae << "] = " << Eid[ae] << "  RelEleIsoNoFsr = " << RelEleIsoNoFsr << std::endl;
        // check size of Eid
        if (DEBUG)
            std::cout << "Line#600: Size of Eid = " << Eid.size() << std::endl;

        if ((Eid[ae] == true) && (RelEleIsoNoFsr < 0.35))
        {
            if (DEBUG)
                std::cout << "Line#605: Eid[" << ae << "] = " << Eid[ae] << "  RelEleIsoNoFsr = " << RelEleIsoNoFsr << "  nTightEle : " << nTightEle << std::endl;
            nTightEle++;
            TightEleindex.push_back(ae);
            nTightEleChgSum += Elechg[ae];
            TightElelep_index.push_back(Lepointer);
            Lepointer++;
            if (DEBUG)
                std::cout << "Line#612: isMC = " << isMC << std::endl;
            if (isMC && year == 2022) // FIXME: Generalise this for all years
            {
                if (DEBUG)
                    std::cout << "Length of Electron_genPartIdx = " << Electron_genPartIdx.size() << " ae = " << ae << std::endl;
                lep_genindex.push_back(Electron_genPartIdx[Electronindex[ae]]);
                if (DEBUG)
                    std::cout << "Line#616: lep_genindex[" << ae << "] = " << lep_genindex[ae] << std::endl;
            }
            else
            {
                lep_genindex.push_back(-1);
            }
        }
        if (DEBUG)
            std::cout << "nTightEle = " << nTightEle << std::endl;
    }

    if (DEBUG)
        std::cout << "Line#632: Size of Muid = " << muid.size() << std::endl;

    for(unsigned int amu=0; amu<muid.size();amu++)
    {
        if (DEBUG)
            std::cout << "amu = " << amu << " Muiso size: " << Muiso.size() << " year: " << year << std::endl;
        float RelIsoNoFsr;
        RelIsoNoFsr = Muiso[amu];
        unsigned int FsrMuonidx;
        if (year == 2022)
            FsrMuonidx = doFsrRecovery_Run3(goodFsrPhotons(), Muonindex[amu], 13);
	else
            FsrMuonidx = doFsrRecovery(Mulist[amu]);

        if (DEBUG) std::cout << "FsrMuonidx = " << FsrMuonidx << std::endl;
        if (isFSR && (FsrMuonidx < 900)){
            TLorentzVector fsrmuon;
            fsrmuon.SetPtEtaPhiM(FsrPhoton_pt[FsrMuonidx],FsrPhoton_eta[FsrMuonidx],FsrPhoton_phi[FsrMuonidx],0);
            if(Mulist[amu].DeltaR(fsrmuon)>0.01){
                RelIsoNoFsr = RelIsoNoFsr - FsrPhoton_pt[FsrMuonidx]/Mulist[amu].Pt();
            }
        }
        if (DEBUG)
            std::cout << "RelIsoNoFsr = " << RelIsoNoFsr << " " << "muid[amu] " << muid[amu] << " size of Muchg = " << Muchg.size() <<  std::endl;
        if((muid[amu]==true)&&(RelIsoNoFsr<0.35)){
            nTightMu++;
            TightMuindex.push_back(amu);
            nTightMuChgSum += Muchg[amu];
            TightMulep_index.push_back(Lepointer);
            Lepointer++;
            if (DEBUG)
                std::cout << "Line#662: isMC = " << isMC <<  " nTightMu : " << nTightMu << std::endl;
            if (isMC && year == 2022) // FIXME: Generalise this for all years
                lep_genindex.push_back(Muon_genPartIdx[Muonindex[amu]]); // FIXME: Check if this is correct
            else
                lep_genindex.push_back(-1);
            if (DEBUG)
                std::cout << "Line#666: lep_genindex[" << amu << "] = " << lep_genindex[amu] << std::endl;
        }
    }
    if (DEBUG)
        std::cout << "nTightMu = " << nTightMu << std::endl;
}

/**
 * This function identifies and selects Z boson candidates from the list of tight leptons.
 * It checks for pairs of leptons that form a Z boson candidate based on their invariant mass
 * and opposite charge selection criteria. The function updates the Z boson candidate list
 * and returns true if at least one Z boson candidate is found, otherwise returns false.
 */
bool H4LTools::findZCandidates(){
    TLorentzVector z1,z2;

    if (nTightEle>=4) {
        cut4e++;
        flag4e = true;
    }
    else if (nTightMu>=4){
        cut4mu++;
        flag4mu = true;
    }
    else if ((nTightMu>=2)&&(nTightEle>=2)){
        cut2e2mu++;
        flag2e2mu = true;
    }

    if (DEBUG)
    {
        std::cout << "nTightEle: " << nTightEle << ", nTightMu: " << nTightMu << std::endl;

        for (unsigned int ke = 0; ke < (TightEleindex.size()); ke++)
        {
            std::cout << "ElelistFsr[TightEleindex[ke]].Pt() = " << ElelistFsr[TightEleindex[ke]].Pt() << std::endl;
        }

        for (unsigned int kmu = 0; kmu < (TightMuindex.size()); kmu++)
        {
            std::cout << "MulistFsr[TightMuindex[kmu]].Pt() = " << MulistFsr[TightMuindex[kmu]].Pt() << std::endl;
        }
    }

    if(TightEleindex.size()>1)
    {
        for(unsigned int ke=0; ke<(TightEleindex.size()-1);ke++)
        {
            for(unsigned int je=ke+1;je<TightEleindex.size();je++)
            {
                // opposite charge requirement
                if ((Elechg[TightEleindex[ke]]+Elechg[TightEleindex[je]]) != 0)
                    continue;
                TLorentzVector Zcan;
                Zcan = ElelistFsr[TightEleindex[ke]] + ElelistFsr[TightEleindex[je]];

                // invariant mass requirement: MZcutdown < MZ < MZcutup
                if (!(Zcan.M() > MZcutdown && Zcan.M() < MZcutup))
                    continue;

                Zlist.push_back(Zcan);
                Zlep1index.push_back(TightEleindex[ke]);
                Zlep2index.push_back(TightEleindex[je]);
                Zlep1lepindex.push_back(TightElelep_index[ke]);
                Zlep2lepindex.push_back(TightElelep_index[je]);
                Zflavor.push_back(11);
                Zlep1pt.push_back(ElelistFsr[TightEleindex[ke]].Pt());
                Zlep2pt.push_back(ElelistFsr[TightEleindex[je]].Pt());
                Zlep1eta.push_back(ElelistFsr[TightEleindex[ke]].Eta());
                Zlep2eta.push_back(ElelistFsr[TightEleindex[je]].Eta());
                Zlep1phi.push_back(ElelistFsr[TightEleindex[ke]].Phi());
                Zlep2phi.push_back(ElelistFsr[TightEleindex[je]].Phi());
                Zlep1mass.push_back(ElelistFsr[TightEleindex[ke]].M());
                Zlep2mass.push_back(ElelistFsr[TightEleindex[je]].M());
                Zlep1ptNoFsr.push_back(Elelist[TightEleindex[ke]].Pt());
                Zlep2ptNoFsr.push_back(Elelist[TightEleindex[je]].Pt());
                Zlep1etaNoFsr.push_back(Elelist[TightEleindex[ke]].Eta());
                Zlep2etaNoFsr.push_back(Elelist[TightEleindex[je]].Eta());
                Zlep1phiNoFsr.push_back(Elelist[TightEleindex[ke]].Phi());
                Zlep2phiNoFsr.push_back(Elelist[TightEleindex[je]].Phi());
                Zlep1massNoFsr.push_back(Elelist[TightEleindex[ke]].M());
                Zlep2massNoFsr.push_back(Elelist[TightEleindex[je]].M());
                Zlep1chg.push_back(Elechg[TightEleindex[ke]]);
                Zlep2chg.push_back(Elechg[TightEleindex[je]]);
            }
        }
    }

    if(TightMuindex.size()>1)
    {
        for(unsigned int kmu=0; kmu<(TightMuindex.size()-1);kmu++)
        {
            for(unsigned int jmu=kmu+1;jmu<TightMuindex.size();jmu++)
            {
                // opposite charge requirement
                if ((Muchg[TightMuindex[kmu]]+Muchg[TightMuindex[jmu]]) != 0)
                    continue;
                TLorentzVector Zcan;
                Zcan = MulistFsr[TightMuindex[kmu]] + MulistFsr[TightMuindex[jmu]];

                // invariant mass requirement : MZcutdown < MZ < MZcutup
                if (!(Zcan.M() > MZcutdown && Zcan.M() < MZcutup))
                    continue;
                Zlist.push_back(Zcan);
                Zlep1index.push_back(TightMuindex[kmu]);
                Zlep2index.push_back(TightMuindex[jmu]);
                Zlep1lepindex.push_back(TightMulep_index[kmu]);
                Zlep2lepindex.push_back(TightMulep_index[jmu]);
                Zflavor.push_back(13);
                Zlep1pt.push_back(MulistFsr[TightMuindex[kmu]].Pt());
                Zlep2pt.push_back(MulistFsr[TightMuindex[jmu]].Pt());
                Zlep1eta.push_back(MulistFsr[TightMuindex[kmu]].Eta());
                Zlep2eta.push_back(MulistFsr[TightMuindex[jmu]].Eta());
                Zlep1phi.push_back(MulistFsr[TightMuindex[kmu]].Phi());
                Zlep2phi.push_back(MulistFsr[TightMuindex[jmu]].Phi());
                Zlep1mass.push_back(MulistFsr[TightMuindex[kmu]].M());
                Zlep2mass.push_back(MulistFsr[TightMuindex[jmu]].M());
                Zlep1ptNoFsr.push_back(Mulist[TightMuindex[kmu]].Pt());
                Zlep2ptNoFsr.push_back(Mulist[TightMuindex[jmu]].Pt());
                Zlep1etaNoFsr.push_back(Mulist[TightMuindex[kmu]].Eta());
                Zlep2etaNoFsr.push_back(Mulist[TightMuindex[jmu]].Eta());
                Zlep1phiNoFsr.push_back(Mulist[TightMuindex[kmu]].Phi());
                Zlep2phiNoFsr.push_back(Mulist[TightMuindex[jmu]].Phi());
                Zlep1massNoFsr.push_back(Mulist[TightMuindex[kmu]].M());
                Zlep2massNoFsr.push_back(Mulist[TightMuindex[jmu]].M());
                Zlep1chg.push_back(Muchg[TightMuindex[kmu]]);
                Zlep2chg.push_back(Muchg[TightMuindex[jmu]]);
            }
        }
    }

    for (unsigned int znofsr = 0; znofsr<Zlist.size(); znofsr++){
        TLorentzVector Zlep1nofsr,Zlep2nofsr,Zcannofsr;
        Zlep1nofsr.SetPtEtaPhiM(Zlep1ptNoFsr[znofsr],Zlep1etaNoFsr[znofsr],Zlep1phiNoFsr[znofsr],Zlep1massNoFsr[znofsr]);
        Zlep2nofsr.SetPtEtaPhiM(Zlep2ptNoFsr[znofsr],Zlep2etaNoFsr[znofsr],Zlep2phiNoFsr[znofsr],Zlep2massNoFsr[znofsr]);
        Zcannofsr = Zlep1nofsr + Zlep2nofsr;
        Zlistnofsr.push_back(Zcannofsr);
    }
    Zsize = Zlist.size();

    if (DEBUG)
    {
        std::cout << "Zsize: " << Zsize << std::endl;
        std::cout << "Zlep1pt size: " << Zlep1pt.size() << std::endl;
        std::cout << "Zlep2pt size: " << Zlep2pt.size() << std::endl;
        // if size is greater than 0, cout the Zlep1pt and Zlep2pt
        if (Zlep1pt.size() > 0)
            std::cout << "Zlep1pt: " << Zlep1pt[0] << std::endl;
        if (Zlep2pt.size() > 0)
            std::cout << "Zlep2pt: " << Zlep2pt[0] << std::endl;
    }

    return Zsize > 0;
}


bool H4LTools::ZZSelection_4l(){

    bool foundZZCandidate = false;
    //std::cout << " Inside the 4l loop in .cc file" << std::endl;
    if(!findZCandidates()){
        return foundZZCandidate;
    }
    if((nTightMu+nTightEle)<4){
        return foundZZCandidate;
    }

    if((abs(nTightEleChgSum)+abs(nTightMuChgSum))>(nTightMu+nTightEle-4)){
        return foundZZCandidate;
    }
    if(Zsize<2){
        return foundZZCandidate;
    }
    //Find ZZ candidate
    std::vector<int> Z1CanIndex;
    std::vector<int> Z2CanIndex;
    int ghosttag = 0, QCDtag=0, lepPtTag = 0;
    for (unsigned int m=0; m<(Zlist.size()-1); m++){
        for (unsigned int n=m+1; n<Zlist.size(); n++){
            if (Zflavor[m]==Zflavor[n]){
               if ((Zlep1index[m] == Zlep1index[n])||(Zlep2index[m] == Zlep1index[n])) continue;  //non-overlapping
               if ((Zlep1index[m] == Zlep2index[n])||(Zlep2index[m] == Zlep2index[n])) continue;
            }
            if (Zlist[m].DeltaR(Zlist[n])<0.02) continue; //ghost removal
            ghosttag++;
            bool nPassPt20;
            int nPassPt10;
            nPassPt20 = (Zlep1pt[m]>20) || (Zlep2pt[m]>20) || (Zlep1pt[n]>20) || (Zlep2pt[n]>20);
            nPassPt10 = 0;
            if (Zlep1pt[m]>10) nPassPt10 += 1;
            if (Zlep2pt[m]>10) nPassPt10 += 1;
            if (Zlep1pt[n]>10) nPassPt10 += 1;
            if (Zlep2pt[n]>10) nPassPt10 += 1;
            if (nPassPt10 < 2) continue;
            if (nPassPt20 == false) continue; //lep Pt requirements
            lepPtTag++;
            if ((Zlep1chg[m]+Zlep1chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep1ptNoFsr[m],Zlep1etaNoFsr[m],Zlep1phiNoFsr[m],Zlep1massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep2ptNoFsr[n],Zlep2etaNoFsr[n],Zlep2phiNoFsr[n],Zlep2massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;  //QCD Supressionas
            }
            if ((Zlep1chg[m]+Zlep2chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep1ptNoFsr[m],Zlep1etaNoFsr[m],Zlep1phiNoFsr[m],Zlep1massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep2ptNoFsr[n],Zlep2etaNoFsr[n],Zlep2phiNoFsr[n],Zlep2massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }
            if ((Zlep2chg[m]+Zlep1chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep2ptNoFsr[m],Zlep2etaNoFsr[m],Zlep2phiNoFsr[m],Zlep2massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep1ptNoFsr[n],Zlep1etaNoFsr[n],Zlep1phiNoFsr[n],Zlep1massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }
            if ((Zlep2chg[m]+Zlep2chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep2ptNoFsr[m],Zlep2etaNoFsr[m],Zlep2phiNoFsr[m],Zlep2massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep1ptNoFsr[n],Zlep1etaNoFsr[n],Zlep1phiNoFsr[n],Zlep1massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }
            QCDtag++;
            if ((Zlist[m].M()<40) && (Zlist[n].M()<40))  continue; //Z1 mass

            TLorentzVector zZ1,zZ2;
            if (fabs(Zlist[m].M()-Zmass)<fabs(Zlist[n].M()-Zmass)){
                zZ1 = Zlist[m];
                zZ2 = Zlist[n];
            }
            else{
                zZ1 = Zlist[n];
                zZ2 = Zlist[m];
            }

            bool passSmartCut = true;
            if (Zflavor[m]==Zflavor[n]){
                TLorentzVector Za,Zb,lepM1,lepM2,lepN1,lepN2;
                int lepM1chg,lepM2chg,lepN1chg,lepN2chg;
                lepM1.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m]);
                lepM2.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m]);
                lepN1.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n]);
                lepN2.SetPtEtaPhiM(Zlep2pt[n],Zlep2eta[n],Zlep2phi[n],Zlep2mass[n]);
                lepM1chg = Zlep1chg[m];
                lepM2chg = Zlep2chg[m];
                lepN1chg = Zlep1chg[n];
                lepN2chg = Zlep2chg[n];
                if(lepM1chg == lepN1chg){
                    Za = lepM1 + lepN2;
                    Zb = lepN1 + lepM2;
                }
                else{
                    Za = lepM1 + lepN1;
                    Zb = lepN2 + lepM2;
                }
                if (fabs(Za.M()-Zmass)<fabs(Zb.M()-Zmass)){
                    if ( (fabs(Za.M()-Zmass)<abs(zZ1.M()-Zmass)) && (Zb.M()<12) ) passSmartCut=false;
                }

                else{
                    if ( (fabs(Zb.M()-Zmass)<fabs(zZ1.M()-Zmass)) && (Za.M()<12) ) passSmartCut=false;
                }
            }
            if (passSmartCut==false) continue ;
            if (zZ1.M()+zZ2.M()<MZZcut) continue;
            foundZZCandidate = true;
            if(Zlist[m].M()>Zlist[n].M()){
                Z1CanIndex.push_back(m);
                Z2CanIndex.push_back(n);
            }
            else{
                Z1CanIndex.push_back(n);
                Z2CanIndex.push_back(m);
            }


        }
    }
    if(ghosttag){
        if (flag2e2mu) cutghost2e2mu++;
        if (flag4e) cutghost4e++;
        if (flag4mu) cutghost4mu++;
    }
    if(lepPtTag){
        if (flag2e2mu) cutLepPt2e2mu++;
        if (flag4e) cutLepPt4e++;
        if (flag4mu) cutLepPt4mu++;
    }
    if(QCDtag){
        if (flag2e2mu) cutQCD2e2mu++;
        if (flag4e) cutQCD4e++;
        if (flag4mu) cutQCD4mu++;
    }
    if(foundZZCandidate == false){
        return foundZZCandidate;
    }
    if (flag2e2mu) cutZZ2e2mu++;
    if (flag4e) cutZZ4e++;
    if (flag4mu) cutZZ4mu++;
    int Z1index,Z2index;
    Z1index = Z1CanIndex[0];
    Z2index = Z2CanIndex[0];
    float Z2Ptsum;
    Z2Ptsum = Zlep1pt[Z2index] + Zlep2pt[Z2index];
    if(Z1CanIndex.size()>1){
        for(unsigned int iz=0;iz<Z1CanIndex.size();iz++){
            if (Z1index==Z1CanIndex[iz]){
                if((Zlep1pt[Z2CanIndex[iz]] + Zlep2pt[Z2CanIndex[iz]])>Z2Ptsum){
                    Z1index = Z1CanIndex[iz];
                    Z2index = Z2CanIndex[iz];
                    Z2Ptsum = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                }
            }
            if(fabs(Zlist[Z1CanIndex[iz]].M()-Zmass)<fabs(Zlist[Z1index].M()-Zmass)){
                Z1index = Z1CanIndex[iz];
                Z2index = Z2CanIndex[iz];
                Z2Ptsum = Zlep1pt[Z2index] + Zlep2pt[Z2index];
            }
        }
    }


    Z1 = Zlist[Z1index];
    Z2 = Zlist[Z2index];

    Z1nofsr = Zlistnofsr[Z1index];
    Z2nofsr = Zlistnofsr[Z2index];
    ZZsystem = Z1+Z2;
    ZZsystemnofsr = Z1nofsr+Z2nofsr;

    /*if(abs(ZZsystemnofsr.M()-ZZsystem.M())>0.000001){
        std::cout<<"FSR works "<<abs(ZZsystemnofsr.M()-ZZsystem.M())<<std::endl;
        std::cout<<"FSR: "<<ZZsystem.M()<<" noFSR:"<<ZZsystemnofsr.M()<<std::endl;
    }*/

    float massZZ;
    if (isFSR) massZZ = ZZsystem.M();
    else massZZ = ZZsystemnofsr.M();
    if ((massZZ>HiggscutDown)&&(massZZ<HiggscutUp)){
        if (flag2e2mu) cutm4l2e2mu++;
        if (flag4e) cutm4l4e++;
        if (flag4mu) cutm4l4mu++;
    }

    unsigned int jet1index, jet2index;
    jet1index = 99;
    jet2index = 99;
    if(jetidx.size()>0)
    {
        if(jetidx.size()==1)
        {
            jet1index = jetidx[0];
        }
        if(jetidx.size()==2)
        {
            jet1index = jetidx[0];
            jet2index = jetidx[1];
            if(Jet_pt[jetidx[1]]>Jet_pt[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }
        }
        if(jetidx.size()>2)
        {
            jet1index = jetidx[0];
            jet2index = jetidx[1];
            if(Jet_pt[jetidx[1]]>Jet_pt[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }
            for (unsigned int pj=2;pj<jetidx.size();pj++){
                if((Jet_pt[jetidx[pj]]>jet1index)&&(Jet_pt[jetidx[pj]]>jet2index)){
                    jet1index = jetidx[pj];
                }
                if(Jet_pt[jetidx[pj]]>jet2index){
                    jet2index = jetidx[pj];
                }
            }
        }
    }
    TLorentzVector Jet1,Jet2;
    SimpleParticleCollection_t associated;
    if(jetidx.size()>0){
        Jet1.SetPtEtaPhiM(Jet_pt[jet1index],Jet_eta[jet1index],Jet_phi[jet1index],Jet_mass[jet1index]);
        associated.push_back(SimpleParticle_t(0, Jet1));
        pTj1 = Jet1.Pt();
        etaj1 = Jet1.Eta();
        phij1 = Jet1.Phi();
        mj1 = Jet1.M();
        if(jetidx.size()>1){
           Jet2.SetPtEtaPhiM(Jet_pt[jet2index],Jet_eta[jet2index],Jet_phi[jet2index],Jet_mass[jet2index]);
           associated.push_back(SimpleParticle_t(0, Jet2));
           pTj2 = Jet2.Pt();
           etaj2 = Jet2.Eta();
           phij2 = Jet2.Phi();
           mj2 = Jet2.M();
        }
    }



    SimpleParticleCollection_t daughters;
    TLorentzVector Lep1,Lep2,Lep3,Lep4;

    Lep1.SetPtEtaPhiM(Zlep1pt[Z1index],Zlep1eta[Z1index],Zlep1phi[Z1index],Zlep1mass[Z1index]);
    Lep2.SetPtEtaPhiM(Zlep2pt[Z1index],Zlep2eta[Z1index],Zlep2phi[Z1index],Zlep2mass[Z1index]);
    Lep3.SetPtEtaPhiM(Zlep1pt[Z2index],Zlep1eta[Z2index],Zlep1phi[Z2index],Zlep1mass[Z2index]);
    Lep4.SetPtEtaPhiM(Zlep2pt[Z2index],Zlep2eta[Z2index],Zlep2phi[Z2index],Zlep2mass[Z2index]);

    if ((Zflavor[Z1index]==11)&&(Zflavor[Z2index]==11)) RecoFourEEvent=true;
    if ((Zflavor[Z1index]==13)&&(Zflavor[Z2index]==13)) RecoFourMuEvent=true;
    if ((Zflavor[Z1index]==11)&&(Zflavor[Z2index]==13)) RecoTwoETwoMuEvent=true;
    if ((Zflavor[Z1index]==13)&&(Zflavor[Z2index]==11)) RecoTwoMuTwoEEvent=true;
    lep_Hindex[0] = Zlep1lepindex[Z1index];
    lep_Hindex[1] = Zlep2lepindex[Z1index];
    lep_Hindex[2] = Zlep1lepindex[Z2index];
    lep_Hindex[3] = Zlep2lepindex[Z2index];
    pTL1 = Lep1.Pt();
    etaL1 = Lep1.Eta();
    phiL1 = Lep1.Phi();
    massL1 = Lep1.M();
    pTL2 = Lep2.Pt();
    etaL2 = Lep2.Eta();
    phiL2 = Lep2.Phi();
    massL2 = Lep2.M();
    pTL3 = Lep3.Pt();
    etaL3 = Lep3.Eta();
    phiL3 = Lep3.Phi();
    massL3 = Lep3.M();
    pTL4 = Lep4.Pt();
    etaL4 = Lep4.Eta();
    phiL4 = Lep4.Phi();
    massL4 = Lep4.M();


    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep1chg[Z1index], Lep1));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep2chg[Z1index], Lep2));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep1chg[Z2index], Lep3));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep2chg[Z2index], Lep4));
    me_0plus_JHU=999.0; me_qqZZ_MCFM=999.0; p0plus_m4l=999.0; bkg_m4l=999.0; D_bkg_kin=999.0; D_bkg=999.0;
    D_bkg_kin_vtx_BS=999.0;

    p0minus_VAJHU=999.0; pg1g4_VAJHU=999.0; Dgg10_VAMCFM=999.0; D_g4=999.0; D_g1g4=999.0; D_0m=999.0; D_CP=999.0; D_0hp=999; D_int=999.0;D_L1=999.0; D_L1_int=999.0; D_L1Zg=999.0; D_L1Zgint=999.0;
    p0plus_VAJHU=9999.0; p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen=999.0; pDL1_VAJHU=999.0; pD_L1Zgint=999.0; p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen=999.0; p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen=999.0, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen=999.0, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen=999.0, p0plus_VAJHU=999.0;

    mela->setInputEvent(&daughters, &associated, 0, 0);
    mela->setCurrentCandidateFromIndex(0);
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(me_0plus_JHU, true);
    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_VAJHU, true);
    // additional probabilities   GG_SIG_ghg2_1_ghz2_1_JHUGen
    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0plus_VAJHU, true);

    // p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, Couplings:ghg2=1,0;ghz1=1,0;ghz2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz2_1_JHUGen
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_2][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, true);    //FIXME

    p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen -= p0plus_VAJHU+me_0plus_JHU;


    // p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen
    // Couplings:ghg2=1,0;ghz1_prime2=10000,0
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1_PRIME2][0]=10000.;
    // (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, true);    //FIXME

	// p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen/1e8, ghg2=1,0;ghz1=1,0;ghz1_prime2=10000,0, Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen
	mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1_PRIME2][0]=10000.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, true);    //FIXME
    p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen -= p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen+me_0plus_JHU;

	// p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, ghg2=1,0;ghzgs1_prime2=10000,0
	mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_ZA_1_PRIME2][0]=10000.;
    // (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, true);    //FIXME

	// p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, ghg2=1,0;ghzgs1_prime2=10000,0
	mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_ZA_1_PRIME2][0]=10000.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, true);    //FIXME

	p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen -= me_0plus_JHU+p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;

	// pg1g4_VAJHU=0.0;
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(pg1g4_VAJHU, true);

    pg1g4_VAJHU -= me_0plus_JHU+p0minus_VAJHU;

    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela->computeP(me_qqZZ_MCFM, true);

    mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);

    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela->computePM4l(TVar::SMSyst_None, p0plus_m4l);

    mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
    mela->computePM4l(TVar::SMSyst_None, bkg_m4l);
    D_0m = me_0plus_JHU / (me_0plus_JHU + (p0minus_VAJHU * pow(getDg4Constant(massZZ),2)));
	D_CP = pg1g4_VAJHU / (2 * sqrt(me_0plus_JHU * p0minus_VAJHU ));
    D_0hp = me_0plus_JHU / (me_0plus_JHU + (p0plus_VAJHU * pow(getDg2Constant(massZZ),2)));
	D_int = p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * sqrt(me_0plus_JHU * p0plus_VAJHU));
    D_L1 = me_0plus_JHU / (me_0plus_JHU + ((p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(getDL1Constant(massZZ),2)));
	D_L1Zg = me_0plus_JHU / (me_0plus_JHU + ((p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(getDL1ZgsConstant(massZZ),2)));
    mela->resetInputEvent();
    return foundZZCandidate;
}

bool H4LTools::GetExactlyTwoTightLeptons()
{
    if (DEBUG)
        std::cout << "Inside function GetExactlyTwoTightLeptons()" << std::endl;
    bool foundZ1Candidate = false;

    if (!findZCandidates())
    {
        return foundZ1Candidate;
    }
    if (!(nTightMu == 2 || nTightEle == 2))
    {
        return foundZ1Candidate;
    }
    if (DEBUG)
        std::cout << "Number of leptons: (Mu, Ele, Total): " << nTightMu << ", " << nTightEle << ", " << nTightMu + nTightEle << std::endl;
    // Set HZZ2l2qNu_isELE to true if there are 2 tight electrons, false if there are 2 tight muons
    HZZ2l2qNu_isELE = (nTightEle == 2);
    HZZ2l2qNu_cut2l++;

    if (DEBUG)
        std::cout << "nTightEleChgSum: " << nTightEleChgSum << "\tnTightMuChgSum: " << nTightMuChgSum << std::endl;

    // Check if the absolute values of nTightEleChgSum and nTightMuChgSum are not zero
    if (!(nTightEleChgSum == 0 && nTightMuChgSum == 0))
    {
        return foundZ1Candidate;
    }
    HZZ2l2qNu_cutOppositeCharge++;
    HZZ2l2qNu_cutOppositeChargeFlag = true;

    if (DEBUG)
        std::cout << "Zlist size: " << Zlist.size() << std::endl;

    if ((Zlep1pt[0] < HZZ2l2nu_Leading_Lep_pT || Zlep2pt[0] < HZZ2l2nu_SubLeading_Lep_pT))
    {
        return foundZ1Candidate;
    }
    HZZ2l2qNu_cutpTl1l2++;

    // eta < HZZ2l2nu_Lep_eta
    if (fabs(Zlep1eta[0]) > HZZ2l2nu_Lep_eta || fabs(Zlep2eta[0]) > HZZ2l2nu_Lep_eta)
    {
        return foundZ1Candidate;
    }
    HZZ2l2qNu_cutETAl1l2++;

    // Find ZZ candidate
    std::vector<int> Z1CanIndex;
    for (unsigned int m = 0; m < (Zlist.size()); m++)
    {
        Z1CanIndex.push_back(m);
    }

    int Z1index;
    Z1index = Z1CanIndex[0];
    Z1 = Zlist[Z1index];
    Z1nofsr = Zlistnofsr[Z1index];

    // The invariant mass of dilepton system within 15 GeV of the known Z boson mass, ensuring that the pair likely originates from a Z-boson decay
    if (fabs(Z1.M() - Zmass) > HZZ2l2nu_M_ll_Window)
    {
        return foundZ1Candidate;
    }
    HZZ2l2qNu_cutmZ1Window++;

    // Also, the transverse momentum of the dilepton system should be > 55 GeV, indicating boosted Z-boson which is coming from High-mass Higgs boson
    if (Z1.Pt() < HZZ2l2nu_Pt_ll)
    {
        return foundZ1Candidate;
    }
    HZZ2l2qNu_cutZ1Pt++;

    if (DEBUG)
        std::cout << "Found Z1 candidate: " << foundZ1Candidate << std::endl;

    TLorentzVector Lep1, Lep2;
    Lep1.SetPtEtaPhiM(Zlep1pt[Z1index], Zlep1eta[Z1index], Zlep1phi[Z1index], Zlep1mass[Z1index]);
    Lep2.SetPtEtaPhiM(Zlep2pt[Z1index], Zlep2eta[Z1index], Zlep2phi[Z1index], Zlep2mass[Z1index]);
    pTL1 = Lep1.Pt();
    etaL1 = Lep1.Eta();
    phiL1 = Lep1.Phi();
    massL1 = Lep1.M();
    pTL2 = Lep2.Pt();
    etaL2 = Lep2.Eta();
    phiL2 = Lep2.Phi();
    massL2 = Lep2.M();

    if(Lep1.DeltaR(Lep2)<0.3)
        return foundZ1Candidate;

    // If the event passes all the cuts, set the flag `foundZ1Candidate` to true
    foundZ1Candidate = true;

    return foundZ1Candidate;
}

// emu control region (Z1 candidate selection) //
bool H4LTools::GetZ1_emuCR()
{
    if (DEBUG)
        std::cout << "*****$$$$*****Inside function GetZ1_emuCR()" << std::endl;
    bool foundZ1_emuCRCandidate = false;
    if (!((nTightMu == 1) && (nTightEle == 1)))
    {
        return foundZ1_emuCRCandidate;
    }

    HZZemuCR_cut2l++;
    if (DEBUG)
        std::cout << "*****$$$$*****Number of leptons inside emu control region: (Mu, Ele, Total): " << nTightMu << ", " << nTightEle << ", " << nTightMu + nTightEle << std::endl;

    if ((TightEleindex.size() == 1) && (TightMuindex.size() == 1))
    {
        if (DEBUG){
            std::cout << "*****$$$$*****Size of tight ele: " << TightEleindex[0] << "\t " << ElelistFsr[TightEleindex[0]].Pt() << std::endl;
            std::cout << "*****$$$$*****Size of tight mu: " << TightMuindex[0] << "\t " << MulistFsr[TightMuindex[0]].Pt() << std::endl;
	}
        Z_emuCRlep1pt.push_back(ElelistFsr[TightEleindex[0]].Pt());
        Z_emuCRlep2pt.push_back(MulistFsr[TightMuindex[0]].Pt());
        Z_emuCRlep1eta.push_back(ElelistFsr[TightEleindex[0]].Eta());
        Z_emuCRlep2eta.push_back(MulistFsr[TightMuindex[0]].Eta());
        Z_emuCRlep1phi.push_back(ElelistFsr[TightEleindex[0]].Phi());
        Z_emuCRlep2phi.push_back(MulistFsr[TightMuindex[0]].Phi());
        Z_emuCRlep1mass.push_back(ElelistFsr[TightEleindex[0]].M());
        Z_emuCRlep2mass.push_back(MulistFsr[TightMuindex[0]].M());
    }

     if (DEBUG)
         std::cout << "##Zlep1pt,Zlep2pt (emu control region): " << ElelistFsr[TightEleindex[0]].Pt() << ", " << MulistFsr[TightMuindex[0]].Pt() << std::endl;

    Z1 = ElelistFsr[TightEleindex[0]] + MulistFsr[TightMuindex[0]];
    TLorentzVector Lep1, Lep2;
    Lep1 = ElelistFsr[TightEleindex[0]];
    Lep2 = MulistFsr[TightMuindex[0]];

    pTL1 = ElelistFsr[TightEleindex[0]].Pt();
    pTL2 = MulistFsr[TightMuindex[0]].Pt();
    etaL1 = ElelistFsr[TightEleindex[0]].Eta();
    etaL2 = MulistFsr[TightMuindex[0]].Eta();
    phiL1 = ElelistFsr[TightEleindex[0]].Phi();
    phiL2 = MulistFsr[TightMuindex[0]].Phi();
    massL1 = ElelistFsr[TightEleindex[0]].M();
    massL2 = MulistFsr[TightMuindex[0]].M();

    /// pT selection
    if ((pTL1 < HZZ2l2nu_Leading_Lep_pT || pTL2 < HZZ2l2nu_SubLeading_Lep_pT))
    {
        return foundZ1_emuCRCandidate;
    }
    HZZemuCR_cutpTl1l2++;

    if (DEBUG)
        std::cout << "*****$$$$*****Zlep1pt,Zlep2pt (emu control region): " << pTL1 << ", " << pTL2 << std::endl;

    /// eta selection
    if (fabs(etaL1) > HZZ2l2nu_Lep_eta || fabs(etaL2) > HZZ2l2nu_Lep_eta)
    {
        return foundZ1_emuCRCandidate;
    }
    HZZemuCR_cutETAl1l2++;
    if (DEBUG)
        std::cout << "*****$$$$*****Zlep1eta,Zlep2eta (emu control region): " << etaL1 << ", " << etaL2 << std::endl;
    // std::cout << "##HELLO##Z_emu mass: " << Z1_emuCR.M() <<  std::endl;
    // std::cout << "##HELLO#Z_emu Pt: " << Z1_emuCR.Pt() <<  std::endl;

    if (fabs(Z1.M() - Zmass) > 160)
    {
        return foundZ1_emuCRCandidate;                     //UNCOMMENT FOR THE EVENTS IN Z1 MASS WINDOW
    }
    HZZemuCR_cutmZ1Window++;

    if (DEBUG)
        std::cout << "*****$$$$***** Z_emu mass: " << Z1.M() << std::endl;

    /*
    //side band
    if (!(((Z1.M() > 40) && (Z1.M() < 70)) || ((Z1.M() > 110) && (Z1.M() < 200))))
    {
        return foundZ1_emuCRCandidate;                     //UNCOMMENT FOR THE EVENTS IN SIDE BAND REGION
    }
    HZZemuCR_cutmZ1Window_SB++;
    std::cout << "*****$$$$***** Z_emu mass_sideband: " << Z1.M() << std::endl;
    */
    /// pT selection of dilepton
    if (Z1.Pt() < 25)
    {
        return foundZ1_emuCRCandidate;
    }
    HZZemuCR_cutZ1Pt++;
    foundZ1_emuCRCandidate = true;
    if (DEBUG)
        std::cout << "*****$$$$***** Found Z1_emuCR candidate: " << foundZ1_emuCRCandidate << std::endl;

    return foundZ1_emuCRCandidate;
}

bool H4LTools::ZZSelection_2l2q()
{
    if (DEBUG)
    {
        std::cout << "==> Inside ZZSelection_2l2q" << std::endl;
    }
    bool foundZZCandidate = false;

    jetidx = SelectedJets(tighteleforjetidx, tightmuforjetidx);
    FatJetidx = SelectedFatJets(tighteleforjetidx, tightmuforjetidx);

    if (FatJetidx.size() > 0 || jetidx.size()>=2)
    {
        if (FatJetidx.size() > 0)
        {
            for (unsigned int i = 0; i < FatJetidx.size(); i++)
            {
                if (FatJet_PNZvsQCD[FatJetidx[i]] < 0.9) continue;
                if (FatJet_pt[FatJetidx[i]] < 200.0) continue;
                // if (FatJet_msoftdrop[FatJetidx[i]] < 40.0) continue;
                // if (FatJet_msoftdrop[FatJetidx[i]] > 180.0) continue;

                foundZZCandidate = true;
                isBoosted2l2q = true;
                cut2l1J++;
                cut2l1Jor2j++;

                boostedJet_PNScore = FatJet_PNZvsQCD[FatJetidx[i]];
                boostedJet_Index = FatJetidx[i];

                Z2.SetPtEtaPhiM(FatJet_pt[FatJetidx[i]], FatJet_eta[FatJetidx[i]], FatJet_phi[FatJetidx[i]], FatJet_SDmass[FatJetidx[i]]);
            }
        }

        // if (jetidx.size() >= 2 && isBoosted2l2q == false) // FIXME: Only for testing purposes; comment this line and uncomment the next line for real analysis
        if (jetidx.size() >= 2)
        {
            foundZZCandidate = true;
            if (Z2.M() < 40.0 || Z2.M() > 250)
            {
                cut2l2j++;
            }
            cut2l1Jor2j++;

            TLorentzVector Z2_1;
            TLorentzVector Z2_2;
            Z2_1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
            Z2_2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);
            Z2_2j = Z2_1 + Z2_2;

            // Select the two jets with mass closest to Z-boson mass
            float mass_diff = 999.0;
            for (unsigned int i = 0; i < jetidx.size(); i++)
            {
                for (unsigned int j = i+1; j < jetidx.size(); j++) // FIXME: Check if there should be +1 or not
                {
                    TLorentzVector Z2_1;
                    TLorentzVector Z2_2;
                    Z2_1.SetPtEtaPhiM(Jet_pt[jetidx[i]], Jet_eta[jetidx[i]], Jet_phi[jetidx[i]], Jet_mass[jetidx[i]]);
                    Z2_2.SetPtEtaPhiM(Jet_pt[jetidx[j]], Jet_eta[jetidx[j]], Jet_phi[jetidx[j]], Jet_mass[jetidx[j]]);
                    TLorentzVector Z2_2j_temp = Z2_1 + Z2_2;

                    if (fabs(Z2_2j_temp.M() - Zmass) < mass_diff)
                    {
                        mass_diff = fabs(Z2_2j_temp.M() - Zmass);
                        Z2 = Z2_2j_temp;
                        resolvedJet1_Index = jetidx[i];
                        resolvedJet2_Index = jetidx[j];
                    }
                }
            }

            if (DEBUG)
                std::cout << "Z2: Mass based, pT based: " << Z2.Pt() << ",  " << Z2_2j.Pt() << std::endl;
        }

        ZZsystem = Z1 + Z2;
        ZZsystemnofsr = Z1nofsr + Z2; // FIXME: Update this with jet information.

        ZZ_2jsystem = Z1 + Z2_2j;
        ZZ_2jsystemnofsr = Z1nofsr + Z2_2j;

        float massZZ;
        massZZ = ZZsystem.M();

        float massZZ_2j;
        massZZ_2j = ZZ_2jsystem.M();
    }

    if(foundZZCandidate == false){
        return foundZZCandidate;
    }

   return foundZZCandidate;
}

bool H4LTools::ZZSelection_2l2nu()
{
    bool foundZZCandidate = false;
    jetidx = SelectedJets(tighteleforjetidx, tightmuforjetidx);
    if (DEBUG)
        std::cout << "Number of jets: " << jetidx.size() << std::endl;

    // No b-tagged jets
    if (!(HZZ2l2qNu_nTightBtagJets > 0))
    {
        HZZ2l2nu_cutbtag++;
    }

    if (DEBUG)
        std::cout << "Number of b-tagged jets [inside 2l2nu]: (Tight, Med, Loose): " << HZZ2l2qNu_nTightBtagJets << ", " << HZZ2l2qNu_nMediumBtagJets << ", " << HZZ2l2qNu_nLooseBtagJets << std::endl;

    // Get Angle between MET and nearest good jet
    for (unsigned int i = 0; i < jetidx.size(); i++)
    {
        float dPhi = fabs(TVector2::Phi_mpi_pi(Jet_phi[jetidx[i]] - MET_phi));
        if (dPhi < minDeltaPhi)
        {
            minDeltaPhi = dPhi;
        }
    }
    // If there are no jets, set the kinematics of ZZ system
    if (minDeltaPhi < HZZ2l2nu_dPhi_jetMET)
    {
        return foundZZCandidate;
    }

    // As we are not adding cut on the b-tag but we want the cut-flow table with the b-tag cut,
    // we are adding this cut here
    if (!(HZZ2l2qNu_nTightBtagJets > 0))
        HZZ2l2nu_cutdPhiJetMET++;

    if (DEBUG)
    {
        std::cout << "Passed dPhiJetMET cut" << std::endl;
        std::cout << "MET_pt: " << MET_pt << std::endl;
    }

    // As we are not adding cut on the b-tag but we want the cut-flow table with the b-tag cut,
    // we are adding this cut here
    if (HZZ2l2qNu_nTightBtagJets == 0 && MET_pt > 100)
    {
        HZZ2l2nu_cutMETgT100++;
    }

    Z2_met.SetPtEtaPhiE(MET_pt, 0, MET_phi, MET_pt);

    ZZ_metsystem = Z1 + Z2_met;
    ZZ_metsystemnofsr = Z1nofsr + Z2_met;

    // float Pz_nu;
    // float Pz_neutrino;
    // Pz_nu = ((Z1.M()*Z1.M())/4)- (MET_pt*MET_pt);
    // //if (Pz_nu < 0) {
    // std::complex<double> complex_pz(0, std::sqrt(-1 * Pz_nu));
    // Pz_neutrino = std::abs(complex_pz);
    // //}
    // if (DEBUG) {
    // std::cout << "Pz of neutrino in COM frame: " << Pz_neutrino << std::endl;
    // }

    float METZZ_met;
    METZZ_met = ZZ_metsystem.E();

    float MT_2l2nu;
    MT_2l2nu = ZZ_metsystem.Mt();

    // if the event passes all the cuts, set the flag `foundZZCandidate` to true
    foundZZCandidate = true;

    // Fetch number of jets
    HZZ2l2qNu_nJets = jetidx.size();
    if (DEBUG)
        std::cout << "Size of jets: [inside 2l2nu] " << HZZ2l2qNu_nJets << std::endl;

    // Get VBF jets having dEta>4.0 and mjj>500
    // If there are more than one pair of VBF jets, select the pair with highest mjj
    float VBF_jj_mjj = 0.0;
    for (unsigned int i = 0; i < jetidx.size(); i++)
    {
        for (unsigned int j = i + 1; j < jetidx.size(); j++)
        {
            // std::cout << "Inside: i: " << i << ", j: " << j << std::endl;
            TLorentzVector VBF_jet1;
            TLorentzVector VBF_jet2;
            VBF_jet1.SetPtEtaPhiM(Jet_pt[jetidx[i]], Jet_eta[jetidx[i]], Jet_phi[jetidx[i]], Jet_mass[jetidx[i]]);
            VBF_jet2.SetPtEtaPhiM(Jet_pt[jetidx[j]], Jet_eta[jetidx[j]], Jet_phi[jetidx[j]], Jet_mass[jetidx[j]]);
            TLorentzVector VBF_jj = VBF_jet1 + VBF_jet2;

            // FIXME: Move these hardcoded cuts to the YAML configuration file
            if (fabs(VBF_jet1.Eta() - VBF_jet2.Eta()) > 4.0 && VBF_jj.M() > 500.0)
            {
                if (VBF_jj.M() > VBF_jj_mjj)
                {
                    VBF_jj_mjj = VBF_jj.M();
                    HZZ2l2nu_VBFIndexJet1 = jetidx[i];
                    HZZ2l2nu_VBFIndexJet2 = jetidx[j];
                    if (DEBUG)
                        std::cout << "Inside: VBF_jj_mjj: " << VBF_jj_mjj << "\tHZZ2l2nu_VBFIndexJet1: " << HZZ2l2nu_VBFIndexJet1 << "\tHZZ2l2nu_VBFIndexJet2: " << HZZ2l2nu_VBFIndexJet2 << std::endl;
                }
            }
        }
    }

    if (jetidx.size() >= 2 && HZZ2l2nu_VBFIndexJet1 >= 0 && HZZ2l2nu_VBFIndexJet2 >= 0)
    {
        HZZ2l2nu_ifVBF = true;
        TLorentzVector VBF_jet1;
        TLorentzVector VBF_jet2;
        VBF_jet1.SetPtEtaPhiM(Jet_pt[HZZ2l2nu_VBFIndexJet1], Jet_eta[HZZ2l2nu_VBFIndexJet1], Jet_phi[HZZ2l2nu_VBFIndexJet1], Jet_mass[HZZ2l2nu_VBFIndexJet1]);
        VBF_jet2.SetPtEtaPhiM(Jet_pt[HZZ2l2nu_VBFIndexJet2], Jet_eta[HZZ2l2nu_VBFIndexJet2], Jet_phi[HZZ2l2nu_VBFIndexJet2], Jet_mass[HZZ2l2nu_VBFIndexJet2]);
        TLorentzVector VBF_jj = VBF_jet1 + VBF_jet2;
        if (DEBUG)
            std::cout << "Outside: VBF_jj_mjj: " << VBF_jj.M() << "\tHZZ2l2nu_VBFIndexJet1: " << HZZ2l2nu_VBFIndexJet1 << "\tHZZ2l2nu_VBFIndexJet2: " << HZZ2l2nu_VBFIndexJet2 << std::endl;
    }

    return foundZZCandidate;
}

bool H4LTools::GetWW_lnuqq()
{
    if (DEBUG)
        std::cout << "Inside function GetWW_lnuqq()" << std::endl;
    bool foundZ1Candidate = false;
    bool foundZZCandidate = false;

        if (!(nTightMu == 1 || nTightEle == 1))
    {
        return foundZ1Candidate;
    }
    if (DEBUG)
        std::cout << "Number of leptons: (Mu, Ele, Total): " << nTightMu << ", " << nTightEle << ", " << nTightMu + nTightEle << std::endl;
    // Set HZZ2l2qNu_isELE to true if there are 2 tight electrons, false if there are 2 tight muons
    HZZ2l2qNu_isELE = true ? (nTightEle == 1) : false;
    HWWlNu2q_cut1l++;

    if (DEBUG)
        std::cout << "nTightEleChgSum: " << nTightEleChgSum << "\tnTightMuChgSum: " << nTightMuChgSum << std::endl;

    jetidx = SelectedJets(tighteleforjetidx, tightmuforjetidx);
    FatJetidx = SelectedFatJets(tighteleforjetidx, tightmuforjetidx);

    if (FatJetidx.size() > 0 || jetidx.size() >= 2)
    {
        if (FatJetidx.size() > 0)
        {
            for (unsigned int i = 0; i < FatJetidx.size(); i++)
            {
                if (FatJet_PNZvsQCD[FatJetidx[i]] < 0.9)
                    continue;
                if (FatJet_pt[FatJetidx[i]] < 200.0)
                    continue;
                // if (FatJet_msoftdrop[FatJetidx[i]] < 40.0) continue;
                // if (FatJet_msoftdrop[FatJetidx[i]] > 180.0) continue;

                foundZZCandidate = true;
                isBoosted2l2q = true;
                cut2l1J++;
                cut2l1Jor2j++;

                boostedJet_PNScore = FatJet_PNZvsQCD[FatJetidx[i]];
                boostedJet_Index = FatJetidx[i];

                Z2.SetPtEtaPhiM(FatJet_pt[FatJetidx[i]], FatJet_eta[FatJetidx[i]], FatJet_phi[FatJetidx[i]], FatJet_SDmass[FatJetidx[i]]);
            }
        }

        // if (jetidx.size() >= 2 && isBoosted2l2q == false) // FIXME: Only for testing purposes; comment this line and uncomment the next line for real analysis
        if (jetidx.size() >= 2)
        {
            foundZZCandidate = true;
            if (Z2.M() < 40.0 || Z2.M() > 250)
            {
                cut2l2j++;
            }
            cut2l1Jor2j++;

            TLorentzVector Z2_1;
            TLorentzVector Z2_2;
            Z2_1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
            Z2_2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);
            Z2_2j = Z2_1 + Z2_2;

            // Select the two jets with mass closest to Z-boson mass
            float mass_diff = 999.0;
            for (unsigned int i = 0; i < jetidx.size(); i++)
            {
                for (unsigned int j = i + 1; j < jetidx.size(); j++) // FIXME: Check if there should be +1 or not
                {
                    TLorentzVector Z2_1;
                    TLorentzVector Z2_2;
                    Z2_1.SetPtEtaPhiM(Jet_pt[jetidx[i]], Jet_eta[jetidx[i]], Jet_phi[jetidx[i]], Jet_mass[jetidx[i]]);
                    Z2_2.SetPtEtaPhiM(Jet_pt[jetidx[j]], Jet_eta[jetidx[j]], Jet_phi[jetidx[j]], Jet_mass[jetidx[j]]);
                    TLorentzVector Z2_2j_temp = Z2_1 + Z2_2;

                    if (fabs(Z2_2j_temp.M() - Zmass) < mass_diff)
                    {
                        mass_diff = fabs(Z2_2j_temp.M() - Zmass);
                        Z2 = Z2_2j_temp;
                        resolvedJet1_Index = jetidx[i];
                        resolvedJet2_Index = jetidx[j];
                    }
                }
            }

            if (DEBUG)
                std::cout << "L#1732: Z2: Mass based, pT based: " << Z2.Pt() << ",  " << Z2_2j.Pt() << std::endl;
        }

        if (DEBUG)
        {
            std::cout << "L#1737: Size of ElelistFsr: " << ElelistFsr.size() << std::endl;
            std::cout << "L#1738: Size of MulistFsr: " << MulistFsr.size() << std::endl;
            std::cout << "L#1739: Size of TightEleindex: " << TightEleindex.size() << std::endl;
            std::cout << "L#1740: Size of TightMuindex: " << TightMuindex.size() << std::endl;
        }
        // Z1 will be the sum of TLorentzVector of lepton and MET
        TLorentzVector Lep1;
        if (HZZ2l2qNu_isELE)
        {
            Lep1.SetPtEtaPhiM(ElelistFsr[TightEleindex[0]].Pt(), ElelistFsr[TightEleindex[0]].Eta(), ElelistFsr[TightEleindex[0]].Phi(), ElelistFsr[TightEleindex[0]].M());
        }
        else
        {
            Lep1.SetPtEtaPhiM(MulistFsr[TightMuindex[0]].Pt(), MulistFsr[TightMuindex[0]].Eta(), MulistFsr[TightMuindex[0]].Phi(), MulistFsr[TightMuindex[0]].M());
        }
        TLorentzVector MET;
        MET.SetPtEtaPhiM(MET_pt, 0, MET_phi, 0);
        Z1 = Lep1 + MET;

        ZZsystem = Z1 + Z2;
        ZZsystemnofsr = Z1nofsr + Z2; // FIXME: Update this with jet information.

        ZZ_2jsystem = Z1 + Z2_2j;
        ZZ_2jsystemnofsr = Z1nofsr + Z2_2j;

    }

    return foundZZCandidate;
}

float H4LTools::getDg4Constant(float ZZMass){
    return spline_g4->Eval(ZZMass);
}

float H4LTools::getDg2Constant(float ZZMass){
    return spline_g2->Eval(ZZMass);
}

float H4LTools::getDL1Constant(float ZZMass){
    return spline_L1->Eval(ZZMass);
}

float H4LTools::getDL1ZgsConstant(float ZZMass){
    return spline_L1Zgs->Eval(ZZMass);
}
