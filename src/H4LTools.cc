#include "../include/H4LTools.h"
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <vector>


std::vector<unsigned int> H4LTools::goodLooseElectrons2012(){
    std::vector<unsigned int> LooseElectronindex;
    for (unsigned int i=0; i<Electron_pt.size(); i++){
        if (DEBUG)
            std::cout << "Inside goodLooseElectrons2012:: Electron_pt[" << i << "] = " << Electron_pt[i] << std::endl;
        //if ((Electron_pt[i]>elePtcut)&&(fabs(Electron_eta[i])<eleEtacut)){
        if ((Electron_pt[i]>elePtcut)&&(((fabs(Electron_eta[i])<1.4442)||(fabs(Electron_eta[i])>1.5660))&&(fabs(Electron_eta[i])<eleEtacut))){
            LooseElectronindex.push_back(i);
            //std::cout << nElectron << std::endl;
        }
    }

    return LooseElectronindex;
}

std::vector<unsigned int> H4LTools::goodLooseMuons2012(){
    std::vector<unsigned int> LooseMuonindex;
    for (unsigned int i=0; i<Muon_eta.size(); i++){
        if (DEBUG)
            std::cout << "Inside goodLooseMuons2012:: Muon_pt[" << i << "] = " << Muon_pt[i] << std::endl;
        if ((Muon_pt[i]>MuPtcut)&&(fabs(Muon_eta[i])<MuEtacut)&&((Muon_isGlobal[i]||Muon_isTracker[i]||Muon_isPFcand[i])&&(Muon_mediumId[i]))){
            LooseMuonindex.push_back(i);
      //      std::cout << nMuon << std::endl;
        }
    }

    return LooseMuonindex;
}
std::vector<unsigned int> H4LTools::goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex){
    std::vector<unsigned int> bestMuonindex;
    for (unsigned int i=0; i<Muonindex.size(); i++){
        if ((Muon_pt[Muonindex[i]]>MuPtcut)&&(fabs(Muon_eta[Muonindex[i]])<MuEtacut)&&(Muon_isGlobal[Muonindex[i]]||Muon_isTracker[Muonindex[i]])&&(Muon_mediumId[Muonindex[i]])){
            //if (Muon_sip3d[Muonindex[i]]<Musip3dCut){
                if((fabs(Muon_dxy[Muonindex[i]])<MuLoosedxycut)&&(fabs(Muon_dz[Muonindex[i]])<MuLoosedzcut)){
                    bestMuonindex.push_back(Muonindex[i]);
               // }
            }
        }
    }
    return bestMuonindex;
}
std::vector<unsigned int> H4LTools::goodElectrons2015_noIso_noBdt(std::vector<unsigned int> Electronindex){
    std::vector<unsigned int> bestElectronindex;
    for (unsigned int i=0; i<Electronindex.size(); i++){
        if ((Electron_pt[Electronindex[i]])>elePtcut){
           // if(Electron_sip3d[Electronindex[i]]<elesip3dCut){
                if((fabs(Electron_dxy[Electronindex[i]])<eleLoosedxycut)&&(fabs(Electron_dz[Electronindex[i]])<eleLoosedzcut)){
                    bestElectronindex.push_back(Electronindex[i]);
                    //std::cout << nElectron << std::endl;
                //}
            }
        }
    }

    return bestElectronindex;
}
std::vector<bool> H4LTools::passTight_BDT_Id(){
    std::vector<bool> tightid;
    float cutVal,mvaVal;
    cutVal = 1000;
    mvaVal = -1;
    //unsigned nE = (*nElectron).Get()[0];
    for (unsigned int i=0; i<Electron_pt.size(); i++){

        mvaVal = Electron_mvaFall17V2Iso_WP90[i];
//        if(mvaVal > cutVal){
            tightid.push_back(mvaVal);
            //std::cout << nElectron << std::endl;
   //     }
     //   else{
       //     tightid.push_back(false);
     //   }

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
}

void H4LTools::LeptonSelection(){
    looseEle = goodLooseElectrons2012();
    looseMu = goodLooseMuons2012();
    bestEle = goodElectrons2015_noIso_noBdt(looseEle);
    bestMu = goodMuons2015_noIso_noPf(looseMu);
    Electronindex = bestEle;
    Muonindex = bestMu;
    AllEid = passTight_BDT_Id();
    AllMuid = passTight_Id();
    for (unsigned int iuj=0;iuj<looseEle.size();iuj++){
        if(AllEid[looseEle[iuj]]) tighteleforjetidx.push_back(looseEle[iuj]);
    }
    for (unsigned int juj=0;juj<looseMu.size();juj++){
        if(AllMuid[looseMu[juj]]) tightmuforjetidx.push_back(looseMu[juj]);
    }
    jetidx = SelectedJets(tighteleforjetidx,tightmuforjetidx);
    FatJetidx = SelectedFatJets(tighteleforjetidx, tightmuforjetidx);

    for(unsigned int ie=0; ie<Electronindex.size();ie++){
        if(Electron_pdgId[Electronindex[ie]]>0){
            Elechg.push_back(-1);
        }
        else{
            Elechg.push_back(1);
        }
        TLorentzVector Ele;
        if (DEBUG)
            std::cout << "Inside LeptonSelection:: Electron_pt[" << Electronindex[ie] << "] = " << Electron_pt[Electronindex[ie]] << std::endl;
        Ele.SetPtEtaPhiM(Electron_pt[Electronindex[ie]],Electron_eta[Electronindex[ie]],Electron_phi[Electronindex[ie]],Electron_mass[Electronindex[ie]]);
        Elelist.push_back(Ele);
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
        TLorentzVector Mu;
        if (DEBUG)
            std::cout << "Inside LeptonSelection:: Muon_pt[" << Muonindex[imu] << "] = " << Muon_pt[Muonindex[imu]] << std::endl;
        Mu.SetPtEtaPhiM(Muon_pt[Muonindex[imu]],Muon_eta[Muonindex[imu]],Muon_phi[Muonindex[imu]],Muon_mass[Muonindex[imu]]);
        Mulist.push_back(Mu);
        muid.push_back(AllMuid[Muonindex[imu]]);
        Muiso.push_back(Muon_pfRelIso03_all[Muonindex[imu]]);
    }

    ElelistFsr = BatchFsrRecovery(Elelist);
    MulistFsr = BatchFsrRecovery(Mulist);

    for(unsigned int ae=0; ae<Eid.size();ae++){
        float RelEleIsoNoFsr;
        RelEleIsoNoFsr = Eiso[ae];
        if (isFSR){
          unsigned int FsrEleidx;
          FsrEleidx = doFsrRecovery(Elelist[ae]);
          if(FsrEleidx<900){
              TLorentzVector fsrele;
              fsrele.SetPtEtaPhiM(FsrPhoton_pt[FsrEleidx],FsrPhoton_eta[FsrEleidx],FsrPhoton_phi[FsrEleidx],0);
              if (DEBUG)
                std::cout<<"Ele correction: "<< std::endl;
              if(Elelist[ae].DeltaR(fsrele)>0.01){
                RelEleIsoNoFsr = RelEleIsoNoFsr - FsrPhoton_pt[FsrEleidx]/Elelist[ae].Pt();
              }
          }

        }
        if((Eid[ae]==true)&&(RelEleIsoNoFsr<0.35)){
            nTightEle++;
            TightEleindex.push_back(ae);
            nTightEleChgSum += Elechg[ae];
        }

    }

    for(unsigned int amu=0; amu<muid.size();amu++){
        float RelIsoNoFsr;
        RelIsoNoFsr = Muiso[amu];
        if (isFSR){
          unsigned int FsrMuonidx;
          FsrMuonidx = doFsrRecovery(Mulist[amu]);
          if(FsrMuonidx<900){
              TLorentzVector fsrmuon;
              fsrmuon.SetPtEtaPhiM(FsrPhoton_pt[FsrMuonidx],FsrPhoton_eta[FsrMuonidx],FsrPhoton_phi[FsrMuonidx],0);
              if (DEBUG)
                std::cout<<"muon FSR recovered"<<endl;
              if(Mulist[amu].DeltaR(fsrmuon)>0.01){
                RelIsoNoFsr = RelIsoNoFsr - FsrPhoton_pt[FsrMuonidx]/Mulist[amu].Pt();
              }
          }
        }
        if((muid[amu]==true)&&(RelIsoNoFsr<0.35)){
            nTightMu++;
            TightMuindex.push_back(amu);
            nTightMuChgSum += Muchg[amu];
        }
    }


}
bool H4LTools::findZCandidate(){

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

    if(TightEleindex.size()>1){
        for(unsigned int ke=0; ke<(TightEleindex.size()-1);ke++){
            for(unsigned int je=ke+1;je<TightEleindex.size();je++){
                // FIXME: Removed the charge requirement for the Z candidate. It should be propagated correctly for the ZZ->4l analysis
                // if ((Elechg[TightEleindex[ke]]+Elechg[TightEleindex[je]])==0)
                {
                    TLorentzVector Zcan;
                    Zcan = ElelistFsr[TightEleindex[ke]] + ElelistFsr[TightEleindex[je]];
                    if((Zcan.M()>MZcutdown)&&(Zcan.M()<MZcutup)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightEleindex[ke]);
                        Zlep2index.push_back(TightEleindex[je]);
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
        }
    }

    if(TightMuindex.size()>1){
        for(unsigned int kmu=0; kmu<(TightMuindex.size()-1);kmu++){
            for(unsigned int jmu=kmu+1;jmu<TightMuindex.size();jmu++){
                // if ((Muchg[TightMuindex[kmu]]+Muchg[TightMuindex[jmu]])==0)
                {
                    TLorentzVector Zcan;
                    Zcan = MulistFsr[TightMuindex[kmu]] + MulistFsr[TightMuindex[jmu]];
                    if (DEBUG)
                        std::cout << "Zcan.M() = " << Zcan.M() << "\tMZcutdown = " << MZcutdown << "\tMZcutup = " << MZcutup << std::endl;

                    if((Zcan.M()>MZcutdown)&&(Zcan.M()<MZcutup)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightMuindex[kmu]);
                        Zlep2index.push_back(TightMuindex[jmu]);
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

    if (Zsize > 0)
    {
        return true;
    }
    else{
        return false;
    }
}


bool H4LTools::ZZSelection_4l(){

    bool foundZZCandidate = false;
    //std::cout << " Inside the 4l loop in .cc file" << std::endl;
    if(!findZCandidate()){
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

bool H4LTools::GetZ1_2l2qOR2l2nu()
{
    if (DEBUG)
        std::cout << "Inside function GetZ1_2l2qOR2l2nu()" << std::endl;
    bool foundZ1Candidate = false;
    if (!findZCandidate())
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
    if (std::abs(nTightEleChgSum) != 0 && std::abs(nTightMuChgSum) != 0)
    {
        HZZ2l2qNu_cutOppositeCharge++;
        HZZ2l2qNu_cutOppositeChargeFlag = true;
    }

    if (DEBUG)
        std::cout << "Zlist size (Before): " << Zlist.size() << std::endl;
    // if (!(Zlist.size() == 1)) // There should be exactly 1 Z candidate
    // {
    //     return foundZ1Candidate;
    // }
    if (DEBUG)
        std::cout << "Zlist size (After): " << Zlist.size() << std::endl;

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
    foundZ1Candidate = true;
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

    if(Lep1.DeltaR(Lep2)<0.3){
    foundZ1Candidate = false;
    }

    jetidx = SelectedJets(tighteleforjetidx, tightmuforjetidx);
    if (DEBUG)
        std::cout << "Number of jets: " << jetidx.size() << std::endl;
    HZZ2l2qNu_nJets = jetidx.size();

    // count the number of tight, medium and loose b-tagged jets
    for (unsigned int i = 0; i < jetidx.size(); i++) // FIXME: These variables seems to be wrong.
    {
        // Reference: https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/#ak4-b-tagging
        if (DEBUG)
        {
            std::cout << "Jet_btagDeepFlavB[" << jetidx[i] << "]: " << Jet_btagDeepFlavB[jetidx[i]] << std::endl;
        }
        if (Jet_btagDeepFlavB[jetidx[i]] > btag_deepJet_Tight)
            HZZ2l2qNu_nTightBtagJets++;
        if (Jet_btagDeepFlavB[jetidx[i]] > btag_deepJet_Medium)
            HZZ2l2qNu_nMediumBtagJets++;
        if (Jet_btagDeepFlavB[jetidx[i]] > btag_deepJet_Loose)
            HZZ2l2qNu_nLooseBtagJets++;
    }
//    HZZ2l2qNu_cutmZ1Window++;

    return foundZ1Candidate;
}

////// emu control region (Z1 candidate selection) /////////
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
            if (Z2.M() < 40.0 || Z2.M() > 180)
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
    //if (HZZ2l2qNu_nMediumBtagJets > 0)
    //{
        //return foundZZCandidate;
    //}
    //HZZ2l2nu_cutbtag++;
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
    foundZZCandidate = true;
    HZZ2l2nu_cutdPhiJetMET++;
    if (DEBUG)
    {
        std::cout << "Passed dPhiJetMET cut" << std::endl;
        std::cout << "MET_pt: " << MET_pt << std::endl;
    }

    if (MET_pt > 100)
    {
        HZZ2l2nu_cutMETgT100++;
    }

    Z2_met.SetPtEtaPhiE(MET_pt, 0, MET_phi, MET_pt);

    ZZ_metsystem = Z1 + Z2_met;
    ZZ_metsystemnofsr = Z1nofsr + Z2_met;

    float Pz_nu;
    float Pz_neutrino;
    Pz_nu = ((Z1.M()*Z1.M())/4)- (MET_pt*MET_pt);
    //if (Pz_nu < 0) {
        std::complex<double> complex_pz(0, std::sqrt(-1 * Pz_nu));
        Pz_neutrino = std::abs(complex_pz);
    //}
    if (DEBUG) {
    std::cout << "Pz of neutrino in COM frame: " << Pz_neutrino << std::endl;
    }

    float METZZ_met;
    METZZ_met = ZZ_metsystem.E();

    float MT_2l2nu;
    MT_2l2nu = ZZ_metsystem.Mt();
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
