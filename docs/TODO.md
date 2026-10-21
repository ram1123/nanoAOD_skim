# To-do

- [ ] **Electron-muon cross-cleaning** — disabled on 2026-10-20 to match AN-2016/325, which has no
  such cleaning. The code is still there: `H4LTools::FlagElectronMuonOverlap()`, YAML
  `Electron.MuonCleaning` (dR 0.05, muon pT > 5 GeV). Re-enable it by restoring the call at the start of
  `H4LTools::LeptonSelection` (and the summary print in `H4LCppModule.endJob`). Before deciding,
  compare yields with and without it, in particular electrons that are really muon tracks
  feeding the extra-lepton veto and the eμ CR.
