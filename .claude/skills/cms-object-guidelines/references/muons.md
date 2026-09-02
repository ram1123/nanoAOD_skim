# Muons — Stored CMS Recommendations

> ## ⚠ Ported from the H→μμ (copperhead) analysis — re-verify before use
>
> This repo is the **H→ZZ→4l / 2l2q / 2l2nu NanoAOD-tools skim**
> (`post_proc.py` → `modules/H4LCppModule.py` → C++ `src/H4LTools.cc`; cut values in
> `config/Input_<year>.yml`). It is **not** the coffea `src/copperhead_processor.py` /
> `configs/parameters/*.yaml` / stage-1/stage-2 pipeline this file was written for —
> ignore every reference to those paths, `src/corrections/*`, `run_stage*.py`, etc.
> - Eras here: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2023/2024/2025-specific rows below do **not** apply.
> - The CMS-POG recommendation content below is retained as a starting point and
>   **must be re-checked against this analysis's H→ZZ note / HIG group** before it is
>   treated as a requirement.
> - This-repo pointers: cuts in `config/Input_<year>.yml` `Muon:` → `H4LTools::InitializeMucut`; Rochester via nanoAOD-tools `muonScaleRes20XX` (`Muon_corrected_pt`), no beam-spot-constrained pT; FSR recovery in `H4LTools::MuonFsr`. **Primary** selection object in every channel.


Responsible POG: **MUO** (Muon POG). Momentum scale/resolution: MUO Rochester /
scale‑smearing group. Beam‑spot–constrained pT: this analysis's custom NanoAODv12.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S0 | **CMS AN‑2016/325** §4.3 Table 6 (muon "tight arbitration" ID + iso Eq. 2), §4.4 (2l2nu pre‑selection, soft/loose 3rd‑lepton veto) | 2l2nu analysis note → `references/hzz-2l2nu.md` §2 + §2.3 below | 2016/legacy — **[Verify]** for UL/Run 3 |
| S1 | H→µµ analysis note AN‑19‑124, §3.2–3.3 (muon selection, momentum, FSR) | skill‑author library: `AN2019_124_v7.pdf` — **not in repo** | local review 2026‑08‑31 |
| S2 | `CERN-THESIS-2021-201`, §4.2.1 + ch. 5 | skill‑author library — **not in repo** | local review 2026‑08‑31 |
| S3 | Beam‑spot–constrained muons for H→µµ | `MuonswithBSconstraintHmumu2_21.pdf` — **not in repo** | local review 2026‑08‑31 |
| S4 | HIG‑19‑006 TWiki snapshot | `CMS_HIG19006_twiki.pdf` — **not in repo** | local review 2026‑08‑31 |
| C1 | Per‑year selection parameters | `configs/parameters/muon.yaml` | 2026‑08‑31 |
| C2 | Trigger list | `configs/parameters/trigger.yaml` | 2026‑08‑31 |
| C3 | Correction payload paths | `configs/parameters/correction_filelist.yaml` (`roccor_file`, `BS_res_calib_path`) | 2026‑08‑31 |
| C4 | Scale‑factor payloads + map keys | `configs/parameters/SF_filelist.yaml` (`muSFFileList`) | 2026‑08‑31 |
| C5 | Implementation | `src/copperhead_processor.py` muon block (~L900–1130); `src/corrections/{rochester,MuonScaRe,geofit,fsr_recovery,muon_sf}.py` | 2026‑08‑31 |
| C6 | NanoAOD muon branch dictionary; Run 3 sync table | `docs/temp/muon_nanoAOD_docs.md`; `docs/Run3_all_basic_Information.md` | 2026‑08‑31 |
| P1 | MUO POG entry points | `MuonPOG#User_Recommendations`; Run 3 `MuonRun32022`, `MuonRun3_2023` (from PdmV `PdmVRun3Analysis` — see `lumi.md` §9) | via user, 2026‑09‑01 |

Not covered by stored sources → **Authoritative CMS verification required** (consult
P1 for the target era): the MUO numeric ID/iso working‑point recommendation per era,
the reco/tracking‑SF requirement, the Run 3 momentum‑calibration prescription, and the
correctionlib payload version.

Classification tags: **[MUO official]**, **[Analysis‑specific]** (AN‑19‑124 choice),
**[Implementation]** (code detail), **[Verify]** (not established from stored material).

---

## 1. Required context

Run 2 vs Run 3; exact era; data vs MC; NanoAOD version — **this repo uses a custom
NanoAODv12 for Run 2 and 2022/2023 to get `Muon_bsConstrainedPt` / `…Chi2`; NanoAODv15
for 2024+**; intended ID/iso working point; whether trigger SFs are needed.

---

## 2. Baseline selection

### 2.1 As implemented (C1, C5)

Base per‑muon selection, `copperhead_processor.py` ~L950:

| Cut | Value | Classification | Key / note |
|-----|-------|----------------|------------|
| pT | > **20 GeV** (all years) | analysis‑specific | `muon_pt_cut` — applied to `pt_raw` (pre‑Rochester); code carries a `FIXME: Why pt_raw` |
| \|η\| | < **2.4** | detector acceptance | `muon_eta_cut` (on `eta_raw`) |
| ID | **`mediumId`** | MUO cut‑based medium WP | `muon_id`, all years |
| track type | `isGlobal` **or** `isTracker` | AN‑19‑124 Table 3.5 | hardcoded |
| Isolation | `pfRelIso04_all` < **0.25** (after FSR) | loose MUO PF iso | `muon_iso_cut` |

Event‑level: **exactly 2** selected muons, **opposite charge**; number of good primary
vertices > 0; **≥ 1** muon trigger‑matched (§3).

### 2.2 AN‑19‑124 baseline items NOT in the current implementation  **[Verify]**

S1/S2 baseline also lists impact‑parameter cuts: `|dxy| < 0.05 cm`, `|dz| < 0.10 cm`,
`SIP3D < 8`. `copperhead_processor.py` (~L950) does **not** apply these — the
`mu*_dxy` / `dz` / `ip3d` / `sip3d` branches are only written as output columns.
Confirm whether the analysis deliberately drops the IP cuts or whether this is a
regression. (This repo's `config/Input_<year>.yml` `Muon:` **does** carry
`Loosedxycut` / `Loosedzcut` / `Tightdxycut` / `Tightdzcut` = 0.045 / 0.2 →
`InitializeMucut` — check they reach the worker.)

### 2.3 2l2nu — muons as the Z→µµ leg (AN‑2016/325 §4.3–§4.4)  **[HIG / AN‑2016‑325]**

For the **2l2nu** channel (current focus) reference numbers from AN‑2016/325
(2016/legacy — **[Verify]** vs the current UL/Run 2+3 note; `references/hzz-2l2nu.md`
§2):

| Requirement | AN‑2016/325 value | This repo |
|-------------|-------------------|-----------|
| `pT` | > **25 GeV** | `HZZ2l2nu.Leading/SubLeading_Lep_pT = 25` ✓ |
| `\|η\|` | < **2.4** | `HZZ2l2nu.Lep_eta = 2.5` — **[Repo divergence]** (applies 2.5 to muons too) |
| ID | "tight" muon (PF + global/tracker, tight arbitration, Table 6) | repo uses `mediumId` + `isGlobal\|isTracker` (`Muon.` cuts) — verify vs the HZZ recommendation |
| Isolation | tight PF iso, Δβ (Eq. 2: `[I_ch + max(I_nh + I_γ − 0.5·I_ch^PU, 0)]/pT`) | `Muon.Isocut = 0.2` (`pfRelIso04_all`) |
| dilepton | `\|m_{µµ} − 91\| < 15 GeV`, **not** charge‑required, reject if > 2 µ candidates | `M_ll_Window: 0.0` — **not applied** — **[Repo divergence]** |
| Z `pT` | `p_T^{µµ} > 55 GeV` | `Pt_ll: 10.0` — **[Repo divergence]** |
| 3rd‑lepton veto | loose muon `I_rel < 0.2`, **or** soft muon `pT > 3 GeV` | verify definitions |

AN‑2016/325 systematics on the µµ leg: trigger 2%, ID+iso 2% per 2µ event, muon
momentum scale 1% (propagated to MET), soft‑muon‑veto efficiency 96–99%
(`references/hzz-2l2nu.md` §8).

---

## 3. Trigger and trigger matching (C1, C2)

| Era | HLT paths | Matched‑muon offline pT | Classification |
|-----|-----------|--------------------------|----------------|
| 2016preVFP / 2016postVFP | `IsoMu24`, `IsoTkMu24` | > 26 GeV | analysis plateau choice |
| 2017 | `IsoMu27` | > 29 GeV | analysis plateau choice |
| 2018 | `IsoMu24` | > 26 GeV | analysis plateau choice |
| 2022preEE … 2026 | `IsoMu24` | > 26 GeV | analysis plateau choice |

Trigger‑match object requirements (`muon_trigmatch_*`): `tightId`, iso < 0.15, ΔR < 0.4,
matched trigger‑object pT > 24 GeV (27 GeV for 2017). Applied **after** Rochester and the
base selection, **before** FSR recovery (AN‑19‑124 L373).

---

## 4. Momentum corrections — order and payloads

Order (`copperhead_processor.py`): **beam‑spot constraint → Rochester → (trigger match) →
FSR recovery**. Data and simulation get different Rochester operations; never apply MC
smearing to data.

### 4.1 Beam‑spot–constrained pT  **[Analysis‑specific / Implementation]**

- `do_beamConstraint`: use `Muon_bsConstrainedPt`; **revert to the default muon pT if
  `Muon_bsConstrainedChi2 ≥ 30`** (C6).
- Overrides `do_geofit` (`src/corrections/geofit.py`) when both are enabled.
- EBE mass‑resolution / BS calibration JSONs: `BS_res_calib_path` (C3), per year, MC +
  Data. 2025/2026 reuse the 2024 calibration — **placeholder, [Verify]**.
- Treat as a resolution improvement, not a generic MUO selection: verify supported
  track types, data/MC treatment, covariance propagation, Z→µµ validation, and its own
  scale/resolution uncertainties before relying on it.

### 4.2 Rochester / MUO scale‑smearing (`roccor_file`, C3)  **[MUO official payload]**

| Era | Payload | Kind |
|-----|---------|------|
| 2016preVFP / 2016postVFP / 2017 / 2018 (UL) | `data/roch_corr/RoccoR20*UL.txt` | Rochester text |
| 2016 / 2017 / 2018 RERECO | `data/roch_corr/RoccoR20*.txt` | Rochester text |
| 2022preEE / 2022postEE / 2023 / 2023BPix | `data/roch_corr/20*_Summer2*.json` | MUO scale‑smearing JSON |
| 2024 | `data/roch_corr/2024_Summer24.json` (copied from GitLab) | JSON |
| 2025 | `data/roch_corr/2025_muon_scalesmearing_VXBS.json` | JSON |
| 2026 | reuses the 2025 VXBS JSON | **placeholder — no 2026 payload on cvmfs, [Verify]** |

`src/corrections/MuonScaRe.py` implements the Run 3 JSON; `src/corrections/rochester.py`
the Run 2 text.

### 4.3 FSR recovery  **[Analysis‑specific]**

`do_fsr` → `fsr_recoveryV1` (`src/corrections/fsr_recovery.py`), applied **after** trigger
matching. Recovered photon is added to the muon and its energy folded into
`pfRelIso04_all` **before** the isolation cut (AN‑19‑124 L360).

### 4.4 L1 prefiring  **[Implementation]**

`do_l1prefiring_wgts` → `L1PreFiringWeight.{Nom,Up,Dn}` as an event weight (Run 2 in
practice).

---

## 5. Scale factors (`muSFFileList`, C4)  **[MUO official payload]**

Payload **`muon_Z.json.gz`** (`jsonpog-integration` `POG/MUO` for Run 2 UL v9;
CAT `metadata/MUO/*` for Run 3). `src/corrections/muon_sf.py` applies ID·iso·trigger SFs
as event weights.

| Era | ID key | Iso key | Trigger key |
|-----|--------|---------|-------------|
| 2016preVFP / 2016postVFP | `NUM_MediumID_DEN_TrackerMuons` | `NUM_LooseRelIso_DEN_MediumID` | `NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight` |
| 2017 | `NUM_MediumID_DEN_TrackerMuons` | `NUM_LooseRelIso_DEN_MediumID` | `NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight` — *config notes an input‑binning FIXME* |
| 2018 | `NUM_MediumID_DEN_TrackerMuons` | `NUM_LooseRelIso_DEN_MediumID` | `NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight` |
| 2022preEE … 2025 | `NUM_MediumID_DEN_TrackerMuons` | `NUM_LoosePFIso_DEN_MediumID` | `NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium` |
| 2026 | same as 2025 | | **placeholder — MUO 2026 payload dir empty on cvmfs, [Verify]** |

- **Reco/tracking SF**: not applied — **[Verify]** whether the target era's MUO
  recommendation requires a separate reco SF for `mediumId` tracker muons.
- Confirm each SF's `DEN`/`NUM` definition matches the analysis selection (medium ID,
  loose PF iso, `IsoMu24`/`IsoMu27`) for the exact era, and keep stat + syst variations.

---

## 6. Review checklist

1. Era / data‑MC / NanoAOD version identified; custom‑v12 BS branches actually present.
2. pT (20), \|η\| (2.4), `mediumId`, `isGlobal|isTracker`, PF iso (0.25) applied per C1/C5.
3. IP cuts (dxy/dz/sip3d) — decision recorded (currently **not** applied, §2.2).
4. Matched‑muon pT threshold correct for the year (26 vs 29 GeV); HLT path matches C2.
5. Correction order = BS constraint → Rochester → FSR; `bsConstrainedChi2 ≥ 30` fallback.
6. Rochester / scale‑smearing payload matches era + NanoAOD version; 2025/2026
   placeholders understood.
7. FSR applied after trigger match; FSR energy in `pfRelIso04_all` before the iso cut.
8. SF keys (ID/iso/trigger) match the selection and era; 2026 placeholder understood;
   reco‑SF need resolved.
9. Systematics: SF up/down, Rochester/scale‑smearing variations, (Run 2) L1 prefiring,
   BS‑constraint scale/resolution.

---

## 7. Cross‑check vs this repo's config (as of 2026‑08‑31)

| Observation | Detail |
|-------------|--------|
| IP cuts missing | AN‑19‑124 baseline has `|dxy|<0.05`, `|dz|<0.10`, `SIP3D<8`; not applied in `copperhead_processor.py` (~L950) — only stored |
| pT cut on `pt_raw` | `muon_pt_cut` is compared against pre‑Rochester `pt_raw`; code has `FIXME: Why pt_raw` |
| 2025 / 2026 placeholders | Rochester JSON, BS/EBE calibration, and `muon_Z.json.gz` for 2026 all reuse earlier years |
| 2017 trigger SF | `muSFFileList["2017"]` marked `FIXME: input binning error` |

---

## 8. Evidence summary

| Item | POG | Eras | Source | Established? |
|------|-----|------|--------|--------------|
| Selection (pT/η/ID/iso, track type, dimuon, trigger match) | analysis | all | S1, C1, C5 | yes — analysis choice, not a POG number |
| IP cuts (dxy/dz/sip3d) | analysis | all | S1, S2 | in note; **not implemented** — [Verify] |
| Correction order BS→Rochester→FSR | analysis (AN‑19‑124) | all | S1, C5 | yes |
| Rochester / MUO scale‑smearing payloads | MUO | all | C3 | paths yes; 2025/2026 placeholders |
| ID/iso/trigger SF keys | MUO | all | C4 | keys yes; correctness vs era **[Verify]** |
| Reco/tracking SF requirement | MUO | all | — | **Authoritative CMS verification required** |
| `mediumId` = recommended ID for H→µµ in the target era | MUO | all | — | **Authoritative CMS verification required** |
| Trigger pT thresholds vs official HLT plateau per era | MUO/HLT | all | C1 | **Authoritative CMS verification required** |
| Run 3 momentum calibration / correctionlib version | MUO | Run 3 | — | **Authoritative CMS verification required** |

## Last verified

- Local source review: 2026‑08‑31
- Current POG recommendation: pending
