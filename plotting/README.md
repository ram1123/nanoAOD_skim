# 2017 control plots with the AN-2016/325 selection — HZZ 2l2nu, UL2017 NanoAOD v15

Moved here from `PlottingCodes/ControlPlots/Run2017_v15_AN325/` (Anusree's PlottingCodes repo). The
plotting engine `myControlPlots.C` and its helpers (`tdrstyle.C`, `utils.C`, `CMS_lumi.*`, `plotvar_t.h`)
are copies of `PlottingCodes/ControlPlots/`, so this directory runs without that clone. Machinery:
`myControlPlots.C`, `overallEventWeight`, `xsec × lumi / Σ genEventSumw`, 41.48 fb⁻¹.

- The ntuples come from a **re-skim** with the selection synced to AN-2016/325 (branch
  `HZZ_Analysis_2l2nu_dev_gammaCR`), see "Skim changes" below.
- The **AN-2016/325 Sec. 4.4 pre-selection** is applied in every region.

Summary of both ANs: [`../docs/AN_summary_HZZ2l2nu.md`](../docs/AN_summary_HZZ2l2nu.md).

## Pipeline

| step | script | output |
|---|---|---|
| 1. slim + merge per sample, Σ genEventSumw, data de-duplication | [`../scripts/analysis/prepare_inputs.py`](../scripts/analysis/prepare_inputs.py) | `$BASE/inputs/` |
| 2. sample list, cut files, variable header (edit xsecs / cuts / ranges **here**) | `make_config.py` | `HZZ2l2nu_Samples13TeV_2017_v15.txt`, `cuts/`, `controlplotvars_2017.h` (generated, git-ignored) |
| 3. plots per region | `RunMacro_2017_v15.C` → `myControlPlots.C` | `$BASE/plots/` |
| 4. yields + α | `yield_table.py` | `yields_2017.md` |
| 5. provenance | `write_metadata.sh` | `$BASE/{inputs,plots}/_metadata/` |

`run_all.sh` runs all five.

## Where everything is

| what | path |
|---|---|
| skims | `/depot/cms/users/shar1172/HZZ2l2nu_skim/skims/HZZ2l2nu_AN325_2017_MCdata/<sample>/` |
| skim metadata (git state + diff, YAML, sample list, tasks, staged-copy check) | `…/skims/HZZ2l2nu_AN325_2017_MCdata/_metadata/<date_time>/` |
| slim plotting inputs | `/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325/inputs/` |
| **plots** | `/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325/plots/<region>/<category>/` |
| plot / input metadata | `…/controlplots_2017_AN325/{plots,inputs}/_metadata/` (written by `write_metadata.sh`) |
| scripts | this directory |

## Selection

UL / v15 object IDs are kept: MVA-iso WP90 electrons, medium-ID muons, PuppiMET with the φ correction,
UParTAK4 b tagging, and the Run-2 UL TightLepVeto jet ID. Event-level cuts and categories follow AN-2016/325.

**Pre-selection in every region** (AN Sec. 4.4):
- trigger: `Triggers_HZZ2l2nu_SR` (ee/μμ) or `_EMuCR` (eμ);
- 2 tight leptons with p_T > 25;
- **no charge requirement** (AN 4.3);
- no extra loose lepton (skim) and no soft muon with p_T > 3 GeV (`HZZ2l2nu_nExtraSoftMu==0`);
- p_T^ll > 55;
- **min Δφ(jet, p_T^miss) > 0.5** (jets with p_T > 30);
- **|Δφ(Z, p_T^miss)| > 0.5**;
- b-veto with the **UParTAK4 Loose** WP (the AN uses CSVv2 Loose).

**Blinding:** the ee/μμ signal region (on-Z, p_T^miss > 125) is in none of the regions.

| region | flavours | extra cuts | categories |
|---|---|---|---|
| `LowMET` | ee, μμ | \|m_ll − m_Z\| < 15, p_T^miss < 125, b-veto | 0 jet / ≥ 1 jet / VBF |
| `mllSideband` (DY CR) | ee, μμ | \|m_ll − m_Z\| > 15 (m_ll 12–200), p_T^miss > 100, b-veto | 0 jet / ≥ 1 jet / VBF |
| `emu` | eμ | \|m_eμ − m_Z\| < 15, p_T^miss > 125, b-veto (= N_eμ^in of the AN α method) | 0 jet / ≥ 1 jet / VBF |
| `alphaSB` | ee, μμ, eμ | 40 < m_ll < 70 or 110 < m_ll < 200, p_T^miss > 70, **≥ 1 loose b-tag** (AN 5.1) | inclusive |

**Jet categories** use jets with p_T > 30:
- 0 jet: `ifVBF==0 && nJets==0`;
- ≥ 1 jet: `ifVBF==0 && nJets>=1`;
- VBF: `ifVBF==1`, using the AN VBF tag below.

`yield_table.py` also computes α_ll = N_ll / N_eμ in `alphaSB` from data. AN-2016/325 found
α_ee = 0.369 ± 0.006 and α_μμ = 0.683 ± 0.009.

## Skim changes (nanoAOD_skim, `HZZ2l2nu` block of `config/Input_2017.yml`)

| AN-2016/325 | before | now |
|---|---|---|
| VBF tag (4.3): two **leading** jets, \|Δη\| > 4, m_jj > 500, Z between the jets, central-jet veto | any jet pair with \|Δη\| > 4 and m_jj > 500, highest m_jj, no centrality or CJV | AN definition (`VBF:` block) |
| no charge requirement on ee/μμ (4.3) | Z1 only from OS pairs | OS and SS pairs; Z1 closest to m_Z (`RequireOS: false`); `HZZ2l2qNu_cutOppositeChargeFlag` still stored |
| third-lepton veto incl. soft muons with p_T > 3 GeV (4.4.4) | loose e/μ only | the skim keeps the loose e/μ veto as a cut; soft-ID muons with p_T > 3 are **stored, not cut**: `HZZ2l2nu_nExtraSoftMu` (AN definition) and `HZZ2l2nu_nExtraSoftMuPV` (+ \|dz\| < 0.5, \|dxy\| < 0.2). The plots require `nExtraSoftMu==0`. On 2017 VBF M170 signal the AN definition is 89 % efficient (AN quotes 96–99 % in 2016; the soft ID allows \|dz\| < 20 cm, so pileup muons pass). The PV variant is 95 %. |
| α sidebands up to m_ll = 200 GeV (5.1) | ee/μμ m_ll 12–120 (4l `MZcut`) | 12–200 (`MZLepcut`, `M_ll_Window: 110`) |

## Reproduce

```bash
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim   # no cmsenv needed
mkdir -p plotting/logs
setsid nohup bash plotting/run_all.sh > plotting/logs/run_all.log 2>&1 < /dev/null &
# another skim / output area:  SKIM=<skim dir> BASE=<output dir> bash plotting/run_all.sh
```

`prepare_inputs.py` skips samples whose slim file already exists; pass `--force` (or delete them) after a
re-skim.
