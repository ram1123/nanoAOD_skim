# Data and MC samples

All samples are **Run 2 UltraLegacy, NanoAOD v15**:

- MC: `RunIISummer20UL1{6APV,6,7,8}NanoAODv15`, global tag `150X_mc20XX_realistic_v1`.
- Data: `Run20XX*-{HIPM_,}UL20XX_NanoAODv15`.

The sample lists are in `input_data_files/` and use one DAS dataset per line, with
`#` in column 0 for a comment. Both the Slurm and the HTCondor submission scripts
read this format.

| list | content |
| --- | --- |
| `sample_list_AN2016_325_v15_<year>.dat` | **main list**: MC + data, one Slurm submission per year |
| `sample_list_AN2016_325_v15_<year>_data.dat` | the data part only |
| `sample_list_gammaCR_v15_2017.dat` | &gamma;+jets CR (run with `--channels 1g`) |
| `sample_list_v15_<year>[_data].dat` | older, shorter v15 lists |
| `sample_list_AN2016_325_v9_*` | NanoAOD v9 equivalents (not usable on this branch) |

The `AN2016_325` lists are **generated** by `input_data_files/make_AN2016_325_lists.py`,
and re-running it overwrites hand edits. The script maps every process of
AN-2016/325 (Tables 1, 2, 3 and 52) to its UL equivalent and resolves the exact
DAS names with `dasgoclient`. It needs a VOMS proxy. For each primary dataset it
picks:

- the standard NanoAOD processing (JMENano / BTVNano only if nothing else exists);
- the highest `-vN`;
- **every** `_extN` extension. Extensions are independent statistics, so
  normalise with the summed `genEventSumw`.

A process with no DAS match is written as `# NOT IN DAS: <pattern>`.

## Data

| year | primary datasets | eras | integrated luminosity |
| --- | --- | --- | --- |
| 2016preVFP | DoubleMuon, DoubleEG, MuonEG, SingleMuon, SingleElectron, SinglePhoton | B (v1 + v2), C, D, E, F (all `HIPM_UL2016`) | 19.5 fb<sup>&minus;1</sup> |
| 2016postVFP | same six | F, G, H (`UL2016`) | 16.8 fb<sup>&minus;1</sup> |
| 2017 | same six | B&ndash;F | 41.48 fb<sup>&minus;1</sup> |
| 2018 | DoubleMuon, MuonEG, SingleMuon, **EGamma** | A&ndash;D | 59.83 fb<sup>&minus;1</sup> |

These are the standard UL golden-JSON luminosities. The code holds only two of
them: `LUMI_PB = 41480` in `plotting/make_config.py` and `lumi: 59830.` in
`Input_2018.yml`.

Because one event can be in several primary datasets, data must be
**de-duplicated** after skimming. The first occurrence of each
(run, luminosityBlock, event) is kept, in this priority order:
DoubleMuon &rarr; DoubleEG &rarr; MuonEG &rarr; SingleMuon &rarr; SingleElectron
&rarr; SinglePhoton. `prepare_inputs.py` does this for the plots (`DATA_PDS`), and
`scripts/analysis/removeDuplicates.C` for other uses (inputs listed in priority
order).

!!! warning "2018 `EGamma`"
    `DATA_PDS` in `prepare_inputs.py` does not include `EGamma`. A 2018 skim
    directory `EGamma/` would therefore be treated as an MC sample and left out of
    the merged data file. Add it to `DATA_PDS`, after DoubleMuon, before preparing
    2018 inputs.

## MC: dilepton analysis (`sample_list_AN2016_325_v15_<year>.dat`)

Dataset counts include `_ext` extensions.

| group (AN table) | samples | 2016pre | 2016post | 2017 | 2018 |
| --- | --- | :-: | :-: | :-: | :-: |
| **Signal** (AN &sect;3) | `GluGluHToZZTo2L2Nu_M*` and `VBF_HToZZTo2L2Nu_M*`, POWHEG + JHUGen. M = 125, 160&ndash;3000 GeV | 50 | 34 | 50 | 50 |
| Z&rarr;&ell;&ell; Drell-Yan | `DYJetsToLL_M-10to50`, `DYJetsToLL_M-50` (amc@NLO FxFx), `DYJetsToLL_LHEFilterPtZ-*` (6 p<sub>T</sub><sup>Z</sup> bins) | 8 | 8 | 8 | 8 |
| t&#773;t + X | `TTTo2L2Nu` (POWHEG), `TTWJetsToLNu`, `TTZToLLNuNu_M-10` | 3 | 3 | 3 | 3 |
| Single top | tW (inclusive and NoFullyHadronic), t-channel, s-channel lepton/hadron decays, `tZq_ll` | 9 | 9 | 9 | 9 |
| Dibosons | `ZZTo2L2Nu`, `ZZTo4L`, `ZZTo2Q2L`, `ZZTo2Nu2Q`, `WZTo3LNu`, `WZTo2Q2L`, `WZTo1L1Nu2Q`, `WZTo1L3Nu`, `WZTo2Q2Nu`, `WWTo2L2Nu`, `WWTo1L1Nu2Q`, gg&rarr;ZZ&rarr;2&ell;2&nu; (MCFM, 2e2&nu; / 2&mu;2&nu;), VBF&rarr;ZZ&rarr;2&ell;2&nu; continuum, gg&rarr;WW (MCFM) | 17 | 17 | 17 | 17 |
| Tribosons | ZZZ, WZZ, WWZ, WWW | 8 | 8 | 8 | 8 |
| W&rarr;&ell;&nu; | inclusive amc@NLO plus HT-binned MLM (70&ndash;&infin;) | 22 | 17 | 22 | 21 |
| Z&rarr;&nu;&nu; | `ZJetsToNuNu_HT-*` (100&ndash;&infin;), for the instrumental-MET closure | 7 | 7 | 7 | 7 |
| &gamma;+jets | `GJets_HT-*` (40&ndash;&infin;) | 5 | 5 | 5 | 5 |
| Top + &gamma; | `TTGJets`, `TGJets` | 2 | 2 | 2 | 3 |
| Z&gamma; | `ZNuNuGJets_MonoPhoton_PtG-40to130 / -130`, `ZGToLLG_01J_5f` | 3 | 3 | 3 | 3 |
| W&gamma; | `WGToLNuG` (MLM) | 1 | 1 | 1 | 1 |
| QCD HT | `QCD_HT*` (100&ndash;&infin;) | 16 | 16 | 8 | 8 |
| QCD EM-enriched | `QCD_Pt-*_EMEnriched` (15&ndash;&infin;) | 8 | 8 | 8 | 8 |
| QCD &mu;-enriched | `QCD_Pt-*_MuEnrichedPt5` (15&ndash;1000+) | 12 | 12 | 12 | 12 |

!!! warning "Gaps in 2016postVFP signal"
    Some signal mass points are not in DAS at v15 for 2016postVFP:

    - ggH: M = 125, 160&ndash;190, 300, 450, 500, 800, 1500, 2500
    - VBF: M = 190, 250, 300, 450, 2500

    2016 QCD-HT has both `madgraph`- and `madgraphMLM`-named samples, which is why
    it has 16 datasets instead of 8.

### What the 2017 control plots actually use

Not every skimmed sample goes into the stack. `plotting/make_config.py` uses:

| stack group | samples | note |
| --- | --- | --- |
| Z+jets | `DYJetsToLL_M-10to50`, `DYJetsToLL_M-50` | inclusive amc@NLO |
| t&#773;t | `TTTo2L2Nu` | |
| Diboson | WW&rarr;2&ell;2&nu;, WW&rarr;&ell;&nu;qq, WZ&rarr;3&ell;&nu;, WZ&rarr;2q2&ell;, ZZ&rarr;2&ell;2&nu;, ZZ&rarr;2q2&ell;, gg&rarr;ZZ&rarr;2e2&nu; / 2&mu;2&nu; | |
| Single top | s-, t-channel (t, t&#773;), tW (t, t&#773;) | |
| W+jets | HT-binned MLM, 100&ndash;&infin; | the inclusive NLO sample has a few events with weights of &plusmn;145 |
| TTV | TTW, TTZ, TT&gamma;, t&gamma; | |
| Triboson | WWZ, WZZ, ZZZ | |

The plots leave out QCD, GJets, Z&rarr;&nu;&nu; and Z&nu;&nu;&gamma;, which have no
genuine dilepton pair. They also leave out Z&gamma; and W&gamma;, which overlap
with DY and W+jets because no photon overlap removal is applied, and the signal
samples. Cross sections are in `make_config.py`. Each MC sample is normalised as

\[
w_\text{sample} = \frac{\sigma \cdot L}{\sum_\text{all events} \text{genEventSumw}},
\]

and every event carries `overallEventWeight` (see [Corrections](corrections.md#event-weights)).

## MC: &gamma;+jets control region (`sample_list_gammaCR_v15_2017.dat`)

2017 only, because the `GammaCR` YAML block exists only in `Input_2017.yml`.

| role | samples |
| --- | --- |
| data | `SinglePhoton` Run2017B&ndash;F |
| &gamma;+jets (target) | `GJets_HT-*` 40&ndash;&infin; |
| fake photons | `QCD_HT*` 100&ndash;&infin;, `QCD_Pt-*_EMEnriched` |
| genuine MET + &gamma; | `ZNuNuGJets_MonoPhoton` (2 bins), `WGToLNuG`, `TTGJets`, `TGJets` |
| other | `ZGToLLG_01J_5f`, `WJetsToLNu_HT-*` 100&ndash;&infin; |

The plots use QCD-HT only, because the EM-enriched p<sub>T</sub> bins overlap it.
The &gamma;-CR cross sections in `make_config.py` are marked *please check*.
