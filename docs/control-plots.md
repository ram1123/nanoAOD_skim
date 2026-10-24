# Data/MC control plots

The `plotting/` directory makes stacked data/MC control plots and yield tables
from the Slurm skims, with the **AN-2016/325** event selection on top of UL / v15
object IDs. It is set up for **2017** (41.5 fb<sup>&minus;1</sup>) and runs on the
AF host with the shared pixi ROOT. It needs **no CMSSW and no container**.

## Pipeline

| step | script | output |
| --- | --- | --- |
| 0. skim | Slurm, see [Running](running.md#batch-slurm-on-purdue-hammer) | `<SKIM>/<sample>/*_Skim.root` |
| 1. slim + merge per sample, &Sigma; genEventSumw, data de-duplication | `scripts/analysis/prepare_inputs.py` | `$BASE/inputs/<sample>.root`, `data_2017.root`, `sumw.json` |
| 2. sample list, cut files, variable header | `plotting/make_config.py` | `HZZ2l2nu_Samples13TeV_2017_v15.txt`, `cuts/*.txt`, `controlplotvars_2017.h` (generated, git-ignored) |
| 3. plots per region | `plotting/RunMacro_2017_v15.C` &rarr; `myControlPlots.C` | `$BASE/plots/<region>/<tag>/` |
| 4. yields and &alpha; | `plotting/yield_table.py` | `plotting/yields_2017.md` |
| 5. provenance | `plotting/write_metadata.sh` | `$BASE/{inputs,plots}/_metadata/` |

`plotting/run_all.sh` runs steps 2&ndash;5.

!!! warning "Step 1 is commented out in `run_all.sh`"
    The `prepare_inputs.py` line in `run_all.sh` is commented out, so after a new
    skim you must run step 1 yourself, as shown below. Otherwise the plots are made
    from the old slim inputs.

## How to make the dilepton plots

From `nanoAOD_skim/`, in a **fresh shell**. Do not run `cmsenv` or `set_env.sh`:
a CMSSW `LD_LIBRARY_PATH` makes pixi ROOT abort with a cling `--optimize-regalloc`
error.

```bash
export SKIM=/depot/cms/users/$USER/HZZ2l2nu_skim/skims/<submission_name>   # input: the Slurm skim
export BASE=/depot/cms/users/$USER/HZZ2l2nu_skim/controlplots_<year>_<tag> # output area
PY=/work/pixi/global/.pixi/envs/default/bin/python

# 1. slim inputs (skips samples that already exist; add --force after a re-skim)
$PY scripts/analysis/prepare_inputs.py --skimdir $SKIM --outdir $BASE/inputs --threads 8

# 2-5. config, plots for all four regions in parallel, yields, metadata
mkdir -p plotting/logs
setsid nohup bash plotting/run_all.sh > plotting/logs/run_all.log 2>&1 < /dev/null &
```

Set `SKIM` and `BASE` every time. If they are unset, `run_all.sh` silently falls
back to the last 2017 run (`skims/HZZ2l2nu_AN325_2017_MCdata`,
`controlplots_2017_AN325`).

To redo a single region after editing `make_config.py`:

```bash
cd plotting
$PY make_config.py --inputs $BASE/inputs
/work/pixi/global/.pixi/envs/default/bin/root -l -b -q "RunMacro_2017_v15.C(\"LowMET\", \"$BASE/plots\")"
```

Run it with the global pixi ROOT on `PATH` and a clean `LD_LIBRARY_PATH`, as
`run_all.sh` sets up.

## How to make the &gamma;+jets CR plots

These use a separate skim (`--channels 1g`, list `sample_list_gammaCR_v15_2017.dat`)
and are **not** part of `run_all.sh`:

```bash
$PY scripts/analysis/prepare_inputs.py --mode gamma --skimdir $SKIM_GAMMA \
    --outdir $BASE_GAMMA/inputs_gamma --threads 8
cd plotting
$PY make_config.py --inputs $BASE/inputs --inputs_gamma $BASE_GAMMA/inputs_gamma
root -l -b -q "RunMacro_2017_v15.C(\"gammaCR\", \"$BASE_GAMMA/plots\")"
```

`--mode gamma` keeps `passGammaCR_Selection && gamma_trigWindowFired` events and
takes data from `SinglePhoton` only. The macro weights data by
`gamma_prescaleWeight` (`gDataWeight`).

## Normalisation

- **Data:** the de-duplicated merge of all primary datasets, weight 1 (or the
  prescale weight in the &gamma; CR).
- **MC:** each event is weighted by `overallEventWeight` (`gMCWeight`), and each
  sample is scaled by \(\sigma \cdot L / \sum \text{genEventSumw}\). `make_config.py`
  writes that factor as `otherscale = 1/sumw`, with `nMCevents = 1` and
  `nMCNegEvents = 0`, in the 8-column sample list that `myControlPlots.C` reads.
  L = 41 480 pb<sup>&minus;1</sup>.

Cross sections, stack groups and colours are edited **only** in the `MC` and
`GAMMA_MC` lists in `make_config.py` (see [Samples](samples.md#what-the-2017-control-plots-actually-use)).

## Control regions { #control-regions }

**Pre-selection in every dilepton region** (AN &sect;4.4), on top of the skim:

| cut | branch expression |
| --- | --- |
| region flag + trigger | `passZZ2l2nuSelection==1 && Triggers_HZZ2l2nu_SR==1`, or `passZZ2l2nu_emuCR_Selection==1 && Triggers_HZZ2l2nu_EMuCR==1` |
| lepton p<sub>T</sub> | `pTL1>25 && pTL2>25` |
| p<sub>T</sub><sup>&ell;&ell;</sup> > 55 | `pTZ1>55` |
| no extra soft muon | `HZZ2l2nu_nExtraSoftMu==0` |
| min &Delta;&phi;(jet, MET) > 0.5 | `HZZ2l2nu_minDPhi_METAK4>0.5` |
| \|&Delta;&phi;(Z, MET)\| > 0.5 | `HZZ2l2nu_dPhi_ZMET>0.5` |
| no charge requirement | (AN &sect;4.3) |

| region | flavours | extra cuts | categories |
| --- | --- | --- | --- |
| `LowMET` | ee, &mu;&mu; | \|m<sub>&ell;&ell;</sub> &minus; m<sub>Z</sub>\| < 15, p<sub>T</sub><sup>miss</sup> < 125, b-veto | 0 jet / &ge; 1 jet / VBF |
| `mllSideband` (DY CR) | ee, &mu;&mu; | \|m<sub>&ell;&ell;</sub> &minus; m<sub>Z</sub>\| > 15, p<sub>T</sub><sup>miss</sup> > 100, b-veto | 0 jet / &ge; 1 jet / VBF |
| `emu` | e&mu; | \|m<sub>e&mu;</sub> &minus; m<sub>Z</sub>\| < 15, p<sub>T</sub><sup>miss</sup> > 125, b-veto (= N<sub>e&mu;</sub><sup>in</sup> of the &alpha; method) | 0 jet / &ge; 1 jet / VBF |
| `alphaSB` | ee, &mu;&mu;, e&mu; | 40 < m<sub>&ell;&ell;</sub> < 70 or 110 < m<sub>&ell;&ell;</sub> < 200, p<sub>T</sub><sup>miss</sup> > 70, **&ge; 1 loose b-tag** (AN &sect;5.1) | inclusive |
| `gammaCR` | &gamma; | `passGammaCR_Selection`, `gamma_trigWindowFired`, p<sub>T</sub><sup>&gamma;</sup> > 55, min &Delta;&phi;(jet, MET) > 0.5, \|&Delta;&phi;(&gamma;, MET)\| > 0.5, b-veto | 0 jet / &ge; 1 jet / VBF |

- **b-veto:** `HZZ2l2qNu_nLooseBtagJets==0`, the UParTAK4 Loose WP.
- **Jet categories:**
    - 0 jet: `ifVBF==0 && nJets==0`;
    - &ge; 1 jet: `ifVBF==0 && nJets>=1`;
    - VBF: `HZZ2l2nu_ifVBF==1`.
- **Blinding:** the ee/&mu;&mu; signal region (on-Z, p<sub>T</sub><sup>miss</sup> > 125)
  is not in any region.

`yield_table.py` also computes &alpha;<sub>&ell;&ell;</sub> = N<sub>&ell;&ell;</sub> / N<sub>e&mu;</sub>
in `alphaSB` from data. For comparison, AN-2016/325 found
&alpha;<sub>ee</sub> = 0.369 &plusmn; 0.006 and &alpha;<sub>&mu;&mu;</sub> = 0.683 &plusmn; 0.009.

## Plotted variables

- **All categories:**
    - m<sub>&ell;&ell;</sub>, p<sub>T</sub><sup>&ell;&ell;</sup>, &eta;<sup>&ell;&ell;</sup>, &phi;<sup>&ell;&ell;</sup>;
    - lepton p<sub>T</sub> and &eta;, &Delta;R(&ell;&#8321;, &ell;&#8322;);
    - p<sub>T</sub><sup>miss</sup>, &phi;<sup>miss</sup>, m<sub>T</sub><sup>ZZ</sup>, &Delta;&phi;(Z, MET);
    - number of good PVs.
- **Categories with jets:** min &Delta;&phi;(jet, MET), N<sub>jets</sub>,
  N<sub>loose b</sub>, leading-jet p<sub>T</sub> and &eta;.
- **VBF only:** m<sub>jj</sub>, \|&Delta;&eta;<sub>jj</sub>\|, sub-leading jet &eta;.

Ranges and binning are set per region in `VARS`, `JETVARS` and `VBFVARS` in
`make_config.py`, which generates `controlplotvars_2017.h`.

## Where the 2017 outputs are

| what | path |
| --- | --- |
| skims | `/depot/cms/users/shar1172/HZZ2l2nu_skim/skims/HZZ2l2nu_AN325_2017_MCdata/<sample>/` |
| slim inputs | `/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325/inputs/` |
| plots | `/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325/plots/<region>/<tag>/` |
| yields | `plotting/yields_2017.md` |

## Porting to another year

1. Add `EGamma` to `DATA_PDS` in `prepare_inputs.py` (2018) and pass
   `--data_file data_<year>.root`.
2. In `make_config.py`, change `LUMI_PB` and the data file name, and check the
   cross sections and sample names in `MC`.
3. In `RunMacro_2017_v15.C`, change the `lumi_13TeV` label.
4. Check that the year's `Triggers_HZZ2l2nu_SR` / `_EMuCR` branches exist in the
   skim.
