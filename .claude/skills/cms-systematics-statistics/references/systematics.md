# Systematic Uncertainties — X/H&rarr;ZZ&rarr;2&ell;2&nu; (AN-2016/325 &sect;8, Appendix A)

> ## Scope and provenance — read first
>
> Transcribed from **CMS AN-2016/325** &sect;8 (main text), &sect;2.2 (K-factors), and
> **Appendix A** (per-process / per-category / per-mass tables, Tables 24&ndash;51,
> Figures 52&ndash;55). 13 TeV, 35.9 fb&#8315;&sup1;, 2016, CMSSW_8_0_25.
>
> - **Legacy 2016 note.** The *set* of nuisances and *how each enters the fit*
>   are the reference. The *numeric sizes* (2.6% lumi, 15% Top/WW, per-category
>   JES tables, &hellip;) are 2016 and **must be re-derived for the target era**
>   (Run 2 UL / Run 3). Mark any size taken only from here **[Verify]**.
> - Object-level detail (the SF payloads, the JES source split, the lumi value and
>   its correlation scheme) lives in `cms-object-guidelines/references/*` — this
>   file says only **which nuisance, what type, what size, how correlated**.
> - Channel selection / categorization definitions:
>   `cms-object-guidelines/references/hzz-2l2nu.md`.
> - This-repo statistical scaffold: `references/statistics.md` &sect;5 and
>   `higgs_combine/make_datacard.cpp`.

Classification tags: **[Instrumental]**, **[Theory]**, **[Data-driven-bkg]**,
**[AN-2016-325]** (size/method from the note), **[Verify]** (needs a current-era
source), **[Repo divergence]** (datacard differs from the note).

---

## 1. Two families (AN &sect;8)

| Family | Members |
|--------|---------|
| **Instrumental** | luminosity, trigger, lepton ID+iso, lepton momentum scale, JES / JER / unclustered-MET scale, b-jet-veto efficiency, 3rd-(soft-)lepton-veto efficiency, pileup |
| **Theoretical** | signal & ZZ/WZ cross section (scale &mu;_R/&mu;_F + PDF+&alpha;_s), jet-bin categorization (Stewart&ndash;Tackmann), NLO EWK K-factor on q&#773;qZZ / WZ, NNLO QCD K-factor, MC statistics |
| **Data-driven background** | non-resonant (&alpha;-method), Z+jets (&gamma;+jets method) |

The values quoted for "gluon fusion and VBF" in the note are for **m_H = 1 TeV**
unless a per-mass table is cited (Appendix A gives 200&ndash;3000 GeV).

---

## 2. Instrumental uncertainties (AN &sect;8.1)  **[Instrumental]**

| Source | 2016 size | Type | Applies to | Notes |
|--------|-----------|------|-----------|-------|
| **Luminosity** | **2.6%** | `lnN` | all MC (signal + MC bkg) | 2016 value. UL/Run 3: use the LUM recommendation + its cross-year `lumi_1/2/3` scheme — see `cms-object-guidelines/references/lumi.md` |
| **Trigger** | **2%** (di-e), **2%** (di-&mu;) | `lnN` | all MC | di-e efficiency slightly &lt; 100% &rarr; 2%; di-&mu; from &epsilon; vs p_T(Z) / opening-angle study |
| **Lepton ID + isolation** | **1% per lepton** &rarr; **2% per 2e** and **2% per 2&mu;** event | `lnN` | all MC | from the tag-and-probe method spread (&sect;7.2) |
| **Lepton momentum scale** | &mu; **1%**; e **0.6%** (barrel) / **1.5%** (endcap) | `shape` | all MC | shift nominal energy &plusmn;1&sigma;, **propagate to MET**, recompute `M_T` and the category. Per-category electron-scale table: AN Table 18 (e.g. ZZ `=0j` 0.10%, WZ `=0j` 0.33%, ggH `=0j` 0.05%); muon effect negligible |
| **JES** | per-process / per-category, AN Table 19 | `shape` | all MC | vary nominal JES &plusmn;1&sigma; (JetMET prescription); recompute MET, `M_T`, jet category, b-tag. Large in **VBF** (e.g. ZZ VBF `~35%`, WZ VBF `~25%`, ggH VBF `~20%`) |
| **JER** | AN Table 19 | `shape` | all MC | same recompute chain; small vs JES (`~0.1&ndash;2.6%`) |
| **Unclustered-MET scale** | AN Table 19 (`uE_T^miss`) | `shape` | all MC | vary the unclustered component, add the clustered part back, re-assess MET/`M_T` |
| **b-jet-veto efficiency** | **2&ndash;4%** | `lnN`/`shape` | MC-driven only (signal, WZ, ZZ) | envelope: `b`/`c` tag prob. **down 2%**, light-jet tag prob. **up 11%**, &eta;/p_T-dependent, per BTV prescription (&sect;7.3) |
| **3rd-(soft-muon)-lepton veto** | veto efficiency **96&ndash;99%** | `lnN` | MC-driven (ZZ, WZ, signal) | assign the (1 &minus; &epsilon;) spread |
| **Pileup** | from &plusmn;4.6% on &sigma;_minbias (69.2 mb) | `shape` | all MC | `pu` up/down weight; see `cms-object-guidelines/references/pileup.md` |

The "recompute MET &rarr; `M_T` &rarr; category &rarr; b-tag for each &plusmn;1&sigma;"
loop is the reason JES/JER/unclustered-MET/lepton-scale are **shape** nuisances,
not `lnN`: a shift moves events between the 0j / &ge;1j / VBF bins.

---

## 3. Theoretical uncertainties (AN &sect;8.2)  **[Theory]**

### 3.1 MC statistics
Per-process, per-category `lnN` from the finite simulated sample
(AN Table 20 — e.g. ZZ `=0j` 1.7%, WZ `=0j` 8.1%, ggH `=0j` 1.4%, qqH VBF 2.2%).
In modern Combine this is the **autoMCStats / Barlow&ndash;Beeston-lite** bin-by-bin
treatment, not one `lnN`.

### 3.2 PDF
All signal + background samples use **NNPDF**; uncertainty per the **PDF4LHC**
recommendation (`thpdf`). Also propagated to the jet-bin migration.

### 3.3 Jet-bin categorization (ggF 0 / 1 / &ge;2 jets)  **[Theory — Stewart&ndash;Tackmann]**
The scale uncertainties of the **inclusive** 0/1/2-jet cross sections are treated
as **independent** (Stewart&ndash;Tackmann, refs [40,41]); the **exclusive** 0-jet bin
then carries a **correlated** piece. Compute the &mu;_R/&mu;_F variation per event at
gen level (largest of up/down of &mu;_R, &mu;_F, or both together); for the 0-jet bin:

```
Delta sigma^{0j} = sqrt( (Delta sigma^{>=0j})^2 + (Delta sigma^{>=1j})^2 )
```

Per-mass signal scale table: AN Table 21 (e.g. m=200 `=0j` **104%**, m=1000 `=0j`
43%, VBF `~20%` across masses).

### 3.4 ZZ / WZ cross section
`Delta sigma_tot = Delta_PDF + Delta_alpha_s + (mu_R, mu_F)` with the scales
varied by 0.5 and 2 (MadGraph5_aMC@NLO), max relative shift taken; the 0-jet bin
from the same `sqrt(...)` combination as &sect;3.3. Results: AN Table 22
(ZZ `=0j` **6.7%**, `&ge;1j` 6.0%, VBF **41%**; WZ `=0j` 10.1%, `&ge;1j` 6.0%,
VBF 40%) — the `qcdscalezz` / `qcdscalewz` nuisances.

### 3.5 NLO EWK K-factor on q&#773;qZZ and WZ  **[Theory — &rho;-method]**
The EWK correction is negative (`~ -4%` &rarr; `~ -10%` at high `m_ZZ` / `p_T,Z`);
forced to 1 for `m_ZZ < 2 m_Z`. Its uncertainty uses the event recoil

```
rho = | sum_i vec p_T^i | / sum_i | vec p_T^i |     (i = the 4 leptons / their Z's)
```

```
delta = | (1 - K_QCD^NLO) (1 - K_EWK^NLO) |   if rho <  0.3   (missing alpha*alpha_s diagrams; maximal when QCD & EWK pull the same way; K_QCD^NLO = 15.99/9.89 from ref [10])
delta = | 1 - K_EWK^NLO |                     if rho >= 0.3   (100% of the EWK correction; ~25% of events)
```

Carried as the `thewk` shape nuisance (AN Fig. 43 shows the up/down `M_T`
templates). WZ EWK: virtual part as for ZZ; photon-induced part &lt; 1% (LUXqed
PDF) &rarr; neglected. NNLO-QCD-K-factor uncertainty: "still to be implemented"
in the note (&sect;8.2) &rarr; **[Verify]** for the current analysis.

---

## 4. Data-driven background uncertainties (AN &sect;8.1)  **[Data-driven-bkg]**

| Background | Nuisance(s) | Size / method |
|-----------|-------------|---------------|
| **Non-resonant** (Top/WW, via the e&mu; &alpha;-method, `references/statistics.md` &sect;4 / `hzz-2l2nu.md` &sect;5) | `topwww` (`lnN`, note uses **13%** in Appendix A / **15%** in &sect;8.1 text), plus `stat` on the sideband counts | total driven by the stat. uncertainty on `&alpha;_e`, `&alpha;_&mu;` + the MC-closure bias (&le; 13% at MET = 125). Shape (`M_T`) taken from e&mu; data. |
| **Z+jets** (instrumental MET, &gamma;+jets method) | `zllinstrmet` (**shape**) + `normalizationzllinstrmet` (**very large `lnN`**, Appendix A Table 25: `+130%&hellip;+300%`) + `stat` | three physical components: (a) photon-`p_T`&rarr;Z-`p_T` reweighting closure `< 10%`; (b) genuine-MET-subtraction (vary W&gamma;/Z&gamma;/W+jets &plusmn;1&sigma;, incl. their EWK corrections, quad-sum bin-by-bin); (c) photon-sample statistics. AN Fig. 42. |

---

## 5. Combine nuisance names (AN Appendix A, Tables 24&ndash;51)  **[AN-2016-325]**

The exact nuisance strings the 2016 datacard used, per process. Use these as the
*naming reference*; re-derive the sizes.

| Nuisance | Type | Meaning | Applied to |
|----------|------|---------|-----------|
| `lumi` | `lnN` | integrated luminosity (2.6%) | all MC |
| `trigger` (`eff` in text) | `lnN` | trigger efficiency 2% / 2% | all MC |
| `effe` / `effm` | `lnN` | electron / muon ID+iso, 2% per 2e / 2&mu; | all MC |
| `e` | `shape` | electron energy scale | all MC |
| `scalem` | `shape` | muon momentum scale (tiny) | all MC |
| `scalej` | `shape` | jet energy **scale** (JES) | all MC |
| `resj` | `shape` | jet energy **resolution** (JER) | all MC |
| `scaleumet` | `shape` | unclustered-MET scale | all MC |
| `resrhoe` | `shape` | &rho; / energy-density term in the electron isolation/energy — **meaning inferred**: AN Appendix A lists the name but never defines it; **[Verify]** | all MC |
| `effb` | `shape` | b-jet-veto efficiency (2&ndash;4%) | signal, WZ, ZZ, ZVV |
| `lepveto` | `shape` | 3rd-(soft-muon)-lepton veto efficiency | signal, WZ, ZZ, ZVV |
| `pu` | `shape` | pileup reweighting | all MC |
| `stat` | `shape` | MC statistics of that process/category (&rarr; use autoMCStats now) | all MC + data-driven |
| `qcdscaleggh` | `lnN` | ggF signal &mu;_R/&mu;_F + jet-bin (Table 21) | signal ggH |
| `qcdscalezz` / `qcdscalewz` | `lnN` | ZZ / WZ scale+PDF+&alpha;_s (Table 22) | ZZ / WZ |
| `thpdf` | `shape` | PDF (PDF4LHC / NNPDF) | all MC |
| `thalphas` | `shape` | &alpha;_s | all MC |
| `thewk` | `shape` | NLO EWK K-factor, &rho;-method (&sect;3.5) | ZZ, WZ (q&#773;q-induced) |
| `gse` | `shape` | GEN-level scale / weight envelope (large only at m &ge; 3 TeV) | signal |
| `normalizationzllinstrmet` | `lnN` | Z+jets (Instr. MET) normalization | Instr. MET |
| `zllinstrmet` | `shape` | Z+jets (Instr. MET) shape | Instr. MET |
| `topwww` | `lnN` | non-resonant (Top/W/WW) &alpha;-method, 13&ndash;15% | Top/W/WW |

Shape summary plots: AN Figs. 52 (irreducible ZZ/WZ/ZVV), 53 (Instr. MET +
non-resonant), 54 (ggH signal), 55 (qqH signal).

---

## 6. Cross-check vs this repo (`higgs_combine/make_datacard.cpp`)  **[Repo divergence]**

The current datacard generator is an **early scaffold**, not the AN model.

> **The three macros are not wired to each other.** `make_shapes.cpp` writes
> `input_root_file_combine_v4.root`; `make_datacard.cpp` reads
> `input_root_file_combine_v3.root`; `make_shapes_final.cpp` is a **third**,
> larger variant (per-W-jet bins, a `category` struct). Before reviewing any
> datacard, establish which shape file is current and regenerate it — a stale
> `_v3` is almost certainly what `make_datacard.cpp` picks up.

| AN-2016/325 &sect;8 | `higgs_combine/make_datacard.cpp` |
|---------------------|-----------------------------------|
| ~20 nuisances, mostly `shape`, per-process / per-category / per-mass | **two** lines only: `lumi lnN` and `alpha shape` |
| `lumi lnN` = **1.026** (2.6%), all MC | `lumi lnN` = **0.84** for every process — wrong value **and** an `lnN` &lt; 1 is a ~16% *down-only* kink, almost certainly a bug |
| shape systematics from real &plusmn;1&sigma; recomputations | `make_shapes.cpp` writes `_alphaUp` / `_alphaDown` as **exact clones** of nominal &rarr; the `alpha` nuisance is currently null |
| 6 categories (`ee`/`&mu;&mu;` &times; `=0j`/`&ge;1j`/`VBF`) combined; ggF vs VBF limits | **one** bin `bin1`, no category split |
| signal = high-mass X (200&ndash;3000 GeV) &times; &Gamma; (5/10/100) &times; interference (Eq. 14) | signal = single `signal_alpha` from `signal_m125.root` (SM 125 GeV) |
| lumi 35.9 fb&#8315;&sup1; (2016) | `make_shapes.cpp` `lumi = 59547 pb&#8315;&sup1;` (2018) |
| non-resonant from e&mu; **data** (&alpha;-method) | Top/single-top/tt+X taken straight from **MC** |
| Z+jets from &gamma;+jets **data** | DY taken straight from **MC** (`dy_alpha`), no Instr.-MET nuisance |
| `M_T` per AN Eq. 6 | branch `HZZ2l2nu_ZZmT`, 100 bins 0&ndash;1000 GeV (definition — see `hzz-2l2nu.md` &sect;3) |

Treat each row as an **analysis-specific inconsistency / missing systematic** to
raise with the analyst — the scaffold is expected to grow, but none of this is
documented.

---

## 7. Review checklist

1. Every AN &sect;8 nuisance present in the datacard, or its absence justified.
2. `lumi` is a small `>1` `lnN` for the era (not `0.84`); correlated across
   categories; cross-year scheme if multiple years.
3. Instrumental shape nuisances (`scalej`, `resj`, `scaleumet`, `e`, `pu`) are
   built from **real** &plusmn;1&sigma; template recomputations, not clones, and the
   recompute chain re-bins MET &rarr; `M_T` &rarr; category &rarr; b-tag.
4. `effb`, `lepveto` applied to MC-driven processes only (signal, WZ, ZZ, ZVV);
   not to the data-driven Top/WW or Z+jets.
5. `qcdscale*` + `thpdf` + `thalphas` + jet-bin Stewart&ndash;Tackmann combination
   (&sect;3.3) present for signal and ZZ/WZ.
6. `thewk` uses the &rho;-method (&sect;3.5), split by `&rho; &lessgtr; 0.3`.
7. Non-resonant = `topwww` (13&ndash;15%) on an e&mu;-**data** template; Z+jets =
   `zllinstrmet` (+ its normalization) on a &gamma;+jets-**data** template — not MC.
8. MC-stat handled by autoMCStats, not a per-process `lnN`.
9. Correlation: object SF/scale nuisances correlated across categories & (usually)
   years; MC-stat uncorrelated; `topwww` / Z+jets per their control regions.
10. Every size that traces only to AN-2016/325 flagged **[Verify]** for the era.

---

## 8. Evidence summary

| Item | Source | Established? |
|------|--------|--------------|
| Instrumental nuisance list + 2016 sizes | AN &sect;8.1, Tables 18&ndash;19 | yes — 2016; **[Verify]** current era |
| Theory nuisance list, jet-bin S&amp;T, ZZ/WZ xsec | AN &sect;8.2, Tables 20&ndash;22 | yes — 2016 |
| EWK-K-factor &rho;-method | AN &sect;2.2, &sect;8.2, refs [10,12,13] | yes — method |
| NNLO-QCD-K-factor uncertainty | AN &sect;8.2 | **not implemented in the note** — [Verify] |
| Combine nuisance names | AN Appendix A, Tables 24&ndash;51 | yes — 2016 datacard naming |
| Data-driven bkg uncertainties (13&ndash;15% / 3-component) | AN &sect;5.1, &sect;5.2, &sect;8.1 | yes — method |
| Repo datacard vs AN (&sect;6) | `higgs_combine/make_datacard.cpp` + `make_shapes.cpp` | scaffold only — divergences listed |

## Last verified

- AN-2016/325 &sect;8 + Appendix A transcription: skill creation.
- `higgs_combine/` cross-check: against the repo macros at skill-creation time
  (ROOT not run — macros read for structure only).
- Current UL / Run 2+3 2&ell;2&nu; systematics model: **not consulted — [Verify]**.
