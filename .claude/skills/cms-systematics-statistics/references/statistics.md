# Statistical Analysis — X/H&rarr;ZZ&rarr;2&ell;2&nu; (AN-2016/325 &sect;3, &sect;6, &sect;9)

> ## Scope and provenance — read first
>
> Transcribed from **CMS AN-2016/325**: &sect;3 (signal simulation + MELA
> reweighting), &sect;6 (optimization + shape-based search), &sect;9 (final yields,
> interference model Eq. 14, upper limits), Figs. 44&ndash;51, Table 23. 13 TeV,
> 35.9 fb&#8315;&sup1;, 2016.
>
> - **Legacy 2016 note.** The *method* (shape fit on `M_T`, MELA reweighting for
>   width/interference, CL_s limits per category) is the reference. Mass grid,
>   yields, and limit numbers are 2016 — **[Verify]** for the current analysis.
> - Nuisance parameters and their sizes: `references/systematics.md`.
> - `M_T` definition, selection, categorization: `cms-object-guidelines/references/hzz-2l2nu.md`.
> - The repo's implementation is `higgs_combine/` (standalone ROOT macros) — an
>   early scaffold; &sect;5 below maps it against this note.

Classification tags: **[AN-2016-325]**, **[Analysis-specific]**, **[Repo
divergence]**, **[Verify]**.

---

## 1. The fit observable — transverse mass `M_T` (AN &sect;4.3 Eq. 6, &sect;6.2)

- **Shape-based** maximum-likelihood fit to the **`M_T`** distribution — the
  reconstructed transverse mass of the ZZ&rarr;2&ell;2&nu; system:

  ```
  M_T^2 = [ sqrt(pT_ll^2 + m_ll^2) + sqrt(MET^2 + m_Z^2) ]^2
          - [ vec(pT_ll) + vec(MET) ]^2
  ```

  (Z&rarr;&nu;&nu; leg assigned the PDG Z mass.)
- **No optimization on `M_T` itself** (it *is* the fit variable). The one cut that
  *was* optimized on top of the selection is **`MET > 125 GeV`** — scanned
  50&ndash;140 GeV in 5 GeV steps per mass / width / category; significance is flat,
  so a single common `MET &ge; 125 GeV` is used everywhere (AN &sect;6.1).
- Template shapes:
  - **signal, WZ, ZZ** &rarr; MC simulation;
  - **Z+jets (Instr. MET)** &rarr; &gamma;+jets **data** (`references/systematics.md` &sect;4);
  - **non-resonant (Top/WW/&hellip;)** &rarr; e&mu; **data**, both **shape and yield**
    predicted from the e&mu; control region via the &alpha;-method.
- Fit range / binning in the note's plots: `M_T` up to ~1000&ndash;3000 GeV
  depending on the mass point (AN Figs. 45&ndash;48). Repo `make_shapes.cpp` uses
  **100 bins, 0&ndash;1000 GeV**, branch `HZZ2l2nu_ZZmT`.

---

## 2. Signal simulation and MELA reweighting (AN &sect;3)  **[Analysis-specific]**

- Production: **POWHEG** (ggF + VBF), Higgs decay ZZ&rarr;2&ell;2&nu; by **JHUGen**;
  NLO QCD in production, **no** interference in the raw sample.
- Mass grid **200&ndash;3000 GeV** (100 GeV steps to 1 TeV, 500 GeV steps to 3 TeV);
  widths **&Gamma; = 5, 10, 100 GeV**.
- Each POWHEG sample is generated with **&sigma; = 1 pb** and then **rescaled by a
  matrix-element (MELA) reweighting** that installs the correct width **and** the
  three interference terms:
  - **Int1** = H &times; continuum `gg&rarr;ZZ&rarr;2&ell;2&nu;`;
  - **Int2** = SM h(125) &times; H;
  - **Int3** = h(125) &times; continuum.
  Int1+Int2 raise the low-`M_T` tail and deplete the high side; Int3 removes the
  high-side depletion &rarr; net enhancement.
- Normalization: the reweighted continuum `M_T` is rescaled to the **real
  continuum** (LO &times; NNLO K-factors); the resulting scale factor is applied to
  **all three** shapes (pure resonance / "Background" / total), then NNLO
  K-factors on top. See `references/systematics.md` &sect;3 for the K-factor
  uncertainties.
- MELA package: JHUGenMELA (this repo pins `external/JHUGenMELA` at `v2.4.2`; the
  4l `D_*` discriminants are documented in
  `cms-object-guidelines/references/mela-discriminants.md` — **not** used in
  2&ell;2&nu;, which fits `M_T` directly).

---

## 3. Signal-strength &amp; interference model (AN &sect;9.2 Eq. 14)  **[AN-2016-325]**

Because the signal interferes with the SM continuum and the 125 GeV Higgs, the
per-mass cross section is parameterized in the signal strength &mu; as:

```
sigma_gg->(H->)VV(mu, m_VV) =
    (mu - sqrt(mu)) * sigma_{gg->H->VV}(m_VV)      # pure heavy resonance
  +  sqrt(mu)       * sigma_{gg->(H->)VV}(m_VV)     # "Background" = continuum + h(125)
  + (mu - sqrt(mu)) * sigma_{gg->VV}(m_VV)          # continuum-only interference piece
```

so that at &mu; = 1 the sum is the full SM+H prediction and at &mu; = 0 it reduces to
the background. MELA produces each of the three &sigma; templates for a fixed mass
with the correct relative scaling; the datacard multiplies them by the
&mu;-dependent coefficients above (a physics model / `PhysicsModel` in Combine
terms). **[Verify]** whether the current analysis keeps this exact
parameterization.

- Limits are set **separately for ggF and VBF** production (no assumption on their
  ratio), using the `M_T` shape per jet-bin category (AN &sect;9.2, Fig. 48).
- Widths tested: &Gamma; = 5 / 10 / 100 GeV (Fig. 49).

---

## 4. Non-resonant prediction — the &alpha; (different-flavour) method (AN &sect;5.1)

Feeds both the yield and the `M_T` shape of the Top/WW/W+jets/Z&rarr;&tau;&tau;
background in the signal region:

```
N^NR_ll   = alpha_l * N^in_emu
alpha_l   = N^out_ll / N^out_emu       (from m_ll sidebands 40-70 & 110-200 GeV,
                                        MET > 70 GeV, >= 1 b-tag)
```

2016 values `alpha_ee ~ 0.37 +/- 0.01`, `alpha_mumu ~ 0.68 +/- 0.01`;
cross-checked with the **k-method** (`N^peak_ll = k_ll * N^peak_emu`,
`k_ll = 0.5 * sqrt(N^peak_ll / N^peak_ll')`). The `M_T` shape is taken from the
e&mu; channel in **data** (MC for VBF where e&mu; stats are too low). Uncertainty:
`references/systematics.md` &sect;4 (`topwww`, 13&ndash;15%).

---

## 5. This repo's `higgs_combine/` — state vs the AN  **[Repo divergence]**

Standalone ROOT macros (CLAUDE.md: `root -l -b -q higgs_combine/make_shapes.cpp`).
Not built by `scram`; not full Combine. **ROOT is not on PATH in a bare session**
(needs `cmsenv`); inputs are hard-coded EOS paths (`/eos/user/a/avijay/...`).

> **The three macros are not connected.** `make_shapes.cpp` &rarr;
> `input_root_file_combine_v4.root`; `make_datacard.cpp` &larr;
> `input_root_file_combine_v3.root`; `make_shapes_final.cpp` is a separate larger
> variant. Running `make_datacard.cpp` as-is reads a **stale `_v3`** file (or
> none). Establish and regenerate the current shape file before touching the
> datacard.

| File | Does | Gap vs AN |
|------|------|-----------|
| `make_shapes.cpp` | fills `HZZ2l2nu_ZZmT` histograms (100 bins, 0&ndash;1000), per background sample, `weight = "puWeight"`, `scale = xsec*lumi/(nMC - 2*nNeg)`; adds `data_obs`, one `signal_alpha` from `signal_m125.root`, and a stitched inclusive `dy_alpha` | `lumi = 59547` (2018, not 2016); shape systematics `_alphaUp/_alphaDown` are **exact clones** of nominal; **single** inclusive histogram — no `ee`/`&mu;&mu;` or `=0j`/`&ge;1j`/`VBF` split; signal is SM m=125, not the 200&ndash;3000 GeV grid; Top / DY taken from **MC**, not the &alpha;- / &gamma;+jets data-driven templates |
| `make_shapes_final.cpp` | larger variant (WJet bins, category struct) | same class of gaps; verify which is current |
| `make_datacard.cpp` | reads `input_root_file_combine_v3.root`, writes `datacard.txt`: `imax 1`, one bin `bin1`, ~24 processes, `shapes * * ... $PROCESS $PROCESS_$SYSTEMATIC`, two nuisances | only `lumi lnN` (value **0.84** — wrong; an `lnN` &lt; 1) and `alpha shape` (currently null, see above); none of the &sect;`systematics.md` nuisances; no category bins; no interference `PhysicsModel` (Eq. 14) |

**To reach the AN model** the scaffold needs: per-category shape inputs; real
&plusmn;1&sigma; systematic templates (`references/systematics.md` &sect;5 names);
data-driven Top/WW (e&mu;) and Z+jets (&gamma;+jets) templates; the high-mass signal
grid with MELA-reweighted width/interference; the Eq. 14 `PhysicsModel`; and the
actual `combine`/`text2workspace.py` limit step (not present in this repo).

---

## 6. Limit procedure (AN &sect;6.2, &sect;9.2)

- **Asymptotic CL_s** 95% upper limits on `sigma(gg->H->ZZ)` and
  `sigma(qq->H->ZZ)` vs `m_H`, per width, expected (&plusmn;1/2&sigma; bands) and
  observed (AN Figs. 49&ndash;51).
- The six categories (`ee`/`&mu;&mu;` &times; `=0j`/`&ge;1j`/`VBF`) are fit jointly;
  per-category limits also shown (Figs. 50&ndash;51).
- The analysis was **blinded** in the note (data replaced by total background in
  Table 23 / `M_T` plots &mdash; "blinded area").
- 2016 result: a narrow scalar excluded over **[300, 3000] GeV** at 95% CL for the
  scenarios tested.

---

## 7. Review checklist

1. Fit variable = `M_T` per Eq. 6; binning/range recorded; `MET >= 125 GeV`
   applied before templating.
2. Templates: signal/WZ/ZZ from MC; Z+jets from &gamma;+jets **data**;
   non-resonant from e&mu; **data** (&alpha;-method) — not MC.
3. Six categories built and combined; ggF and VBF limits kept separate.
4. Signal: high-mass grid, &Gamma; = 5/10/100, MELA reweighting for
   width + Int1/Int2/Int3; Eq. 14 &mu;-parameterization in the datacard.
5. Nuisances from `references/systematics.md` &sect;5 all wired; shape variations are
   real recomputations, not clones.
6. Blinding policy explicit; expected limits with bands; asymptotic CL_s.
7. Repo `higgs_combine/` gaps (&sect;5) tracked as open items.

---

## 8. Evidence summary

| Item | Source | Established? |
|------|--------|--------------|
| Shape fit on `M_T`; `MET >= 125` optimization | AN &sect;6 | yes — method |
| MELA reweighting for width + Int1/2/3 | AN &sect;3 | yes — method |
| Interference &mu;-model Eq. 14 | AN &sect;9.2 | yes — 2016 form; **[Verify]** current |
| &alpha;- / k-method non-resonant prediction | AN &sect;5.1 | yes — method, 2016 &alpha; values |
| CL_s limits per category, ggF/VBF separate, blinded | AN &sect;6.2, &sect;9.2 | yes |
| `higgs_combine/` implementation | repo macros | early scaffold — divergences in &sect;5 |

## Last verified

- AN-2016/325 &sect;3 / &sect;6 / &sect;9 transcription: skill creation.
- `higgs_combine/` cross-check: repo macros read at skill-creation time (ROOT /
  Combine not executed — not available in a bare session).
- Current UL / Run 2+3 2&ell;2&nu; statistical model: **not consulted — [Verify]**.
