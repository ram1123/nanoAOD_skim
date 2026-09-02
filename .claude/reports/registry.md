# CMS Analysis Work Registry

Concise index of durable project work for the **H→ZZ→4l / 2l2q / 2l2nu NanoAOD
skim** (`post_proc.py` / `modules/H4LCppModule.py` / `src/H4LTools.cc`). Not a
conversation log. One row per report; newest first.

| Date | Type | Topic | Outcome | Report |
|---|---|---|---|---|
| 2026-09-02 | Implementation | Systematic uncertainties review + weight branches | Partial — `overallEventWeight` (+PU/prefire up-down) implemented & smoke-tested; PU producer fixed to `puWeight_UL20XX`; lepton/b-tag SF + JES/JER/MET-shape deferred | implementations/2026-09-02_systematics.md |

<!--
Add a row when a report is written under reports/{decisions,investigations,implementations}/.
Keep "Outcome" to one line (status + headline result). Example:
| 2026-09-01 | Investigation | 2l2nu MET-phi correction wiring | Open — corrector built per-event in H4LCppModule; UL16/17/18 only | investigations/2026-09-01_metphi.md |
-->
