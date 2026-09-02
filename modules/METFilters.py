# MET noise / event filters.
#
# Source (authoritative): CMS JME
#   https://cms-jme-jmar.docs.cern.ch/recommendations/met/noise_filters/
#
# ---------------------------------------------------------------------------
# Run 2 UltraLegacy (2016 / 2017 / 2018)                                #run-2
#   Flag_goodVertices
#   Flag_globalSuperTightHalo2016Filter
#   Flag_HBHENoiseFilter
#   Flag_HBHENoiseIsoFilter
#   Flag_EcalDeadCellTriggerPrimitiveFilter
#   Flag_BadPFMuonFilter
#   Flag_BadPFMuonDzFilter
#   Flag_hfNoisyHitsFilter           -- "Optional"; applied here
#   Flag_eeBadScFilter
#   Flag_ecalBadCalibFilter          -- 2017-2018 ONLY (not in the 2016 UL recommendation)
#   Flag_BadChargedCandidateFilter   -- "Not recommended": NOT applied
#
# ---------------------------------------------------------------------------
# Run 3 (2022-2026)                                                     #run-3
#   Flag_goodVertices
#   Flag_globalSuperTightHalo2016Filter
#   Flag_EcalDeadCellTriggerPrimitiveFilter
#   Flag_BadPFMuonFilter
#   Flag_BadPFMuonDzFilter
#   Flag_hfNoisyHitsFilter
#   Flag_eeBadScFilter
#   Flag_ecalBadCalibFilter          -- 2024: use the stored flag.
#                                       2025-2026: recommendation under review.
#                                       2022-2023 PROMPT-reco DATA: do NOT use the stored
#                                       flag (NanoAOD lacks the ECAL rec hits) - apply the
#                                       special event-level selection below instead.
#   HBHE (Flag_HBHENoiseFilter / Flag_HBHENoiseIsoFilter) are NOT required for Run 3.
#
# 2022-2023 prompt-reco ECAL bad crystal (detector ID 838871812), DATA ONLY,
# runs 362433-367144: reject the event if
#     PuppiMET_pt > 100 GeV  AND  >= 1 AK4 jet with
#         Jet_pt  > 50
#         -0.5 < Jet_eta < -0.1
#         -2.1 < Jet_phi < -1.8
#         (Jet_neEmEF > 0.9 or Jet_chEmEF > 0.9)
#     (no Jet_jetId requirement). Not for simulation; not needed for re-reco datasets.
#     Estimated loss of good data events < 0.2%.
# ---------------------------------------------------------------------------

_MET_FILTERS_RUN2 = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_eeBadScFilter",
]

_MET_FILTERS_RUN3 = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_eeBadScFilter",
]

# Set True ONLY when processing 2022/2023 *prompt-reconstruction* DATA (not re-reco,
# not MC). When True the stored Flag_ecalBadCalibFilter is replaced, for runs
# 362433-367144, by the special event-level selection described above.
IS_2022_2023_PROMPT_RECO_DATA = False


def _fail_2022_2023_promptreco_ecal(event):
    """2022-2023 prompt-reco ECAL bad-crystal special handling (data only).
    Returns True if the event must be REJECTED."""
    if not (362433 <= event.run <= 367144):
        return False
    if event.PuppiMET_pt <= 100:
        return False
    for i in range(event.nJet):
        if (event.Jet_pt[i] > 50
                and -0.5 < event.Jet_eta[i] < -0.1
                and -2.1 < event.Jet_phi[i] < -1.8
                and (event.Jet_neEmEF[i] > 0.9 or event.Jet_chEmEF[i] > 0.9)):
            return True
    return False


def passFilters(event, year, isMC=None, debug=False):
    if year in (2016, 2017, 2018):
        flags = list(_MET_FILTERS_RUN2)
        if year in (2017, 2018):
            flags.append("Flag_ecalBadCalibFilter")
        for flag in flags:
            if getattr(event, flag) == 0:
                if debug:
                    print("DEBUG: MET filter FAILED: {}".format(flag))
                return False
            if debug:
                print("DEBUG: {} passed".format(flag))
        return True

    elif year in (2022, 2023, 2024, 2025, 2026):
        for flag in _MET_FILTERS_RUN3:
            if getattr(event, flag) == 0:
                if debug:
                    print("DEBUG: MET filter FAILED: {}".format(flag))
                return False
            if debug:
                print("DEBUG: {} passed".format(flag))

        # ECAL bad calibration.
        if year in (2022, 2023) and IS_2022_2023_PROMPT_RECO_DATA and isMC is False:
            if _fail_2022_2023_promptreco_ecal(event):
                if debug:
                    print("DEBUG: 2022/2023 prompt-reco ECAL special selection FAILED")
                return False
            if debug:
                print("DEBUG: 2022/2023 prompt-reco ECAL special selection passed")
        else:
            # 2024 and (2022/2023 re-reco or MC): use the stored flag.
            # 2025-2026: JME recommendation under review -- [Verify].
            if getattr(event, "Flag_ecalBadCalibFilter") == 0:
                if debug:
                    print("DEBUG: MET filter FAILED: Flag_ecalBadCalibFilter")
                return False
            if debug:
                print("DEBUG: Flag_ecalBadCalibFilter passed")
        return True

    else:
        print("ERROR: Invalid year: {}".format(year))
        exit(1)
