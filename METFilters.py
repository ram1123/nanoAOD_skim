def passFilters(event, year, debug=False):
    # Referece: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2018_2017_data_and_MC_UL
    if year == 2017 or year == 2018:
        # For 2017 and 2018
        if event.Flag_goodVertices == 0:
            return False
        if debug: print("DEBUG: Flag_goodVertices passed")
        if event.Flag_globalSuperTightHalo2016Filter == 0:
            return False
        if debug: print("DEBUG: Flag_globalSuperTightHalo2016Filter passed")
        if event.Flag_HBHENoiseFilter == 0:
            return False
        if debug: print("DEBUG: Flag_HBHENoiseFilter passed")
        if event.Flag_HBHENoiseIsoFilter == 0:
            return False
        if debug: print("DEBUG: Flag_HBHENoiseIsoFilter passed")
        if event.Flag_EcalDeadCellTriggerPrimitiveFilter == 0:
            return False
        if debug: print("DEBUG: Flag_EcalDeadCellTriggerPrimitiveFilter passed")
        if event.Flag_BadPFMuonFilter == 0:
            return False
        if debug: print("DEBUG: Flag_BadPFMuonFilter passed")
        if event.Flag_BadPFMuonDzFilter == 0:
            return False
        if debug: print("DEBUG: Flag_BadPFMuonDzFilter passed")
        if event.Flag_hfNoisyHitsFilter == 0:
            return False
        if debug: print("DEBUG: Flag_hfNoisyHitsFilter passed")
        if event.Flag_BadChargedCandidateFilter == 0:
            return False
        if debug: print("DEBUG: Flag_BadChargedCandidateFilter passed")
        if event.Flag_eeBadScFilter == 0:
            return False
        if debug: print("DEBUG: Flag_eeBadScFilter passed")
        if event.Flag_ecalBadCalibFilter == 0:
            return False
        if debug: print("DEBUG: Flag_ecalBadCalibFilter passed")
        return True
    elif year == 2016:
        # For 2016
        if event.Flag_goodVertices == 0:
            return False
        if debug: print("DEBUG: Flag_goodVertices passed")
        if event.Flag_globalSuperTightHalo2016Filter == 0:
            return False
        if debug: print("DEBUG: Flag_globalSuperTightHalo2016Filter passed")
        if event.Flag_HBHENoiseFilter == 0:
            return False
        if debug: print("DEBUG: Flag_HBHENoiseFilter passed")
        if event.Flag_HBHENoiseIsoFilter == 0:
            return False
        if debug: print("DEBUG: Flag_HBHENoiseIsoFilter passed")
        if event.Flag_EcalDeadCellTriggerPrimitiveFilter == 0:
            return False
        if debug: print("DEBUG: Flag_EcalDeadCellTriggerPrimitiveFilter passed")
        if event.Flag_BadPFMuonFilter == 0:
            return False
        if debug: print("DEBUG: Flag_BadPFMuonFilter passed")
        if event.Flag_BadPFMuonDzFilter == 0:
            return False
        if debug: print("DEBUG: Flag_BadPFMuonDzFilter passed")
        if event.Flag_BadChargedCandidateFilter == 0:
            return False
        if debug: print("DEBUG: Flag_BadChargedCandidateFilter passed")
        if event.Flag_eeBadScFilter == 0:
            return False
        if debug: print("DEBUG: Flag_eeBadScFilter passed")
        if event.Flag_hfNoisyHitsFilter == 0:
            return False
        if debug: print("DEBUG: Flag_hfNoisyHitsFilter passed")
        return True
    else:
        print("ERROR: Invalid year: {}".format(year))
        exit(1)
