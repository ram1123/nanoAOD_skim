import uproot
import awkward as ak
import numpy as np

def extract_2l2nu_sync_info(input_file, output_file):
    """
    Extracts relevant 2l2nu information from the NanoAOD file, applies selection cuts,
    and writes to a text file for synchronization.
    """
    # Open the file and access the Events tree
    with uproot.open(input_file) as file:
        events = file["Events"]

        # Load relevant branches
        branches = [
            "run",
            "luminosityBlock",
            "event",
            "pTL1", "etaL1", "phiL1",  # Lepton 1
            "pTL2", "etaL2", "phiL2",  # Lepton 2
            "MET_pt", "MET_phi",       # Missing transverse energy
            "massZ1", "pTZ1",          # Z boson properties
            "massZ2", "pTZ2",          # Z boson properties
            "HZZ2l2nu_ZZmT", "HZZ2l2nu_ZZpT",  # ZZ system properties
            "HZZ2l2qNu_nJets",         # Jet multiplicity
            "HZZ2l2nu_minDPhi_METAK4", # Delta Phi(Jet, MET)
            "passZZ2l2nuSelection"   # Preselection flag
            # "btagWeight_DeepCSVB"      # b-tagging weight
        ]

        # Check if branches exist
        available_branches = events.keys()
        missing_branches = [branch for branch in branches if branch not in available_branches]
        if missing_branches:
            print(f"Warning: Missing branches: {missing_branches}")
            branches = [branch for branch in branches if branch in available_branches]

        # Load available branches
        data = events.arrays(branches, library="ak")

    # Apply event selection cuts
    selected_events = (
        (data["passZZ2l2nuSelection"] == 1)  &
        (data["pTL1"] > 25) & (abs(data["etaL1"]) < 2.4) &
        (data["pTL2"] > 25) & (abs(data["etaL2"]) < 2.4) &
        (data["pTZ1"] > 55) & (abs(data["massZ1"] - 91) < 15) &
        (data["MET_pt"] > 100)
        # (data.get("HZZ2l2nu_minDPhi_METAK4", ak.ones_like(data["run"])) > 0.5) &  # Default to no cut if branch missing
        # (data.get("btagWeight_DeepCSVB", ak.zeros_like(data["run"])) < 0.5)       # Default to no veto if branch missing
    )

    # Select data passing cuts
    filtered_data = data[selected_events]

    # Write the filtered data to a sync file
    with open(output_file, "w") as f:
        f.write("# Run:Lumi:Event | L1_pt, eta, phi | L2_pt, eta, phi | MET_pt, phi | Z_mass, Z_pT\n")
        for run, lumi, event, pt1, eta1, phi1, pt2, eta2, phi2, met_pt, met_phi, mass_z1, pt_z1 in zip(
            filtered_data["run"],
            filtered_data["luminosityBlock"],
            filtered_data["event"],
            filtered_data["pTL1"], filtered_data["etaL1"], filtered_data["phiL1"],
            filtered_data["pTL2"], filtered_data["etaL2"], filtered_data["phiL2"],
            filtered_data["MET_pt"], filtered_data["MET_phi"],
            filtered_data["massZ1"], filtered_data["pTZ1"]
        ):
            f.write(
                f"{run}:{lumi}:{event} | "
                f"{pt1:.2f}, {eta1:.2f}, {phi1:.2f} | "
                f"{pt2:.2f}, {eta2:.2f}, {phi2:.2f} | "
                f"{met_pt:.2f}, {met_phi:.2f} | "
                f"{mass_z1:.2f}, {pt_z1:.2f}\n"
            )

    print(f"Filtered synchronization data written to {output_file}")



# Example usage
input_file = "skimmed_nano.root"
output_file = "2l2nu_sync_info.txt"
extract_2l2nu_sync_info(input_file, output_file)
