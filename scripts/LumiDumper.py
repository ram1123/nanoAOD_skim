import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from collections import defaultdict
import json

class LumiDumper(Module):
    def __init__(self):
        self.lumi_dict = defaultdict(list)

    def analyze(self, event):
        run = getattr(event, "run")
        lumi = getattr(event, "luminosityBlock")

        if lumi not in self.lumi_dict[run]:
            self.lumi_dict[run].append(lumi)

        return True

    def endJob(self):
        # Convert lumi sections into JSON format
        for run, lumis in self.lumi_dict.items():
            self.lumi_dict[run] = [[min(lumis), max(lumis)]]

        # Dump the dictionary into a JSON file
        with open('lumi.json', 'w') as f:
            json.dump(self.lumi_dict, f)
