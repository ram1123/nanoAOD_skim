import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
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
            lumis = sorted(lumis)
            ranges = []
            start = lumis[0]
            prev = lumis[0]

            for l in lumis[1:]:
                if l == prev + 1:
                    prev = l
                else:
                    ranges.append([start, prev])
                    start = l
                    prev = l
            ranges.append([start, prev])

            self.lumi_dict[run] = ranges

        with open('lumi.json', 'w') as f:
            json.dump(self.lumi_dict, f, indent=2)

# Input ROOT file
input_file = "/eos/user/a/avijay/Mergedrootfiles_hzz_new/file_26march/data/data_2018_noDuplicates.root"

# Run PostProcessor
p = PostProcessor(
    ".",
    [input_file],
    modules=[LumiDumper()],
    noOut=True
)

p.run()