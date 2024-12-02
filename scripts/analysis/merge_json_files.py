# get all files from directory /eos/user/r/rasharma/nanoAOD_ntuplesRun2018_v9/SingleMuon/Run2018B/231105_170505/*.json

import os
import sys
import glob
import json
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--inputdir', action='store', dest='inputdir', default='/eos/user/r/rasharma/nanoAOD_ntuplesRun2018_v9/SingleMuon/Run2018B/231105_170505/', help='input directory')
parser.add_argument('--outputdir', action='store', dest='outputdir', default='/eos/user/r/rasharma/nanoAOD_ntuplesRun2018_v9/SingleMuon/Run2018B/231105_170505/', help='output directory')
parser.add_argument('--outputfile', action='store', dest='outputfile', default='output.root', help='output file name')
parser.add_argument('--inputfile', action='store', dest='inputfile', default='input.json', help='input file name')
parser.add_argument('--outputlog', action='store', dest='outputlog', default='output.log', help='output log file name')

args = parser.parse_args()

inputdir = args.inputdir
outputdir = args.outputdir
outputfile = args.outputfile
inputfile = args.inputfile
outputlog = args.outputlog

# print("inputdir = ", inputdir)
# print("outputdir = ", outputdir)
# print("outputfile = ", outputfile)
# print("inputfile = ", inputfile)
# print("outputlog = ", outputlog)

# get all files from directory
#files = glob.glob(inputdir + "*.json")
files = glob.glob(inputdir + "*.json")

# print("files = ", files)

command = "mergeJSON.py "
for file in files:
    # print("file = ", file)
    command += file + " "
    # print("command = ", command)

print("command = {}".format(command))
#os.system(command + " --output=total.json")
os.system(command + " --output=total.json")
