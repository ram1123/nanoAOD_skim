# Example Usage:
#
# python Utils/nanoAOD_condor_resubmit.py -d condor_logs/SkimNanoAOD_2022_ZXCR/240312_135155/ -s /eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCR/ -i submit_condor_jobs_lnujj_SkimNanoAOD_2022_ZXCR.txt -n 1
#

import os
from optparse import OptionParser
import ROOT
from ROOT import TFile
import re
import subprocess

DEBUG = False

def files_to_remove(files):
    print('Checking for corrupted files')
    filelist_to_remove = []
    for file in files:
        if DEBUG: print("Checking file: {}".format(file))
        try:
            tfile = TFile.Open(file);
        except:
            pass
        if tfile:
            if (tfile.IsZombie()):
                filelist_to_remove.append(file.replace("_Skim",""))
        else:
            print('File could not be opened, adding it to missing files')
            filelist_to_remove.append(file.replace("_Skim",""))

    print('Number of corrupted files: {}'.format(len(filelist_to_remove)))
    if DEBUG: print("Files to remove: {}".format(filelist_to_remove))
    return filelist_to_remove

def list_files(file_name):
    with open(file_name) as f:
        datafile = f.readlines()
        file_list = []
        for line in datafile:
            if "USER_" in line:
                start_pt = line.find("output")
                end_pt = line.find(".root", start_pt + 1)
                res = line[start_pt:end_pt+5]
                if "/" in res :
                    res = res[res.find("/")+1:]
                file_list.append(res)
    return file_list

def get_files_from_jdl(path_jdl):
    """Open the JDL file and grab the list of root files present in the JDL file

    Args:
        path_jdl (str): JDL file with its path

    Returns:
        list: list having list of all root files present in the JDL file
    """
    with open(path_jdl) as myfile:
        content = myfile.read()

    flist = re.findall("[a-z0-9A-Z-]+.root", content)
    return flist

def list_root(directory):
    """Get the list of root files present in the directory

    Args:
        directory (str): path of directory from where we need to get the list of root files

    Returns:
        list: list of lists having root files and root files with their path
    """
    flist = []
    flistWithPath = []
    for root, directories, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith("_Skim.root"):
                flist.append(filename.replace("_Skim",""))
                flistWithPath.append(os.path.join(root,filename)) # Append file name with path
    return flist,flistWithPath

def submit_missing(InputJdlFile,resubmit=True):
    bashCommand = "condor_submit {}".format(InputJdlFile)
    if resubmit :
        print('Resubmitting now!')
        os.system(bashCommand)
    else :
        print('Ready to resubmit, please set resubmit to True if you are ready : ')
        print(bashCommand)

def get_condor_job_details(job_id):
    """
    Grabs a list of all Condor jobs associated with a specific job ID,
    then extracts the specified arguments from each job's command.

    Parameters:
    - job_id: The ID of the Condor jobs to filter by.

    Returns:
    A list of tuples, each containing the extracted arguments for each job matching the job ID.
    """
    # Execute the condor_q command and capture its output
    command = 'condor_q -nobatch'
    result = subprocess.check_output(command, shell=True)

    # Decode result to convert bytes to str (for Python 3 compatibility)
    result_str = result.decode('utf-8')
    if DEBUG: print("{}\n {}\n{}".format("="*51,result_str,"="*51))

    # Initialize an empty list to hold the extracted details
    job_details = []

    # Use regular expression to parse each line of the condor_q output
    pattern = re.compile(r'({}\.\d+)'.format(job_id) + r'.*?sh\s+.*?\s+(.*?/eos/.*?)\s.*?([a-f0-9\-]+)$')
    if DEBUG: print("{}\n {}\n{}".format("="*51,pattern,"="*51))

    for line in result_str.split('\n'):
        match = pattern.search(line)
        if match:
            # Extract the 3rd last argument and the last argument
            path_arg = match.group(2)
            last_arg = match.group(3)
            job_details.append((path_arg, last_arg))

    if DEBUG: print("{}\n {}\n{}".format("="*51,job_details,"="*51))
    return job_details

def prepare_runJobs_missing(FailedJobRootFile,InputJdlFile,CondorLogDir,EOSDir,Resubmit_no):
    if DEBUG: print("FailedJobRootFile: {}".format(FailedJobRootFile))
    if DEBUG: print("InputJdlFile: {}".format(InputJdlFile))
    if DEBUG: print("CondorLogDir: {}".format(CondorLogDir))
    if DEBUG: print("EOSDir: {}".format(EOSDir))

    bashCommand = "cp {0}  original_{0}".format(InputJdlFile)
    if DEBUG: print("copy command: {}".format(bashCommand))
    os.system(bashCommand)

    outjdl_fileName = InputJdlFile.replace(".txt", "_resubmit_"+str(Resubmit_no)+".txt")
    outjdl_file = open(outjdl_fileName,"w")

    # Update the new txt file name in the JDL file
    bashCommand = "sed -i 's/{0}/{1}/g' {2}".format(InputJdlFile, outjdl_fileName, InputJdlFile.replace(".txt",".jdl"))
    if DEBUG: print("sed command: {}".format(bashCommand))
    os.system(bashCommand)

    for RootFiles in FailedJobRootFile:
        if DEBUG: print("Root file to look for in stdout files: {}".format(RootFiles))
        bashCommand = "grep {} {}/*.stdout".format(RootFiles, CondorLogDir)
        if DEBUG: print("grep command: {}".format(bashCommand))
        grep_stdout_files = os.popen(bashCommand).read()
        if DEBUG: print("{}\n{}\n{}".format("="*51,grep_stdout_files,"="*51))

        # Regular expression to match paths ending with .stdout
        stdout_file_pattern = re.compile(r'\S+\.stdout')

        # Search for the .stdout file path in the output
        match = stdout_file_pattern.search(grep_stdout_files)
        # print("===")
        # print("grep command: {}".format(bashCommand))
        # print("{}".format(grep_stdout_files))
        # print("Match: {}".format(match))

        OldRefFile = ""
        if match:
            stdout_file_path = match.group()
            if DEBUG: print(stdout_file_path.strip())
            OldRefFile = stdout_file_path.strip().split("/")[-1].replace(".stdout","").split("_")[-1]
        else:
            if DEBUG: print("No .stdout file path found in the output.")
            OldRefFile = ""
        if DEBUG: print("OldRefFile: {}".format(OldRefFile))

        grepCommand_GetJdlInfo = 'grep  "{}" {}'.format(RootFiles, InputJdlFile)
        if DEBUG: print(grepCommand_GetJdlInfo)
        grep_condor_jdl_part = os.popen(grepCommand_GetJdlInfo).read()
        if DEBUG: print("=="*51)
        if DEBUG: print(grep_condor_jdl_part)
        updateString = grep_condor_jdl_part.replace('$(Process)',OldRefFile+'_$(Process)'+ '_resubmit_' +Resubmit_no)
        if DEBUG: print("=="*51)
        if DEBUG: print(updateString)
        if DEBUG: print("=="*51)
        outjdl_file.write(updateString)
    outjdl_file.close()
    return outjdl_fileName

def main():
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir",default="",help="log directory path")
    parser.add_option("-s", "--stage-dest", dest="stage_dest",help="directory output files were staged to")
    parser.add_option("-i", "--input", dest="input",default="all_root.jdl",help="input jdl file with all root files present")
    parser.add_option("-r", "--resubmit", action="store_true",  dest="resubmit",default=False,help="resubmit")
    parser.add_option("-n", "--resubmit_no", dest="resubmit_no",default=1,help="resubmit counter")
    parser.add_option("-c", "--condor_job_id", dest="condor_job_id",default="",help="condor job id")

    (options, args) = parser.parse_args()

    if options.stage_dest is not None:
        stageDir = os.path.abspath(options.stage_dest)
    else:
        stageDir = dir

    # Get list of root files from JDL file, this does not have the "_Skim" with the root files
    full_output =  get_files_from_jdl(options.input)
    print("Lenght of root files from JDL file: {}".format(len(full_output)))
    if DEBUG: print("Root files from JDL file: {}".format(full_output))

    # Get list of root files from the directory, this has the "_Skim" with the root files
    # So, the list named "present_output" will have the root files without "_Skim" in the end
    present_output, present_output_WithPath =  list_root(stageDir)
    print("Lenght of root files from EOS directory: {}".format(len(present_output)))
    if DEBUG: print("Root files from EOS directory: {}".format(present_output))
    if DEBUG: print("Root files with path from EOS directory: {}".format(present_output_WithPath))

    # not_finished will contain the list of root files which are not present in the directory
    not_finished = list(set(full_output) - set(present_output))
    print('Number of missing files : {}'.format(len(not_finished)))
    if DEBUG: print('Missing files : {}'.format(not_finished))

    # corrupted_files will contain the list of root files which are present in the directory but are corrupted
    corrupted_files = files_to_remove(present_output_WithPath)
    print('Number of corrupted files : {}'.format(len(corrupted_files)))
    if DEBUG: print('Corrupted files : {}'.format(corrupted_files))

    # Update the list "not_finished" with the list of corrupted files
    not_finished += corrupted_files
    print('Number of missing files (after including the corrupted files) : {}'.format(len(not_finished)))
    if DEBUG: print('Missing files (including the corrupted files) : {}'.format(not_finished))

    # If the condor jobs are still running then remove the files over which condor jobs are running from the list "not_finished"
    if options.condor_job_id:
        details = get_condor_job_details(options.condor_job_id)
        for detail in details:
            # print("{}/{}_Skim.root".format(detail[0], detail[1]))
            not_finished = list(set(not_finished) - set(["{}.root".format(detail[1])]))

        print('Number of missing files (after removing the files over which condor jobs are running) : {}'.format(len(not_finished)))
        if DEBUG: print('Missing files (after removing the files over which condor jobs are running) : {}'.format(not_finished))

    jdlfile = prepare_runJobs_missing(not_finished,options.input,options.dir,stageDir,str(options.resubmit_no))
    print(jdlfile)
    print('Submitting missing jobs : ')
    submit_missing(options.input.replace(".txt",".jdl"),options.resubmit)

if __name__ == "__main__":
    main()
