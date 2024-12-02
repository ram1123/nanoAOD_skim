import os
import sys
from datetime import datetime

CURRENT_DATETIME = datetime.now()

class FileHelper():
  """docstring for FileHelper"""
  def __init__(self, InputDirectoryPath, OutputDirectoryPath, ifeos=False):
    self.ifeos = ifeos
    self.InputDirectoryPath = InputDirectoryPath
    self.OutputDirectoryPath = OutputDirectoryPath
    self.AllSubDirectoryLists = []
    self.ListOfAllInputRootFiles = []
    self.ListOfAllOutputRootFiles = []

  def list_files(self):
    flist = []
    flistWithPath = []
    for root, directories, filenames in os.walk(self.InputDirectoryPath):
      self.AllSubDirectoryLists.append(root)
      # print root, directories, filenames
      for filename in filenames:
        if filename.endswith(".root"):
          fileWithPath = os.path.join(root,filename)  # Get file name with path
          self.ListOfAllInputRootFiles.append(fileWithPath)

  def CreateSameHirarcyDirsSameAsListFiles(self):
    if self.InputDirectoryPath[-1] == "/":
      ReplaceString = self.InputDirectoryPath.split("/")[-2]
    else:
      ReplaceString = self.InputDirectoryPath.split("/")[-1]
    DirectoryPathToReplace = self.InputDirectoryPath.replace(ReplaceString,"")
    for directories in self.AllSubDirectoryLists:
      NewDirecotyName = directories.replace(DirectoryPathToReplace,(self.OutputDirectoryPath+"/").replace("//","/"))
      os.system("eos root://cmseos.fnal.gov mkdir "+NewDirecotyName)

  def GetOutPutFileList(self):
    if self.InputDirectoryPath[-1] == "/":
      ReplaceString = self.InputDirectoryPath.split("/")[-2]
    else:
      ReplaceString = self.InputDirectoryPath.split("/")[-1]
    DirectoryPathToReplace = (self.InputDirectoryPath.replace(ReplaceString,"")).replace("//","/")
    self.ListOfAllOutputRootFiles = [files.replace(DirectoryPathToReplace,(self.OutputDirectoryPath+"/").replace("//","/")) for files in self.ListOfAllInputRootFiles]
