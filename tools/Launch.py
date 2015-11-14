import urllib
import string
import os
import sys
from ROOT import TFile
import LaunchOnCondor

FarmDirectory = "condor"
JobName = "cuba"
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

path = "/home/fynu/swertz/tests_MEM/MEMcpp/"
data_path = "/home/fynu/swertz/storage/Delphes/condorDelphes/"

process = "TTbar_qCut50"

dir_path = process + "/condor/output/output_selected.root"

file_path = os.path.join(data_path, dir_path)

file = TFile(file_path)
tree = file.Get("t")

start_evt = 0
end_evt = -1
evt_per_job = 100
max_evt= tree.GetEntries()
file.Close()
i = 0

while start_evt < max_evt:
	end_evt += evt_per_job
	LaunchOnCondor.SendCluster_Push(["BASH", os.path.join(path, "tools/condor.sh"), file_path, "output_weighted_" + process + "_" + str(i) + ".root", os.path.join(path, "../binnedTF/TF_generator/Control_plots_hh_TF.root"), start_evt, end_evt])
	start_evt = end_evt+1
	i += 1

LaunchOnCondor.SendCluster_Submit()
