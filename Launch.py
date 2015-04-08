import urllib
import string
import os
import sys
import LaunchOnCondor

FarmDirectory = "condor"
JobName = "test"
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

path = "/home/fynu/swertz/tests_MEM/MEMcpp/"

start_evt = 0
end_evt = -1
evt_per_job = 50
max_evt= 2000
i = 0

while start_evt < max_evt:
	end_evt += evt_per_job
	LaunchOnCondor.SendCluster_Push(["BASH", path + "/condor.sh", path + "/data/ttbar_weighted_noTF_isr0_pdfMtop_noMCoPerms_2000evt.root", "ttbar_noTF_isr0_pdfMtop_noMCoPerms_as013_50000p_750c_50s_" + str(i) + ".root", start_evt, end_evt])
	start_evt = end_evt+1
	i += 1

LaunchOnCondor.SendCluster_Submit()
