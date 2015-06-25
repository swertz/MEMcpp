import urllib
import string
import os
import sys
import LaunchOnCondor

FarmDirectory = "condor"
JobName = "cuba_Htt"
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

path = "/home/fynu/swertz/tests_MEM/MEMcpp/"

start_evt = 0
end_evt = -1
evt_per_job = 10
max_evt= 2000
i = 0

while start_evt < max_evt:
	end_evt += evt_per_job
	#LaunchOnCondor.SendCluster_Push(["BASH", os.path.join(path, "condor.sh"), os.path.join(path, "data/ttbar.root"), "ttbar_DMEM_MTTbar_noTF_try0_" + str(i) + ".root", start_evt, end_evt])
	#LaunchOnCondor.SendCluster_Push(["BASH", os.path.join(path, "condor.sh"), os.path.join(path, "data/ttbar_weighted_binnedTF_correctWidths_isr0_pdfMtop_noMCoPerms_sobol_10000evt.root"), "ttbar_DMEM_MTTbar_binnedTF_isr0_pdfMtop_noMCoPerms_as013_cuba_sobol_bugFix_" + str(i) + ".root", os.path.join(path, "../binnedTF/TF_generator/Control_plots_hh_TF.root"), start_evt, end_evt])
	LaunchOnCondor.SendCluster_Push(["BASH", os.path.join(path, "condor.sh"), "/home/fynu/amertens/Delphes/delphes/Htt_v3.root", "ttbar_DMEM_MTTbar_binnedTF_isr0_pdfMtop_noMCoPerms_as013_cuba_sobol_bugFix_Htt_" + str(i) + ".root", os.path.join(path, "../binnedTF/TF_generator/Control_plots_hh_TF.root"), start_evt, end_evt])
	start_evt = end_evt+1
	i += 1

LaunchOnCondor.SendCluster_Submit()
