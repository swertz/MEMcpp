import urllib
import string
import os
import sys
import LaunchOnCondor

FarmDirectory = "condor"
JobName = "cuba"
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

try:
    if sys.argv[1] == "del":
        print "Clean condor working directories."
        os.system("rm " + os.path.join(FarmDirectory, "inputs/*"))
        os.system("rm " + os.path.join(FarmDirectory, "logs/*"))
except IndexError:
    pass

path = "/home/fynu/swertz/tests_MEM/MEMcpp/"

start_evt = 0
end_evt = -1
evt_per_job = 20
max_evt= 10000
i = 0

while start_evt < max_evt:
	end_evt += evt_per_job
	#LaunchOnCondor.SendCluster_Push(["BASH", os.path.join(path, "condor.sh"), "/home/fynu/amertens/Delphes/delphes/Htt_v3.root", "ttbar_DMEM_MTTbar_noTF_Htt_v3_maxL_" + str(i) + ".root", start_evt, end_evt])
	LaunchOnCondor.SendCluster_Push(["BASH", os.path.join(path, "condor.sh"), os.path.join(path, "data/ttbar.root"), "ttbar_DMEM_MTTbar_noTF_try3_maxL_" + str(i) + ".root", start_evt, end_evt])
	#LaunchOnCondor.SendCluster_Push(["BASH", os.path.join(path, "condor.sh"), os.path.join(path, "data/ttbar_weighted_noTF_isr0_pdfMtop_noMCoPerms_sobol_2000evt.root"), "ttbar_noTF_isr0_pdfMtop_noMCoPerms_as013_cuba_sobol_noSmooth_jacthres_raise30000_" + str(i) + ".root", start_evt, end_evt])
	#LaunchOnCondor.SendCluster_Push(["BASH", path + "/condor.sh", path + "/results/ttbar_noTF_isr0_pdfMtop_noMCoPerms_as013_cuba_Mersenne_noSmooth_50000.root", "ttbar_noTF_isr0_pdfMtop_noMCoPerms_as013_cuba_Mersenne_noSmooth_50000_bis_" + str(i) + ".root", start_evt, end_evt])
	start_evt = end_evt+1
	i += 1

LaunchOnCondor.SendCluster_Submit()
