#module load root/5.34.09-sl6_gcc44 
setup-root
module load lhapdf/6.1
module load gcc/gcc-4.9.1-sl6_amd64

export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/home/fynu/amertens/scratch/Reconstruction/LHAPDF-6.1.4/
export LIBRARY_PATH=$LIBRARY_PATH:/home/fynu/swertz/soft/Cuba-4.2/:/home/fynu/swertz/storage/Delphes/Delphes-3.1.2/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fynu/swertz/soft/Cuba-4.2/:/home/fynu/swertz/storage/Delphes/Delphes-3.1.2/
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/home/fynu/swertz/storage/Delphes/Delphes-3.1.2/:/usr/include/boost141/:/home/fynu/swertz/soft/Cuba-4.2/
