setup-root

if [ "lhapdf" != *`module list`* ];
then
	module load lhapdf/6.1
fi

export LHAPDF_DATA_PATH=/home/fynu/amertens/scratch/Reconstruction/LHAPDF-6.1.4/
