#!/bin/sh
#SBATCH --job-name=sigma_male_grindset          # Job name
#SBATCH --time=48:00:00                         # Time limit hrs:min:sec
#SBATCH --mem=16gb                              # Total memory limit
#SBATCH ----threads-per-core=16                 # nThreads
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hb1115@york.ac.uk           # Where to send mail
#SBATCH --account=chem-structglyco-2019

module purge

module load devel/CMake/3.16.4-GCCcore-9.3.0
module load lang/Python/3.8.2-GCCcore-9.3.0

GCC="$(which gcc)"
GPLUSPLUS="$(which g++)"
echo "GCC: $GCC"
echo "GCC: $GPLUSPLUS"
GFORTRAN="$(which gfortran)"
echo "GFORTRAN: $GFORTRAN"

export CC=$GCC
export CXX=$GPLUSPLUS
export FC=$GFORTRAN

mainDir=$PWD
cd $mainDir
projectSigmaDir=$mainDir/project_sigma

source $mainDir/privateerpython/bin/activate
source $mainDir/ccp4.envsetup-sh

cd $projectSigmaDir
cd pdb_mirror

rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb

cd $projectSigmaDir
python clustering_analysis/prepare_glycans_for_clustering.py
cd $mainDir
ls
git status