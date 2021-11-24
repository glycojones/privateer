#!/bin/sh
#SBATCH --time=00:30:00                # Time limit hrs:min:sec
#SBATCH --mem=1000                     # Total memory limit
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hb1115@york.ac.uk   # Where to send mail
#SBATCH --account=chem-structglyco-2019

module purge

module load devel/CMake/3.16.4-GCCcore-9.3.0
module load lang/Python/3.8.2-GCCcore-9.3.0

cmake --version || "CMake is not found!. Install CMake then re-run this script " || exit 3
gcc --version  || "GCC is not found!. Install GCC then re-run this script " || exit 3

GCC="$(which gcc)"
GPLUSPLUS="$(which g++)"
echo "GCC: $GCC"
echo "GCC: $GPLUSPLUS"
GFORTRAN="$(which gfortran)"
echo "GFORTRAN: $GFORTRAN"
Threads="$(nproc --all)"

export CC=$GCC
export CXX=$GPLUSPLUS
export FC=$GFORTRAN

mainDir=$PWD
dependencyDir=$mainDir/dependencies

export LDFLAGS="-L$mainDir/dependencies/lib -L$mainDir/dependencies/lib64"
export CPPFLAGS="-I$mainDir/dependencies/include"

source $mainDir/privateerpython/bin/activate

cd $dependencyDir
if [[ ! -d $dependencyDir/bzr ]]; then
cd $dependencyDir
mkdir bzr
cd bzr
wget https://launchpad.net/bzr/2.7/2.7.0/+download/bzr-2.7.0.tar.gz
tar -zxvf bzr-2.7.0.tar.gz
cd bzr-2.7.0
python2 setup.py install --home ~
fi

if [[ ! -f bzr ]]; then
echo "Bazaar installation ... falied. We can not continue the rest of the installation steps."
exit 3
fi

cd $mainDir

python $mainDir/setup.py install