cmake --version || "CMake is not found!. Please install CMake and try again" || exit 3
cmake --version || "nproc is not found!. Please install nproc and try again " || exit 3

OS_NAME="$(uname)"
CPU_NAME="$(uname -m)"

if [[ "${OS_NAME}" == 'Darwin' ]]; then
  echo "MacOS detected"
  which -s brew
  if [[ $? != 0 ]] ; then
    # Install Homebrew
    echo "Homebrew is required on MacOS but could not detected - please install it and try again"
    exit 3
  else
    echo "Homebrew found!"
  fi

  HOMEBREW_PREFIX="$(brew info gcc | grep Cellar | awk '{print $1}')"

  if [[ "${CPU_NAME}" == 'arm64' ]]; then # identical to intel for now, but we might want to add other things here
    echo "Configuring compilers for Apple Silicon"
    GCC_VERSION="$(brew info gcc | grep Cellar | awk -F"/gcc/" '{print $2}' | awk -F"." '{print $1}')"
    GCC=$HOMEBREW_PREFIX/bin/gcc-$GCC_VERSION
    GPLUSPLUS=$HOMEBREW_PREFIX/bin/g++-$GCC_VERSION
    GFORTRAN=$HOMEBREW_PREFIX/bin/gfortran
    Threads="$(nproc --all)"
  else
    echo "Configuring compilers for intel CPU"
    GCC_VERSION="$(brew info gcc | grep Cellar | awk -F"/gcc/" '{print $2}' | awk -F"." '{print $1}')"
    GCC=$HOMEBREW_PREFIX/bin/gcc-$GCC_VERSION
    GPLUSPLUS=$HOMEBREW_PREFIX/bin/g++-$GCC_VERSION
    GFORTRAN=$HOMEBREW_PREFIX/bin/gfortran
    Threads="$(nproc --all)"
  fi
else # assuming GNU/Linux
  echo "Configuring compilers for GNU/Linux"
  GCC="$(which gcc)"
  GPLUSPLUS="$(which g++)"
  GFORTRAN="$(which gfortran)"
  Threads="$(nproc --all)"
fi

export CC=$GCC
export CXX=$GPLUSPLUS
export FC=$GFORTRAN

$CC --version  || echo "gcc not found. Please install gcc then re-run this script " || exit 3
$CXX --version  || echo "g++ not found. Please install g++ then re-run this script " || exit 3
$FC --version  || echo "gfortran not found. Please install gfortran then re-run this script " || exit 3

mainDir=$PWD
dependencyDir=$mainDir/dependencies

export LD_LIBRARY_PATH=$mainDir/dependencies/lib:$LD_LIBRARY_PATH
export LDFLAGS="-L$mainDir/dependencies/lib -L$mainDir/dependencies/lib64"
export CPPFLAGS="-I$mainDir/dependencies/include"


cd $dependencyDir
# Clipper requires fftw2 for now

if [[ ! -f include/fftw.h ]]; then
#  cd $dependencyDir

  if [[ ! -d fftw ]]; then
    mkdir fftw
    cd fftw
    wget ftp://ftp.fftw.org/pub/fftw/fftw-2.1.5.tar.gz
    tar -zxvf fftw-2.1.5.tar.gz
  else
    cd fftw
  fi

  cd fftw-2.1.5
  curl https://git.savannah.gnu.org/cgit/config.git/plain/config.guess --output config.guess
  curl https://git.savannah.gnu.org/cgit/config.git/plain/config.sub --output config.sub
  CC=$GCC CXX=$GPLUSPLUS ./configure CXXFLAGS='-g -O2 -w -std=c++11' CCFLAGS='-g -O2 -w' --prefix=$dependencyDir --enable-single --enable-float --enable-shared F77=gfortran
  make
  make install
fi

cd $dependencyDir
if [[ ! -f include/fftw.h ]]; then
echo "fftw2 installation FAILED. Cannot continue without fftw2."
exit 3
fi

if [[ ! -d include/mmdb2 ]]; then
cd $dependencyDir
if [[  -d mmdb2 ]]; then
rm -rf mmdb2
fi
mkdir mmdb2
cd mmdb2
wget https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/dependencies/mmdb2-2.0.19.tar.gz
tar -zxvf mmdb2-2.0.19.tar.gz
cd mmdb2-2.0.19
CC=$GCC CXX=$GPLUSPLUS ./configure CXXFLAGS='-g -O2 -w -std=c++11' CFLAGS='-g -O2 -w' --prefix=$dependencyDir --disable-Werror --enable-silent-rules --enable-shared --disable-static
make
make install
fi
cd $dependencyDir
if [[ ! -d include/mmdb2 ]]; then
echo "mmdb2 installation ... falied. We can not continue the rest of the installation steps."
exit 3
fi

if [[ ! -d include/ccp4 ]]; then
cd $dependencyDir
if [[  -d libccp4 ]]; then
rm -rf libccp4
fi
mkdir libccp4
cd libccp4
bzr checkout http://fg.oisin.rc-harwell.ac.uk/anonscm/bzr/libccp4/trunk
cd trunk
CC=$GCC CXX=$GPLUSPLUS ./configure CXXFLAGS='-g -O2 -w -std=c++11' CCFLAGS='-g -O2 -w' \--prefix=$dependencyDir --disable-Werror --enable-silent-rules --enable-shared --disable-static
make
make install
fi
cd $dependencyDir
if [[ ! -d include/ccp4 ]]; then
echo "CCP4 (libccp4) installation ... falied. We can not continue the rest of the installation steps."
exit 3
fi

if [[ ! -d share/ccp4srs ]]; then
cd $dependencyDir/share
if [[  -d ccp4srs ]]; then
rm -rf ccp4srs
fi
mkdir ccp4srs
cd ccp4srs
tar -zxvf $dependencyDir/ccp4srs-data-20180406.tar.gz --directory $dependencyDir/share/ccp4srs
fi
cd $dependencyDir

if [[ ! -d include/ccp4srs ]]; then
cd $dependencyDir
if [[  -d ccp4srs ]]; then
rm -rf ccp4srs
fi
mkdir ccp4srs
cd ccp4srs
bzr checkout http://fg.oisin.rc-harwell.ac.uk/anonscm/bzr/ccp4srs/trunk
cd trunk
 CC=$GCC CXX=$GPLUSPLUS ./configure CXXFLAGS='-g -O2 -w -std=c++11' CCFLAGS='-g -O2 -w' --prefix=$dependencyDir --disable-Werror --enable-silent-rules --enable-shared --disable-static
make
make install
fi
cd $dependencyDir
if [[ ! -d include/ccp4srs ]]; then
echo "ccp4srs installation ... falied. We can not continue the rest of the installation steps."
exit 3
fi

if [[ ! -d lib/data ]]; then
cd $dependencyDir/lib
if [[  -d data ]]; then
rm -rf data
fi
mkdir data
cp $dependencyDir/syminfo.lib $dependencyDir/symop.lib $dependencyDir/symop_old.lib data/
cd data
mkdir monomers
bzr checkout http://fg.oisin.rc-harwell.ac.uk/anonscm/bzr/monomers/trunk/
cd trunk
mv * $dependencyDir/lib/data/monomers
fi
cd $dependencyDir/lib/data
rm -rf trunk
cd $dependencyDir

if [[ ! -d include/clipper ]]; then
cd $dependencyDir
if [[  -d clipper ]]; then
rm -rf clipper
fi
mkdir clipper
cd clipper
bzr checkout http://fg.oisin.rc-harwell.ac.uk/anonscm/bzr/clipper/trunk/
cd trunk
patch -p0 < $dependencyDir/localfftw2_clipper.patch
 CC=$GCC CXX=$GPLUSPLUS ./configure CXXFLAGS='-g -O2 -w -std=c++11' \
 CCFLAGS='-g -O2 -w' \
 --enable-mmdb=$dependencyDir \
 --enable-mmdb \
 --enable-minimol \
 --enable-cif \
 --enable-ccp4=$dependencyDir \
 --enable-ccp4 \
 --prefix=$dependencyDir \
 --disable-Werror \
 --enable-shared \
 --enable-cns \
 --with-fftw2-prefix=$dependencyDir
make
make install
fi
cd $dependencyDir
if [[ ! -d include/clipper ]]; then
echo "clipper installation ... falied. We can not continue the rest of the installation steps."
exit 3
fi

if [[ -f include/fftw.h ]]; then
echo "fftw is installed. Remove include/fftw.h and fftw folder if you want to re-install"
fi
if [[  -d include/ccp4srs ]]; then
echo "ccp4srs is installed. Remove include/ccp4srs and ccp4srs folder  if you want to re-install"
fi
if [[  -d include/mmdb2 ]]; then
echo "mmdb2 is installed. Remove include/mmdb2 and mmdb2 folder if you want to re-install"
fi
if [[  -d include/mmdb ]]; then
echo "mmdb is installed. Remove include/mmdb and mmdb folder if you want to re-install"
fi
if [[  -d include/clipper ]]; then
echo "clipper is installed. Remove include/clipper and clipper folder  if you want to re-install"
fi


if [[ "$OSTYPE" == "linux-gnu" ]]; then
echo "GCC used in compiling: (Please note that using different GCC versions might cause segmentation fault error or a compiler error)"
if [[ -f lib/libfftw.so ]]; then
echo "fftw "
readelf -p .comment lib/libfftw.so | grep  "GCC" | sort --unique
fi
if [[ -f lib/libccp4c.so ]]; then
echo "libccp4c "
readelf -p .comment lib/libccp4c.so | grep  "GCC" | sort --unique
fi
if [[ -f lib/libmmdb.so ]]; then
echo "mmdb "
readelf -p .comment lib/libmmdb.so | grep  "GCC" | sort --unique
fi
if [[ -f lib/libmmdb2.so ]]; then
echo "mmdb2 "
readelf -p .comment lib/libmmdb2.so | grep  "GCC" | sort --unique
fi
if [[ -f lib/libclipper-core.so  ]]; then
echo "clipper "
readelf -p .comment lib/libclipper-core.so | grep  "GCC" | sort --unique
fi
fi
