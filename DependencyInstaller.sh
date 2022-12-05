cmake --version || "CMake is not found!. Install CMake then re-run this script " || exit 3
gcc --version  || "GCC is not found!. Install GCC then re-run this script " || exit 3

if [[ $(uname -m) == 'arm64' ]]; then
  echo "Configuring compilers for Apple Silicon"
  GCC="/opt/homebrew/bin/gcc-12"
  GPLUSPLUS="/opt/homebrew/bin/g++-12"
  echo "GCC: $GCC"
  echo "GCC: $GPLUSPLUS"
  GFORTRAN="/opt/homebrew/bin/gfortran"
  echo "GFORTRAN: $GFORTRAN"
  Threads="$(nproc --all)"

  export CC=$GCC
  export CXX=$GPLUSPLUS
  export FC=$GFORTRAN
else 
  GCC="$(which gcc)"
  GPLUSPLUS="$(which g++)"
  echo "GCC: $GCC"
  echo "GCC: $GPLUSPLUS"
  GFORTRAN="$(which gfortran)"
  echo "GFORTRAN: $GFORTRAN"
  Threads="$(nproc --all)"
fi

export CC=$GCC
export CXX=$GPLUSPLUS
export FC=$GFORTRAN

mainDir=$PWD
dependencyDir=$mainDir/dependencies

export LD_LIBRARY_PATH=$mainDir/dependencies/lib:$LD_LIBRARY_PATH
export LDFLAGS="-L$mainDir/dependencies/lib -L$mainDir/dependencies/lib64"
export CPPFLAGS="-I$mainDir/dependencies/include"


cd $dependencyDir
# Clipper requires fftw2 for now

if [[ ! -f include/fftw.h ]]; then
  cd $dependencyDir

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
