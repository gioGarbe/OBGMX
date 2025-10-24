#!/bin/bash

# directory where the patched OpenBabel will be installed
# change as needed
# you will need sudo access to this directory
DEFAULT_OB_HOME=/usr/local/ob-gg-2.3.2

# number of threads to use when building the library
# (set to the number of CPU cores, I assume 8 is a good guess in 2025)
BUILD_THREADS=8

#
# NO NEED TO CHANGE ANYTHING BELOW HERE
#

echo ""
echo "INSTALLING OBGMX"
echo ""
read -e -p "Installation directory [enter for ${DEFAULT_OB_HOME}]: " OB_HOME
OB_HOME=${OB_HOME:-$DEFAULT_OB_HOME}     # use the default if the line is empty

if [ -e ${OB_HOME} ]; then
    echo "${OB_HOME} already exists. Aborting installation."
    exit 1
fi

SUDO=""    
if ! mkdir -p ${OB_HOME}; then
    echo "Cannot create ${OB_HOME}, need root privileges. Please insert password."
    SUDO='sudo'
    if ! $SUDO mkdir -p ${OB_HOME}; then
        echo "Wrong password. Aborting installation."
        exit 2
    fi
fi

OB_VERSION="openbabel-2.3.2"
OB_SRC_NAME=${OB_VERSION}".tar.gz"
OB_SRC_URL="https://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz/download"

if [ ! -e ${OB_SRC_NAME} ]; then
    echo "Downloading "${OB_SRC_NAME}
    curl -L -o ${OB_SRC_NAME} ${OB_SRC_URL}
fi

# uncompress openbabel 2.3.2 and cd into the source directory
tar zxvf ${OB_SRC_NAME}
cd ${OB_VERSION}

# apply the patch
cat ../ob-2.3.2-gg-2025.patch | patch -p1

# build and install the library
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=${OB_HOME} \
      -DCMAKE_CXX_STANDARD=11  \
      -DCMAKE_CXX_STANDARD_REQUIRED=ON \
      -DBUILD_SHARED=OFF ../

make -j ${BUILD_THREADS}
echo ""

echo -n "INSTALLING PATCHED OPENBABEL TO "${OB_HOME}
if [ ! -z $SUDO ]; then    
    echo " | root PASSWORD WILL BE REQUIRED."
else
    echo ""
fi
$SUDO make install

# get the CXX compiler that has been used
CXX=`grep CMAKE_CXX_COMPILER CMakeCache.txt|head -n1 | cut -d \= -f 2`
if [[ "$CXX" == *"Xcode"* ]]; then
    # macOS (Tahoe) fix
    CXX=/usr/bin/g++
fi

echo ""
echo "COMPILING OBGMX USING "${CXX}
echo ""
# return to the original directory
cd ../../
# make obgmx
${CXX} -O2 -o obgmx obgmx.cpp -I${OB_HOME}/include/openbabel-2.0 -L${OB_HOME}/lib -lopenbabel -lz

# clean the build directory
echo "Cleaning the build directory"
rm -rf ${OB_VERSION}
