#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge3.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniforge3/bin:$PATH"
bash miniforge3.sh -b -p $HOME/miniforge3
conda env create --file=kep3_devel.yml -q -p $deps_dir
source activate $deps_dir
conda install lcov -y
conda list

export CXXFLAGS="$CXXFLAGS --coverage"

mkdir build
cd build

cmake -G "Ninja" ../ -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_BUILD_TYPE=Debug -Dkep3_BUILD_TESTS=yes -DBoost_NO_BOOST_CMAKE=ON

cmake --build . -- -v

ctest -j4 -VV

# Create lcov report
lcov --capture --directory . --output-file coverage.info

set +e
set +x