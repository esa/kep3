#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O mambaforge.sh
export deps_dir=$HOME/local
export PATH="$HOME/mambaforge/bin:$PATH"
bash mambaforge.sh -b -p $HOME/mambaforge
mamba env create -f kep3_devel.yml -q -p $deps_dir
source activate $deps_dir

# First we build and install the kep3 library
mkdir build
cd build
cmake -G "Ninja" ../ -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_BUILD_TYPE=Release -Dkep3_BUILD_TESTS=no -Dkep3_BUILD_BENCHMARKS=no -Dkep3_BUILD_PYTHON_BINDINGS=no -DBoost_NO_BOOST_CMAKE=ON
cmake --build . --target=install --config=Release -- -j 2
# Then we build and install pykep
cd ..
mkdir build_pykep
cd build_pykep
cmake -G "Ninja" ../ -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_BUILD_TYPE=Release -Dkep3_BUILD_TESTS=no -Dkep3_BUILD_BENCHMARKS=no -Dkep3_BUILD_PYTHON_BINDINGS=yes -DBoost_NO_BOOST_CMAKE=ON
cmake --build . --target=install --config=Release -- -j 2

python -c "import pykep.test; pykep.test.run_test_suite()"

set +e
set +x
