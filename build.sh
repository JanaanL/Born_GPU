#!/bin/bash
rm -rf build
mkdir build
cd build
cmake -DgenericIO_DIR=/opt/genericIO/lib -DCMAKE_INSTALL_PREFIX=/opt/Born ..
make
make install
