#!/bin/bash
rm -rf build
mkdir build
cd build
export CXX=g++
cmake -DgenericIO_DIR=/sep/bob/genericIO/lib -DCMAKE_INSTALL_PREFIX=/sep/bob/Born ..
make
make install
