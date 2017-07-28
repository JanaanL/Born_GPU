#!/bin/bash
rm -rf build
mkdir build
cd build
export CXX=g++
cmake -DgenericIO_DIR=/data/sep/bob/test/genericIO/lib -DCMAKE_INSTALL_PREFIX=/data/sep/bob/test/Born ..
make
make install
