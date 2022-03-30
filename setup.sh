#!/bin/bash

mkdir deps
cd deps
git submodule add -b master https://github.com/nmwsharp/polyscope.git
git submodule add -b master https://github.com/nmwsharp/geometry-central.git
git submodule update --init --recursive
cd ..
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
