#!/bin/bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release "$@" ..
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)