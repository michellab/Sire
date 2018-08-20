#!/bin/bash

cd cmake-3.6.0
./configure --prefix=$HOME/local
make -j 4 install

