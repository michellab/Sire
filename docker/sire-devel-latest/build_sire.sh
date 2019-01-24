#!/bin/bash
unset PYTHONPATH
cd Sire
./compile_sire.sh | tee ../compile_sire.log
