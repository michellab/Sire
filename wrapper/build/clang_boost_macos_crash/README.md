Minimal example that shows a crash with clang++ on MacOS X86

Compile using 

# ./compile.sh

Run the test using

# python test.py

Clean up using

# ./clean.sh

Assumes you have a conda environment activated with clang, clangxx and boost_python installed

Remove the "-Os" in compile.sh and the test will run without error.

