:: Sire build script for Windows.

:: Build and install Sire.
python setup.py install --skip-deps

:: Exit immediately on error.
if errorlevel 1 exit 1

:: Remove the build files to free up space.
rmdir /s/q build
