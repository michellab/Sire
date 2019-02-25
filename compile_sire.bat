@echo off

set getopt_prm_count=1

set argc=0
for %%x in (%*) do set /A argc+=1
echo %*&echo.

set _myvar=%*

rem Loop through all command line arguments one at a time
:varloop
set isparam=1
for /f "tokens=1*" %%a in ('echo %_myvar%') DO (
    set getopt_prm=%%a
    set _myvar=%%b
    call :paramtype

    rem shift along arguments and rerun loop
    if not "%%b"=="" goto varloop
)
goto :main

:paramtype
rem If first character starts with a - or / it must be an option
if /i "%getopt_prm:~0,1%"=="-" call :option
if /i "%getopt_prm:~0,1%"=="/" call :option 
if /i "%isparam%"=="1" call :param
exit /B 0

:option
    set isparam=0
    rem Set the Equal Index to the position of the colon.  0 means none was found
    for /f %%j in ('findstring %getopt_prm% :') do set getopt_eq_idx=%%j

    rem If the index is GE 0 then we must have a colon in the option.
    if /i "%getopt_eq_idx%"=="0" (call :nocolon) else (call :colon)
    exit /B 0

        :colon
            rem set the OPTION value to the stuff to the right of the colon
            set /a getopt_prm_name_end=%getopt_eq_idx%-2
            call set getopt_prm_name=%%getopt_prm:~1,%getopt_prm_name_end%%%
            call set getopt_prm_value=%%getopt_prm:~%getopt_eq_idx%%%
            set OPTION_%getopt_prm_name%=%getopt_prm_value%
            exit /B 0

        :nocolon
            rem This is a flag, so simply set the value to 1
            set getopt_prm_name=%getopt_prm:~1%
            set getopt_prm_value=1
            set OPTION_%getopt_prm_name%=%getopt_prm_value%
            exit /B 0

:param
    rem There was no / or - found, therefore this must be a paramater, not an option
    set PARAM_%getopt_prm_count%=%getopt_prm%
    set PARAM_0=%getopt_prm_count%
    set /a getopt_prm_count+=1
    exit /B 0

:main
if defined OPTION_h set OPTION_help=1
if defined OPTION_help (
    echo compile_sire.bat /h or /help shows help.
    echo compile_sire.bat /install:\path\to\sire.app will install Sire in \path\to\sire.app
    echo compile_sire.bat /clean completely cleans the build directory.
    exit /B 0
)
if not defined OPTION_install set OPTION_install=C:\sire.app
if not defined OPTION_G set "OPTION_G=Visual Studio 15 2017 Win64"
set "GENERATOR=%OPTION_G%"
set "INSTALL_SIRE_DIR=%OPTION_install%"
if defined OPTION_clean (
    echo rmdir /S /Q build/miniconda.sh build\corelib build\wrapper
    echo del Miniconda3-*-Windows-*.exe
    rmdir 2> NUL /S /Q build/miniconda.sh build\corelib build\wrapper
    del 2> NUL Miniconda3-*-Windows-*.exe
    echo ...all clean
    exit /B 0
)

set "MINICONDA_VERSION=4.5.12"
if [%processor_architecture%]==[AMD64] set MINICONDA_ARCH=x86_64
if [%processor_architecture%]==[x86] set MINICONDA_ARCH=x86
if not defined MINICONDA_ARCH (
    echo Unsupported architecture.
    exit /B 1
)
echo Running an install under Windows
set MINICONDA_INSTALLER=Miniconda3-%MINICONDA_VERSION%-Windows-%MINICONDA_ARCH%.exe
set MINICONDA=https://repo.continuum.io/miniconda/%MINICONDA_INSTALLER%
if exist "%INSTALL_SIRE_DIR%\" (
    echo Install directory "%INSTALL_SIRE_DIR%" already exists. Assuming that miniconda is already installed here.
) else (
    if not exist build\%MINICONDA_INSTALLER% (
        echo ** Downloading miniconda from %MINICONDA%... **
        powershell -Command "Invoke-WebRequest '%MINICONDA%' -OutFile 'build\%MINICONDA_INSTALLER%'"
    )
    echo ** Installing miniconda **
    build\%MINICONDA_INSTALLER% /InstallationType=JustMe /AddToPath=0 /RegisterPython=0 /NoRegistry=1 /S /D=%INSTALL_SIRE_DIR%
)
if exist "%INSTALL_SIRE_DIR%\python.exe" (
    echo ** Running the Python install script... **
    echo ** "%INSTALL_SIRE_DIR%\python.exe" build\build_sire.py -G "%GENERATOR%" **
    call "%INSTALL_SIRE_DIR%\Scripts\activate.bat
    "%INSTALL_SIRE_DIR%\python.exe" build\build_sire.py -G "%GENERATOR%"
    call conda.bat deactivate
    set PROMPT=$P$G
    exit /B 0
) else (
    echo ** FATAL **
    echo ** Cannot find "%INSTALL_SIRE_DIR%\python.exe" **
    echo ** Something went wrong with the miniconda install! **
    echo ** Remove "%INSTALL_SIRE_DIR%", then run compile_sire.bat --clean, then try again **
    exit /B 1
)
