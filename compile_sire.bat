@echo off

set OPTIND=0
set OPTION_h=
set OPTION_help=
set OPTION_install=
set OPTION_G=
set OPTION_clean=

:_loop_args
    set _getopt_arg=%1
    if defined _getopt_arg (
        call :parse_arg %_getopt_arg%
        if not errorlevel 1 (
            shift /1
            set /a OPTIND+=1
            goto :_loop_args
        )
    )
set _getopt_arg=
set _getopt_opt=
set _getopt_value=
set _getopt_bool=
goto :main

:parse_arg
set _getopt_opt=
set _getopt_value=%~1
set _getopt_bool=1
if not "%_getopt_value:~0,1%"=="/" goto :_parse_arg_ret
set _getopt_value=%_getopt_value:~1%
:_loop_arg_chars
    if not defined _getopt_value goto :_parse_arg_ret
    if "%_getopt_value:~0,1%"==":" (
        set _getopt_bool=
        set _getopt_value=%_getopt_value:~1%
        goto :_parse_arg_ret
    ) else (
        set _getopt_opt=%_getopt_opt%%_getopt_value:~0,1%
        set _getopt_value=%_getopt_value:~1%
    )
    goto :_loop_arg_chars
:_parse_arg_ret
if not defined _getopt_opt exit /B 1
if defined _getopt_bool (
    call :set_opt %_getopt_opt% %_getopt_bool%
) else (
    call :set_opt %_getopt_opt% %_getopt_value%
)
exit /B 0

:set_opt
set OPTION_%1=%~2
exit /B 0

:main
if defined OPTION_h set OPTION_help=1
if defined OPTION_help (
    echo compile_sire.bat /h or /help shows help. 1>&2
    echo compile_sire.bat /install:\path\to\sire.app will install Sire in \path\to\sire.app 1>&2
    echo compile_sire.bat /clean completely cleans the build directory. 1>&2
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
