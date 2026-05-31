@echo off
setlocal enabledelayedexpansion

echo ============================================
echo   Install & Check Project Dependencies
echo   MATLAB, WSL (Ubuntu) - DACE, nlhomann-Json
echo ============================================
echo.

:MATLAB
echo Searching for MATLAB installation...

set FOUND_MATLAB=
set MATLAB_SEARCH_DIRS="C:\Program Files\MATLAB" "C:\Program Files (x86)\MATLAB"

for %%D in (%MATLAB_SEARCH_DIRS%) do (
    if exist %%D (
        for /d %%V in (%%D\R*) do (
            if exist "%%V\bin\matlab.exe" (
                set FOUND_MATLAB=%%V\bin\matlab.exe
            )
        )
    )
)

if "%FOUND_MATLAB%"=="" (
    echo [WARNING] MATLAB not found
    echo Please install MATLAB manually from https://www.mathworks.com/
    echo.
    pause
    goto :WSL
) 
    echo MATLAB found at: %FOUND_MATLAB%
    echo Please, ensure you have a valid MATLAB license.
    echo.
    goto :WSL


REM ==========================================================
REM  Check WSL Ubuntu installation and CMAKE installation
REM ==========================================================
:WSL
set "ROOT=%~dp0"
set "LOG=%ROOT%install_log.txt"
del "%LOG%" >nul

REM Detect WSL user/home
for /f "delims=" %%u in ('wsl -d Ubuntu echo $USER') do set "WSLUSER=%%u"
if "%WSLUSER%"=="" (
  echo [ERROR] Could not detect WSL user. >> "%LOG%"
  echo [ERROR] Could not detect WSL user. Ensure Ubuntu WSL is installed -- Download it from the Windows store, then re-run this batch script
  pause
  exit /b 1
)
set "WSLHOME=\\wsl$\Ubuntu\home\%WSLUSER%"
echo Detected WSL user: %WSLUSER%
echo WSL home: %WSLHOME%
echo Log file: %LOG%
echo.

REM Helper macro for WSL execution
set "WSL=cmd /c wsl -d Ubuntu bash -lc"
choice /M "Install ubuntu updates now?"
    if errorlevel 2 (
        echo Skipped Ubuntu updates by user.
        echo.
        goto :CMAKE
    )
    echo Installing Ubuntu updates...
    %WSL% "sudo apt-get update" >> "%LOG%"
    goto :CMAKE

:CMAKE
echo ==== Checking CMake installation ====

REM Check if cmake is already installed
%WSL% "cmake --version" >> "%LOG%"
if errorlevel 1 (
    echo   CMake not found 
    choice /M "Install CMake now?"
    if errorlevel 2 (
        echo Skipped CMake installation by user.
        echo.
        goto :DACE
    )
    
    echo   Installing CMake...
    %WSL% "sudo apt-get install -y cmake" >> "%LOG%"
    if errorlevel 1 (
        echo [ERROR] Failed to install CMake
        echo.
        pause
        exit /b 2
    )
    echo   CMake installed
    echo.
) else (
    echo   CMake already installed
    echo.
)

REM Install build essentials (gcc, g++, make)
echo ==== Installing build essentials (GCC, G++, Make, BLAS, unzip) ====
%WSL% "sudo apt-get install -y build-essential unzip" >> "%LOG%" 
if errorlevel 1 (
    echo [ERROR] Failed to install build essentials
    echo.
    pause
    exit /b 2
) 
    echo Build essentials installed
    echo.
    goto :DACE


REM ==========================================================
REM ==    DACE detection & installation (Linux procedure)
REM ==    Requires -DWITH_ALGEBRAICMATRIX=ON for Astrotools
REM ==========================================================

:DACE
echo ==== Checking DACE installation (/usr/local) ====

REM Header
%WSL% "sh -c 'test -f /usr/local/include/dace/dace.h || test -f /usr/include/dace/dace.h'" >> "%LOG%"
set DACE_HDR=%ERRORLEVEL%

REM Library
%WSL% "sh -c '[ -f /usr/local/lib/libdace_s.a ] || [ -f /usr/local/lib/libdace.so ] || [ -f /usr/local/lib/libdace.a ]'" >> "%LOG%"
set DACE_LIB=%ERRORLEVEL%

if "%DACE_HDR%"=="0" if "%DACE_LIB%"=="0" (
    echo   DACE already installed
    echo.
    goto :JSON
)  
    echo   DACE not found 
    choice /M "Install DACE now?"
    if errorlevel 2 (
        echo Skipped DACE installation by user.
        echo.
        goto :JSON
    )

    REM Remove and clone fresh
    %WSL% "cd /home/%WSLUSER% && git clone https://github.com/dacelib/dace.git dace" >> "%LOG%"

    REM Configure (with algebraic matrix support)
    %WSL% "cmake -S /home/%WSLUSER%/dace -B /home/%WSLUSER%/dace-build -DWITH_ALGEBRAICMATRIX=ON" >> "%LOG%"
    if errorlevel 1 goto :DACE_FAIL

    REM Build
    %WSL% "cmake --build /home/%WSLUSER%/dace-build -j" >> "%LOG%"
    if errorlevel 1 goto :DACE_FAIL

    REM Install into /usr/local
    %WSL% "sudo cmake --install /home/%WSLUSER%/dace-build" >> "%LOG%"
    if errorlevel 1 goto :DACE_FAIL

    %WSL% "sudo ldconfig" >> "%LOG%"

    echo   DACE installed
    echo.
    goto :JSON

:DACE_FAIL
echo [ERROR] DACE installation failed. See %LOG%
pause
exit /b 10

REM ==========================================================
REM == 3) nlohmann-Json detection and installation
REM ==========================================================
:JSON
REM 2026: Code Improved using Claude (Sonnet 4.6)
echo ==== Checking nlohmann-Json ====

REM Header check (allow /usr/include or /usr/local/include)
%WSL% "sh -lc 'test -f /usr/include/nlohmann/json.hpp || test -f /usr/local/include/nlohmann/json.hpp'" >> "%LOG%"
set JSON_HDR=%ERRORLEVEL%

if "%JSON_HDR%"=="0" (
    echo   nlohmann-Json is already installed
    echo.
    goto :DONE
) 
    echo   nlohmann-Json not found 
    choice /M "Install nlohmann-Json now?"
    if errorlevel 2 (
        echo Skipped nlohmann-Json installation by user.
        echo.
        goto :DONE
    )
    
    echo   Installing nlohmann-Json...
    %WSL% "sudo apt-get install -y nlohmann-json3-dev" >> "%LOG%"
    if errorlevel 1 (
        echo [ERROR] Failed to install nlohmann-Json
        echo.
        pause
        exit /b 30
    )
    echo   nlohmann-Json installed 
    echo.

:DONE
echo ============================================
echo   All dependencies checked and installed 
echo   See log: %LOG%
echo ============================================
pause
exit /b 0

