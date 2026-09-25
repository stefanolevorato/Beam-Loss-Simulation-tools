@echo off
setlocal EnableExtensions EnableDelayedExpansion

set "PROJECT_DIR=%~dp0"
pushd "%PROJECT_DIR%" >nul
if errorlevel 1 (
    echo ERROR: Could not open the project folder:
    echo        %PROJECT_DIR%
    set "EXIT_CODE=1"
    goto :finish
)

set "VENV_DIR=%PROJECT_DIR%.venv"
for %%I in ("%VENV_DIR%") do set "VENV_DIR=%%~fI"
set "VENV_PYTHON=%VENV_DIR%\Scripts\python.exe"
set "ACTIVATE_BAT=%VENV_DIR%\Scripts\activate.bat"

if not exist "%VENV_PYTHON%" (
    if exist "%VENV_DIR%" (
        echo ERROR: The project .venv folder exists but its Python executable is missing.
        echo        The environment was not deleted or overwritten.
        echo        Repair it manually, or rename/remove "%VENV_DIR%" and run this file again.
        set "EXIT_CODE=1"
        goto :finish
    )

    echo Project .venv not found. Creating it now...
    where py >nul 2>&1
    if not errorlevel 1 (
        py -3 -m venv "%VENV_DIR%"
    ) else (
        where python >nul 2>&1
        if errorlevel 1 (
            echo ERROR: Neither the Python Launcher ^(py^) nor python was found on PATH.
            echo        Install Python for Windows with the Python Launcher, then retry.
            set "EXIT_CODE=1"
            goto :finish
        )
        python -m venv "%VENV_DIR%"
    )
    if errorlevel 1 (
        echo ERROR: Could not create the project virtual environment.
        set "EXIT_CODE=1"
        goto :finish
    )
    if not exist "%VENV_PYTHON%" (
        echo ERROR: The virtual environment was created without Scripts\python.exe.
        echo        Check the Python installation and retry.
        set "EXIT_CODE=1"
        goto :finish
    )
    echo Project .venv created.
)

set "ACTIVE_ENV=__NOT_ACTIVE__"
if defined VIRTUAL_ENV for %%I in ("!VIRTUAL_ENV!") do set "ACTIVE_ENV=%%~fI"
if /I "!ACTIVE_ENV!"=="!VENV_DIR!" (
    echo Project .venv already active.
) else (
    echo Activating project .venv...
    if not exist "%ACTIVATE_BAT%" (
        echo ERROR: The project activation script is missing:
        echo        %ACTIVATE_BAT%
        set "EXIT_CODE=1"
        goto :finish
    )
    call "%ACTIVATE_BAT%"
    if errorlevel 1 (
        echo ERROR: Could not activate the project virtual environment.
        set "EXIT_CODE=1"
        goto :finish
    )
)

if not exist "%PROJECT_DIR%requirements.txt" (
    echo ERROR: Missing requirements.txt in the project folder.
    set "EXIT_CODE=1"
    goto :finish
)

"%VENV_PYTHON%" -c "import numpy, matplotlib" >nul 2>&1
if errorlevel 1 (
    echo Required packages are missing. Installing requirements.txt...
    "%VENV_PYTHON%" -m pip install -r "%PROJECT_DIR%requirements.txt"
    if errorlevel 1 (
        echo ERROR: Dependency installation failed.
        set "EXIT_CODE=1"
        goto :finish
    )
    "%VENV_PYTHON%" -c "import numpy, matplotlib" >nul 2>&1
    if errorlevel 1 (
        echo ERROR: numpy and matplotlib are still unavailable after installation.
        set "EXIT_CODE=1"
        goto :finish
    )
)

echo Starting quartz_cherenkov_gui.py...
"%VENV_PYTHON%" "%PROJECT_DIR%quartz_cherenkov_gui.py"
set "EXIT_CODE=!ERRORLEVEL!"
if not "!EXIT_CODE!"=="0" echo GUI exited with code !EXIT_CODE!.

:finish
if not defined EXIT_CODE set "EXIT_CODE=0"
popd
echo.
if "!EXIT_CODE!"=="0" (
    echo GUI finished with exit code 0.
) else (
    echo Start_GUI.bat finished with error code !EXIT_CODE!.
)
pause
endlocal & exit /b %EXIT_CODE%
