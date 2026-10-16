[CmdletBinding()]
param(
    [switch]$InstallDependencies
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$venv = Join-Path $root ".venv"
$python = Join-Path $venv "Scripts\python.exe"
$checker = Join-Path $root "check_environment.py"
$gui = Join-Path $root "quartz_cherenkov_gui.py"

function Test-VenvPython([string]$PythonPath) {
    if (-not (Test-Path -LiteralPath $PythonPath -PathType Leaf)) { return $false }
    & $PythonPath -c "import numpy, matplotlib, tkinter" *> $null
    return ($LASTEXITCODE -eq 0)
}

if (-not (Test-VenvPython $python)) {
    if (Test-Path -LiteralPath $venv) {
        $stamp = Get-Date -Format "yyyyMMdd-HHmmss"
        $backup = Join-Path $root ".venv.broken-$stamp"
        Move-Item -LiteralPath $venv -Destination $backup
    }
    & py -3 -m venv $venv
    if ($LASTEXITCODE -ne 0) { throw "Could not create the project Python environment." }
}

if ($InstallDependencies) {
    & $python -m pip install -r (Join-Path $root "requirements.txt")
    if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed." }
} elseif (-not (Test-VenvPython $python)) {
    throw "Required packages are missing. Rerun Start_GUI.bat with -InstallDependencies."
}

& $python $checker
if ($LASTEXITCODE -ne 0) { throw "Environment checks failed; the GUI was not started." }

& $python $gui
exit $LASTEXITCODE
