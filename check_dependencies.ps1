$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$python = Join-Path $root ".venv\Scripts\python.exe"
$checker = Join-Path $root "check_environment.py"

if (-not (Test-Path -LiteralPath $python -PathType Leaf)) {
    Write-Error "Project environment is missing. Run Start_GUI.bat first."
    exit 1
}

& $python $checker
exit $LASTEXITCODE
