# Avvia Quartz Cherenkov usando il virtual environment locale

$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
Set-Location $ProjectDir

$ActivateScript = Join-Path $ProjectDir ".venv\Scripts\Activate.ps1"
$PythonScript = Join-Path $ProjectDir "quartz_cherenkov_gui.py"

if (-not (Test-Path $ActivateScript)) {
    Write-Host "ERRORE: virtual environment non trovato:" -ForegroundColor Red
    Write-Host $ActivateScript
    Read-Host "Premi Invio per uscire"
    exit 1
}

if (-not (Test-Path $PythonScript)) {
    Write-Host "ERRORE: file quartz_cherenkov_gui.py non trovato:" -ForegroundColor Red
    Write-Host $PythonScript
    Read-Host "Premi Invio per uscire"
    exit 1
}

# Attiva il virtual environment
. $ActivateScript

# Usa il backend Matplotlib compatibile con la policy Windows
$env:MPLBACKEND = "Agg"

Write-Host "Virtual environment attivo:" -ForegroundColor Green
Write-Host $env:VIRTUAL_ENV
Write-Host "Avvio Quartz Cherenkov GUI..." -ForegroundColor Green

python $PythonScript

$ExitCode = $LASTEXITCODE

if ($ExitCode -ne 0) {
    Write-Host ""
    Write-Host "Il programma è terminato con codice: $ExitCode" -ForegroundColor Red
    Read-Host "Premi Invio per chiudere"
}

exit $ExitCode