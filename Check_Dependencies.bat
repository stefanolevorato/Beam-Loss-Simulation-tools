@echo off
setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0check_dependencies.ps1"
exit /b %errorlevel%
