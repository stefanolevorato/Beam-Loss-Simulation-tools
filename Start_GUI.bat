@echo off
setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0setup_and_start.ps1" %*
exit /b %errorlevel%
