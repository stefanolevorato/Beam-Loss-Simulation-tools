Windows setup for the complete v72 simulator

Keep every file from v72_gamma_event_list_geometry.zip together in one folder. The
launcher always uses that folder's .venv and runs check_environment.py before
starting the GUI. If an existing .venv Python cannot start, the launcher keeps
the old folder as .venv.broken-<timestamp> and creates a new environment with
the Windows py -3 launcher. It does not remove the old environment.

Start_GUI.bat runs setup_and_start.ps1. Dependency installation is opt-in so
an existing working installation can launch without contacting package
servers. If dependencies need to be installed, open PowerShell in this folder
and run:

  .\setup_and_start.ps1 -InstallDependencies

That command installs requirements.txt, including optional ROOT output
packages uproot and awkward. numpy and matplotlib are required for the GUI.
Tkinter comes with the selected Windows Python installation and is not a pip
package. Python 3 and the Windows py launcher are required.

The full environment/regression checker is a required launch stage. If it
fails, the GUI does not start; read the failed check and resolve it first.
Check_Dependencies.bat runs the same checker without opening the GUI. It does
not run pip's global consistency checker, which can report unrelated packages.

The scripts were checked for structure in the release sandbox but were not
executed on Windows. The simulator's visual GUI also requires an interactive
Windows desktop session.
