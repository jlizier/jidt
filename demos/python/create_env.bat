@ECHO OFF
REM
REM  Java Information Dynamics Toolkit (JIDT)
REM  Copyright (C) 2022, Joseph T. Lizier
REM  
REM  This program is free software: you can redistribute it and/or modify
REM  it under the terms of the GNU General Public License as published by
REM  the Free Software Foundation, either version 3 of the License, or
REM  (at your option) any later version.
REM  
REM  This program is distributed in the hope that it will be useful,
REM  but WITHOUT ANY WARRANTY; without even the implied warranty of
REM  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
REM  GNU General Public License for more details.
REM  
REM  You should have received a copy of the GNU General Public License
REM  along with this program.  If not, see <http://www.gnu.org/licenses/>.
REM

REM Create a python environment (stored in folder %folder%) with jpype1, numpy and scipy installed

REM Name of folder to use and python commands -- change if required:
set folder=jpype_env
set pythonCmd=python
set pipCmd=pip

REM First make sure that the virtualenv package is installed.
%pipCmd% show virtualenv >nul 2>&1
if %errorlevel% == 0 (
    echo virtualenv already installed, proceeding
) else (
    echo installing virtualenv with %pipCmd% ...
    %pythonCmd% -m pip install --user virtualenv
    REM %errorlevel% doesnt seem to return as expect from the above, so checking success via pip: 
    %pipCmd% show virtualenv >nul 2>&1
    if %errorlevel% neq 0 (
        echo pip install of virtualenv failed
        exit /b 1
    ) else (
        echo pip install of virtualenv succeeded
    )
)

REM Create a python environment (stored in folder %folder%)
%pythonCmd% -m venv %folder%
if %errorlevel% neq 0 (
    REM Virtual environment creation did not work:
    echo Virtual environment creation did not work. Do you need to pip install virtualenv? >&2
    exit /b 2
) else (
    echo Virtual environment created in %folder%
)

REM enter the environment
call %folder%\Scripts\activate.bat
if %errorlevel% neq 0 (
    echo Virtual environment unable to be activated
    exit /b 3
) else (
    echo Python environment started and activated.
    echo Beginning pip installations for the environment
)

REM install jpype1 and numpy (does not matter if they are already installed)
%pipCmd% install jpype1
%pipCmd% install numpy

echo.
echo jpype1 and numpy installed - you have a functional installation.
echo.
echo Now trying scipy, matplotlib and jupyter, but they are optional...
echo.

%pipCmd% install scipy
%pipCmd% install matplotlib
%pipCmd% install jupyter

echo.
echo scipy, matplotlib and jupyter installed
echo.

echo.
echo In Powershell activate the environment via calling: %folder%\Scripts\Activate.ps1
echo Otherwise activate the environment via calling: %folder%\Scripts\activate.bat

deactivate

