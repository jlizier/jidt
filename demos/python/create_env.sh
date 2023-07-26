#!/bin/bash
##
##  Java Information Dynamics Toolkit (JIDT)
##  Copyright (C) 2022, Joseph T. Lizier
##  
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##

# Create a python environment (stored in folder $folder) with jpype1, numpy and scipy installed

# Name of folder to use and python commands
folder=jpype_env
pythonCmd=python3
pipCmd=pip3

# First make sure that the virtualenv package is installed.
$pipCmd show virtualenv > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "virtualenv already installed, proceeding"
else
    echo "installing virtualenv with $pipCmd ...".
    # On ubuntu, one could also install via the main package manager, e.g. sudo apt-get install python3-venv (I think this takes care of the followng anyway, but am unsure)
    $pythonCmd -m pip install --user virtualenv
    if [ $? -ne 0 ]; then
        echo "pip install of virtualenv failed"
        exit 1
    else
        echo "pip install of virtualenv succeeded"
    fi
fi

# Create a python environment (stored in folder $folder)
$pythonCmd -m venv $folder
if [ $? -ne 0 ]; then
    # Virtual environment creation did not work:
    echo "Virtual environment creation did not work." >&2
    echo "If you are on ubuntu you should now run: sudo apt-get install python3-venv" >&2
    echo "Then run this script again" >&2
    exit 2
else
    echo "Virtual environment created in $folder"
fi

# enter the environment
source $folder/bin/activate
if [ $? -ne 0 ]; then
    echo "Virtual environment unable to be activated" >&2
    exit 3
else
    echo "Python environment started and activated."
    echo "Beginning pip installations for the environment"
fi

# install jpype1 and numpy (does not matter if they are already installed)
$pipCmd install jpype1
$pipCmd install numpy

echo
echo "Jpype1 and numpy installed - you have a functional installation."
echo
echo "Now trying scipy, matplotlib and jupyter, but they are optional..."
echo

$pipCmd install scipy
$pipCmd install matplotlib
$pipCmd install jupyter

echo
echo "Scipy, matplotlib and jupyter installed"
echo

deactivate

