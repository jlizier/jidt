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

# Create a python environment (stored in folder ) with jpype1 installed
$pythonCmd -m venv $folder

# enter the environment
source $folder/bin/activate

echo "Python environment started and activated"

# install jpype1, numpy and scipy (does not matter if they are already installed)
$pipCmd install jpype1
$pipCmd install numpy

echo
echo "Jpype1 and numpy installed - you have a functional installation."
echo
echo "Now trying scipy matplotlib, but they are optional..."
echo

$pipCmd install scipy
$pipCmd install matplotlib

echo
echo "Scipy matplotlib installed"
echo

deactivate

