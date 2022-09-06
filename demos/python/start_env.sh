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

# Start the python environment (stored in folder $folder) with jpype1, numpy and scipy installed

# Name of folder to use and python commands
folder=jpype_env

# enter the environment
source $folder/bin/activate

echo
echo "Make sure you called this script as: source start_env.sh"
echo
echo "If you called it like that, you will have your python environment activated."
echo "If you just called ./start_env.sh go back and run again as above"
echo

