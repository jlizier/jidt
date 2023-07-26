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

# Start the python environment (stored in folder $folder) with jpype1, numpy and others installed
# and then launch jupyter

source start_env.sh

echo "Started virtual environment, now starting jupyter ..."

cd ../..
jupyter notebook

