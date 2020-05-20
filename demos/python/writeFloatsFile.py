##
##  Java Information Dynamics Toolkit (JIDT)
##  Copyright (C) 2020, Joseph T. Lizier
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

def writeFloatsFile(filename, array):
	"Write a 2D array of floats to a given file"
	with open(filename, "w") as f:
		# Space separate numbers, one time step per line, each column is a variable
		for row in array:
			# write all rows
			first = True;
			for column in row:
				if (not(first)):
					f.write(" ");
				f.write("%f" % column)
				first = False;
			f.write("\n")
		f.close()
	
