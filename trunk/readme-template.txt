Java Information Dynamics Toolkit
Copyright (C) 2012 Joseph T. Lizier

Version @VERSION@

=============
   License
=============

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=============
  Citation
=============

Please cite your use of this toolkit as:

Joseph T. Lizier, "JIDT: An information-theoretic toolkit for studying the dynamics of complex systems", 2012, https://code.google.com/p/information-dynamics-toolkit/

=============
   Website
=============

Full information on the JIDT (usage, etc) is provided at the project page on google code:

http://code.google.com/p/information-dynamics-toolkit/

=============
Installation
=============

"Full" description of any required installation is at: http://code.google.com/p/information-dynamics-toolkit/wiki/Installation

However, if you are reading this file, you've downloaded a distribution and you're halfway there!

There are no dependencies to download; unless:
 a. You don't have java installed - download it from http://www.java.com/
 b. You wish to build the project using the build.xml script - this requires ant: http://ant.apache.org/
 c. You wish to run the JUnit test cases - this requires JUnit: http://www.junit.org/ - for how to run JUnit with our ant script see http://code.google.com/p/information-dynamics-toolkit/wiki/JUnitTestCases

Then just put the jar in a relevant location in your file structure.

That's it.

=============
Documentation
=============

Javadocs for the toolkit are included in the full distribution at javadocs.
They can also be generated using "ant javadocs" (useful if you are on an SVN view).
Further, they will soon be posted on the web.

A research paper describing the toolkit and its use is in preparation and will be included in the distribution in future.

Further documentation is provided by the Usage examples below.

=============
    Usage
=============

Several sets of demonstration code are distributed with the toolkit:

 a. demos/octave - basic examples on easily using the Java toolkit from Octave or Matlab environments - see description at http://code.google.com/p/information-dynamics-toolkit/wiki/OctaveMatlabExamples
     Note that this also provides useful examples on how to use the Java code in general!
 
 b. demos/octave/CellularAutomata - using the Java toolkit to plot local information dynamics profiles in cellular automata; the toolkit is run under Octave or Matlab - see description at http://code.google.com/p/information-dynamics-toolkit/wiki/CellularAutomataDemos
 
 c. demos/octave/DetectingInteractionLags - brief examples using the transfer entropy to examine source-delay interaction lags. Documentation to come soon; in the interim, see header comments in the .m files.

 d. java/unittests - the JUnit test cases for the Java toolkit are included in the distribution - these case also be browsed to see simple use cases for the toolkit - see description at http://code.google.com/p/information-dynamics-toolkit/wiki/JUnitTestCases

=============

Joseph T. Lizier, @DATE@

