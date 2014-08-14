Java Information Dynamics Toolkit (JIDT)
Copyright (C) 2012 Joseph T. Lizier

Version @VERSION@

JIDT provides a standalone, open source code Java implementation (usable in Matlab, Octave and Python) of information-theoretic measures of distributed computation in complex systems: i.e. information storage, transfer and modification.

This includes implementations for:
- both discrete and continuous-valued variables, principally for the measures transfer entropy, mutual information and active information storage;
- using various types of estimators (e.g. Kraskov-Stögbauer-Grassberger estimators, kernel estimation, linear-Gaussian).
    
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

A research paper describing the toolkit is included in the top level directory -- "InfoDynamicsToolkit.pdf".

Javadocs for the toolkit are included in the full distribution at javadocs.
They can also be generated using "ant javadocs" (useful if you are on an SVN view).
Further, they will soon be posted on the web.

Further documentation is provided by the Usage examples below.

You can also join our email discussion group jidt-discuss at http://groups.google.com/d/forum/jidt-discuss

=============
    Usage
=============

Several sets of demonstration code are distributed with the toolkit:

 a. demos/java -- basic examples on easily using the Java toolkit -- run these from the shell scripts in this directory -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/SimpleJavaExamples

 b. demos/octave -- basic examples on easily using the Java toolkit from Octave or Matlab environments -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/OctaveMatlabExamples
 
 c. demos/python -- basic examples on easily using the Java toolkit from Python -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/PythonExamples

 d. demos/octave/CellularAutomata -- using the Java toolkit to plot local information dynamics profiles in cellular automata; the toolkit is run under Octave or Matlab -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/CellularAutomataDemos
 
 e. demos/octave/SchreiberTransferEntropyExamples -- recreates the transfer entropy examples in Schreiber's original paper presenting this measure; shows the correct parameter settings to reproduce these results  -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/SchreiberTeDemos
 
 f. demos/octave/DetectingInteractionLags -- demonstration of using the transfer entropy with source-destination lags; the demo is run under Octave or Matlab -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/DetectingInteractionLags

 g. demos/java/InterregionalTransfer -- higher level example using collective transfer entropy to infer effective connections between "regions" of data -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/InterregionalTransfer

 h. demos/octave/NullDistributions --  investigating the correspondence between analytic and bootstrapped distributions for TE and MI under null hypotheses of no relationship; the demo is run under Octave or Matlab -- see description at https://code.google.com/p/information-dynamics-toolkit/wiki/NullDistributions

 i. java/unittests -- the JUnit test cases for the Java toolkit are included in the distribution -- these case also be browsed to see simple use cases for the various calculators in the toolkit -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/JUnitTestCases

=============
  Citation
=============

Please cite your use of this toolkit as:

Joseph T. Lizier, "JIDT: An information-theoretic toolkit for studying the dynamics of complex systems", 2012, arXiv:<INSERT_NUMBER>

=============
   Notices
=============

This project includes modified files from the Apache Commons Math library -- http://commons.apache.org/proper/commons-math/
This Apache 2 software is now included as a derivative work in this GPLv3 licensed JIDT project, as per: http://www.apache.org/licenses/GPL-compatibility.html
Notices and license for this software are found in the notices/commons-math directory.

The project includes adapted code from the JAMA project -- http://math.nist.gov/javanumerics/jama/
Notices and license for this software are found in the notices/JAMA directory.

The project includes adapted code from the octave-java package of the Octave-Forge project -- http://octave.sourceforge.net/java/
Notices for this software are found in the notices/JAMA directory.


=============

Joseph T. Lizier, @DATE@

