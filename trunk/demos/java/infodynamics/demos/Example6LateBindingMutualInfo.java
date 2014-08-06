/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.demos;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.ParsedProperties;

/**
 * = Example 6 - Late binding Mutual info calculator =
 *
 * This class is used to demonstrate the manner in which a user
 *  can code to the interfaces defined in infodynamics.measures.continuous,
 *  and dynamically alter the instantiated class at runtime.
 * We demonstrate this using a multivariate mutual information calculation.
 * 
 * This example also demonstrates how to read simple files of arrays of data
 *  with the toolkit, as well as how to dynamically load properties from a 
 *  java properties file. 
 * 
 * @author Joseph Lizier
 *
 */
public class Example6LateBindingMutualInfo {

	/**
	 * @param args One command line argument taken, specifying location of 
	 *  the properties file. This should be example6LateBindingMutualInfo.props
	 *  in the demos/java directory.
	 */
	public static void main(String[] args) throws Exception {
		
		// 0. Preliminaries (reading in the dynamic properties and the data):
		//  a. Read in the properties file defined as the first
		//     command line argument:
		ParsedProperties props = new ParsedProperties(args[0]);
		//  b. Read in the data file, whose filename is defined in the
		//     property "datafile" in our properties file:
		ArrayFileReader afr = new ArrayFileReader(props.getStringProperty("datafile"));
		double[][] data = afr.getDouble2DMatrix();
		//  c. Pull out the columns from the data set which 
		//     correspond to each of variable 1 and 2: 
		int[] variable1Columns = props.getIntArrayProperty("variable1Columns");
		int[] variable2Columns = props.getIntArrayProperty("variable2Columns");
		double[][] variable1 = MatrixUtils.selectColumns(data, variable1Columns);
		double[][] variable2 = MatrixUtils.selectColumns(data, variable2Columns);
		
		// 1. Create a reference for our calculator as
		//  an object implementing the interface type:
		MutualInfoCalculatorMultiVariate miCalc;
		
		// 2. Define the name of the class to be instantiated here:
		String implementingClass = props.getStringProperty("implementingClass");
		
		// 3. Dynamically instantiate an object of the given class:
		//  Part 1: Class.forName(implementingClass) grabs a reference to
		//   the class named by implementingClass.
		//  Part 2: .newInstance() creates an object instance of that class.
		//  Part 3: (MutualInfoCalculatorMultiVariate) casts the return
		//   object into an instance of our generic interface type.
		miCalc = (MutualInfoCalculatorMultiVariate)
				Class.forName(implementingClass).newInstance();
		
		// 4. Start using our MI calculator, paying attention to only
		//  call common methods defined in the interface type, not methods
		//  only defined in a given implementation class.
		// a. Initialise the calculator to use the required number of
		//   dimensions for each variable:
		miCalc.initialise(variable1Columns.length, variable2Columns.length);
		// b. Supply the observations to compute the PDFs from:
		miCalc.setObservations(variable1, variable2);
		// c. Make the MI calculation:
		double miValue = miCalc.computeAverageLocalOfObservations();
		
		System.out.printf("MI calculator %s computed the joint MI as %.5f\n",
				implementingClass, miValue);
	}

}
