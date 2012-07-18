package infodynamics.utils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Vector;

/**
 * Utilities for computations on matrices.
 * All multidimensional matrices are assumed to have consistent lengths in each dimension
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 */
public class MatrixUtils {

	public static double sum(double[] input) {
		double total = 0;
		for (int i = 0; i < input.length; i++) {
			total += input[i];
		}
		return total;
	}
	
	public static double sum(double[] input, int startIndex, int length) {
		double total = 0;
		for (int i = startIndex; i < startIndex + length; i++) {
			total += input[i];
		}
		return total;
	}
	
	public static double sumSpecificIndices(double[] input, int[] indices) {
		double total = 0;
		for (int i = 0; i < indices.length; i++) {
			total += input[indices[i]];
		}
		return total;
	}

	public static double sumSpecificIndices(double[] input, int[][] indices, int columnInIndices) {
		double total = 0;
		for (int i = 0; i < indices.length; i++) {
			total += input[indices[i][columnInIndices]];
		}
		return total;
	}

	public static double sumSpecificIndices(double[] input, int[][] indices, int columnInIndices,
			int indicesOffset) {
		double total = 0;
		for (int i = 0; i < indices.length; i++) {
			total += input[indices[i][columnInIndices] + indicesOffset];
		}
		return total;
	}

	public static double sum(double[][] input) {
		double total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[i].length; j++) {
				total += input[i][j];
			}
		}
		return total;
	}

	public static double sum(double[][] input, int column) {
		double total = 0;
		for (int i = 0; i < input.length; i++) {
			total += input[i][column];
		}
		return total;
	}

	public static int sum(int[] input) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			total += input[i];
		}
		return total;
	}

	public static int sum(int[][] input) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[i].length; j++) {
				total += input[i][j];
			}
		}
		return total;
	}

	/**
	 * Return an array of the sums for each column in the 2D input
	 * 
	 * @param input
	 * @return
	 */
	public static double[] sums(double[][] input) {
		double[] theSums = new double[input[0].length];
		for (int r = 0; r < input.length; r++) {
			for (int c = 0; c < input[r].length; c++) {
				theSums[c] += input[r][c];
			}
		}
		return theSums;
	}
	
	public static int countIf(int[] input, int condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			if (input[i] == condition)
				total++;
		}
		return total;
	}

	public static int countIf(int[][] input, int condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[0].length; j++) {
				if (input[i][j] == condition)
					total++;				
			}
		}
		return total;
	}

	public static int countIf(long[][] input, long condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[0].length; j++) {
				if (input[i][j] == condition)
					total++;				
			}
		}
		return total;
	}

	public static int countIf(int[] input1, int condition1, int[] input2, int condition2) 
		throws Exception {
		
		if (input1.length != input2.length)
			throw new Exception("MatrixUtils.sumIf() - arguments are not of equal length (" +
					input1.length + " != " + input2.length + ")");
		int total = 0;
		for (int i = 0; i < input1.length; i++) {
			if ((input1[i] == condition1) && (input2[i] == condition2))
				total++;
		}
		return total;
	}

	public static int countIf(int[] input1, int condition1, int[] input2, int condition2,
			int[] input3, int condition3) throws Exception {
	
		if ((input1.length != input2.length) || (input1.length != input3.length))
			throw new Exception("MatrixUtils.sumIf() - arguments are not of equal length (" +
					input1.length + " != " + input2.length + " != " + input3.length + ")");
		int total = 0;
		for (int i = 0; i < input1.length; i++) {
			if ((input1[i] == condition1) && (input2[i] == condition2) && (input3[i] == condition3))
				total++;
		}
		return total;
	}
	
	public static int countIf(boolean[] input, boolean condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			if (input[i] == condition)
				total++;
		}
		return total;
	}

	public static double mean(int[] input) {
		return sum(input) / (double) input.length;
	}
	
	public static double mean(double[] input) {
		return sum(input) / (double) input.length;
	}
	
	public static double mean(double[] input, int startIndex, int length) {
		return sum(input, startIndex, length) / (double) length;
	}

	public static double mean(double[][] input) {
		return sum(input) / (double) (input.length * input[0].length);
	}

	/**
	 * Compute the mean along the given column 
	 * 
	 * @param input
	 * @param column
	 * @return
	 */
	public static double mean(double[][] input, int column) {
		return sum(input, column) / (double) input.length;
	}

	/**
	 * Return an array of the means of each column in the 2D input
	 * 
	 * @param input
	 * @return
	 */
	public static double[] means(double[][] input) {
		double[] theMeans = sums(input);
		for (int i = 0; i < theMeans.length; i++) {
			theMeans[i] = theMeans[i] / input.length;
		}
		return theMeans;
	}

	/**
	 * Return an array of the means of each row in the 2D input matrix
	 * 
	 * @param input
	 * @return
	 */
	public static double[] meansOfRows(double[][] input) {
		double[] theMeans = new double[input.length];
		for (int i = 0; i < input.length; i++) {
			theMeans[i] = mean(input[i]);
		}
		return theMeans;
	}

	public static int[][] columnShift(int[][] input, int shiftBy){
		if (shiftBy == 0) {
			return input;
		}
		
		int rows = input.length;
		int columns = input[0].length;
		for ( ; shiftBy < 0; shiftBy += columns) {
			// Using % mod operator to come back to a +ve column value wont work.
			// So we're shifting the shiftBy value (above) until it's in the appropriate
			// range 0 .. columns-1
		}
		int[][] output = new int[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				output[r][(c + shiftBy) % columns] = input[r][c];
			}
		}
		return output;
	}
	
	public static double[][] columnShift(double[][] input, int shiftBy){
		if (shiftBy == 0) {
			return input;
		}
		
		int rows = input.length;
		int columns = input[0].length;
		for ( ; shiftBy < 0; shiftBy += columns) {
			// Using % mod operator to come back to a +ve column value wont work.
			// So we're shifting the shiftBy value (above) until it's in the appropriate
			// range 0 .. columns-1
		}
		double[][] output = new double[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				output[r][(c + shiftBy) % columns] = input[r][c];
			}
		}
		return output;
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows or columns. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @param fromColumn
	 * @param colums
	 * @return Single dimensional array containing the required data
	 */
	public static int[] matrixToArray(int[][] input, int fromRow, int rows, int fromColumn, int columns) {
		int[] output = new int[rows * columns];
		for (int c = 0; c < columns; c++) {
			for (int r = 0; r < rows; r++) {
				output[c * rows + r] = input[r + fromRow][c + fromColumn];
			}
		}
		return output;
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows only. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @return Single dimensional array containing the required data
	 */
	public static int[] matrixToArray(int[][] input, int fromRow, int rows) {
		return matrixToArray(input, fromRow, rows, 0, input[0].length);
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * 
	 * @param input
	 * @return Single dimensional array containing the required data
	 */
	public static int[] matrixToArray(int[][] input) {
		return matrixToArray(input, 0, input.length, 0, input[0].length);
	}
	
	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows or columns. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @param fromColumn
	 * @param colums
	 * @return Single dimensional array containing the required data
	 */
	public static double[] matrixToArray(double[][] input, int fromRow, int rows, int fromColumn, int columns) {
		double[] output = new double[rows * columns];
		for (int c = 0; c < columns; c++) {
			for (int r = 0; r < rows; r++) {
				output[c * rows + r] = input[r + fromRow][c + fromColumn];
			}
		}
		return output;
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows only. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @return Single dimensional array containing the required data
	 */
	public static double[] matrixToArray(double[][] input, int fromRow, int rows) {
		return matrixToArray(input, fromRow, rows, 0, input[0].length);
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * 
	 * @param input
	 * @return Single dimensional array containing the required data
	 */
	public static double[] matrixToArray(double[][] input) {
		return matrixToArray(input, 0, input.length, 0, input[0].length);
	}

	/**
	 * Adds two arrays together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 */
	public static int[] add(int[] input1, int[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		int[] returnValues = new int[input1.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = input1[i] + input2[i];
		}
		return returnValues;
	}

	/**
	 * Adds two arrays together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 */
	public static double[] add(double[] input1, double[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		double[] returnValues = new double[input1.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = input1[i] + input2[i];
		}
		return returnValues;
	}

	/**
	 * Adds two arrays together, returning the result in input1
	 * 
	 * @param input1
	 * @param input2
	 */
	public static void addInPlace(double[] input1, double[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		for (int i = 0; i < input1.length; i++) {
			input1[i] = input1[i] + input2[i];
		}
	}

	/**
	 * Adds the squares of the second array to the first,
	 *  returning the result in input1
	 * 
	 * @param input1
	 * @param input2
	 */
	public static void addSquaresInPlace(double[] input1, double[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		for (int i = 0; i < input1.length; i++) {
			input1[i] = input1[i] + input2[i] * input2[i];
		}
	}

	/**
	 * Adds two matrices together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 * @throws Exception
	 */
	public static int[][] add(int[][] input1, int[][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		int[][] returnValues = new int[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				returnValues[r][c] = input1[r][c] + input2[r][c];
			}
		}
		return returnValues;
	}
	
	/**
	 * Adds two matrices together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 * @throws Exception
	 */
	public static double[][] add(double[][] input1, double[][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		double[][] returnValues = new double[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				returnValues[r][c] = input1[r][c] + input2[r][c];
			}
		}
		return returnValues;
	}

	/**
	 * Adds two matrices together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 * @throws Exception
	 */
	public static double[][][] add(double[][][] input1, double[][][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		int height = input1[0][0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		if (input2[0][0].length != height) {
			throw new Exception("Heights (3rd dim) of arrays are not equal");
		}
		double[][][] returnValues = new double[rows][columns][height];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				for (int h = 0; h < height; h++) {
					returnValues[r][c][h] = input1[r][c][h] + input2[r][c][h];
				}
			}
		}
		return returnValues;
	}

	/**
	 * Subtracts second array from the first
	 * 
	 * @param first
	 * @param second
	 * @return first - second
	 */
	public static double[] subtract(double[] first, double[] second) throws Exception {
		if (first.length != second.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		double[] returnValues = new double[first.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = first[i] - second[i];
		}
		return returnValues;
	}

	/**
	 * Subtracts second array from the first, overwriting the
	 *  values in first
	 * 
	 * @param first
	 * @param second
	 */
	public static void subtractInPlace(double[] first, double[] second) throws Exception {
		if (first.length != second.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		for (int i = 0; i < first.length; i++) {
			first[i] = first[i] - second[i];
		}
	}

	/**
	 * Subtract one matrix from another
	 * 
	 * @param input1
	 * @param input2
	 * @return input1 - input2
	 * @throws Exception
	 */
	public static int[][] subtract(int[][] input1, int[][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		int[][] returnValues = new int[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				returnValues[r][c] = input1[r][c] - input2[r][c];
			}
		}
		return returnValues;
	}

	/**
	 * Return the matrix product a x b
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static double[][] matrixProduct(double[][] a, double[][] b) throws Exception {
		if (a[0].length != b.length) {
			throw new Exception("Number of columns of a " + a[0].length +
					" does not match the number of rows of b " + b.length);
		}
		double[][] result = new double[a.length][b[0].length];
		for (int r = 0; r < result.length; r++) {
			for (int c = 0; c < result[r].length; c++) {
				result[r][c] = 0;
				for (int k = 0; k < a[r].length; k++) {
					result[r][c] += a[r][k] * b[k][c];
				}
			}
		}
		return result;
	}
	
	/**
	 * Duplicates a matrix; handles different number of columns
	 *  for each row
	 * 
	 * @param src
	 * @return
	 */
	public static int[][] duplicateMatrix(int[][] src) {
		int[][] dest = new int[src.length][];
		for (int r = 0; r < src.length; r++) {
			dest[r] = new int[src[r].length];
			System.arraycopy(src[r], 0, dest[r], 0, src[r].length);
		}
		return dest;
	}
	
	/**
	 * Copies all rows and columns between two double arrays
	 * 
	 * @param src
	 * @param dest
	 */
	public static void arrayCopy(double[][] src, double[][] dest) {
		for (int r = 0; r < src.length; r++) {
			System.arraycopy(src[r], 0,	dest[r], 0,	src[r].length);
		}
	}

	/**
	 * Copies the required rows and columns between two 
	 * double arrays
	 * 
	 * @param src
	 * @param srcStartRow
	 * @param srcStartCol
	 * @param dest
	 * @param destStartRow
	 * @param destStartCol
	 * @param rows
	 * @param cols
	 */
	public static void arrayCopy(double[][] src, int srcStartRow, int srcStartCol,
								double[][] dest, int destStartRow, int destStartCol,
								int rows, int cols) {
		
		for (int r = 0; r < rows; r++) {
			System.arraycopy(src[srcStartRow + r], srcStartCol,
					dest[destStartRow + r], destStartCol,
					cols);
		}
	}

	/**
	 * Copies all rows and columns between two int arrays
	 * 
	 * @param src
	 * @param dest
	 */
	public static void arrayCopy(int[][] src, int[][] dest) {
		for (int r = 0; r < src.length; r++) {
			System.arraycopy(src[r], 0,	dest[r], 0,	src[r].length);
		}
	}

	/**
	 * Copies the required rows and columns between two 
	 * double arrays
	 * 
	 * @param src
	 * @param srcStartRow
	 * @param srcStartCol
	 * @param dest
	 * @param destStartRow
	 * @param destStartCol
	 * @param rows
	 * @param cols
	 */
	public static void arrayCopy(int[][] src, int srcStartRow, int srcStartCol,
								int[][] dest, int destStartRow, int destStartCol,
								int rows, int cols) {
		
		for (int r = 0; r < rows; r++) {
			System.arraycopy(src[srcStartRow + r], srcStartCol,
					dest[destStartRow + r], destStartCol,
					cols);
		}
	}
	
	/**
	 * Copies the given source array into the required column number of the destination
	 * @param destination
	 * @param column
	 * @param source
	 */
	public static void copyIntoColumn(int[][] destination, int column, int[] source) throws Exception {
		if (source.length != destination.length) {
			throw new Exception("Destination column is not of the same length as the source (" +
					destination.length + " vs " + source.length + ")");
		}
		for (int r = 0; r < destination.length; r++) {
			destination[r][column] = source[r];
		}
	}

	/**
	 * Copies the given source array into the required column number of the destination
	 * @param destination
	 * @param column
	 * @param source
	 */
	public static void copyIntoColumn(double[][] destination, int column, 
			int destFromRowNumber, double[] source, int sourceFromRowNumber,
			int rows) throws Exception {
		if (sourceFromRowNumber + rows > source.length) {
			throw new Exception("Attempting to copy too many rows " + rows +
					" after the start row " + sourceFromRowNumber +
					" from the source of length " + source.length);
		}
		if (destFromRowNumber + rows > destination.length) {
			throw new Exception("Attempting to copy too many rows " + rows +
					" after the start row " + destFromRowNumber +
					" from the destination of length " + destination.length);
		}
		for (int r = 0; r < rows; r++) {
			destination[r + destFromRowNumber][column] = source[r + sourceFromRowNumber];
		}
	}

	/**
	 * Copies the given source array into the required column number of the destination
	 * @param destination
	 * @param column
	 * @param source
	 */
	public static void copyIntoColumn(double[][] destination, int column, double[] source) throws Exception {
		if (source.length != destination.length) {
			throw new Exception("Destination column is not of the same length as the source (" +
					destination.length + " vs " + source.length + ")");
		}
		for (int r = 0; r < destination.length; r++) {
			destination[r][column] = source[r];
		}
	}

	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply each column
	 *  by a different power of the base.
	 */
	public static int[] computeCombinedValues(int separateValues[][], int base) throws Exception {
		// Number of columns (second index) is sizeof first element
		int columns = separateValues[0].length;
		
		return computeCombinedValues(separateValues, columns, base);
		/*
		// rows = first index specifies rows
		int rows = separateValues.length;
		int[] combinedValues = new int[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			int combinedRowValue = 0;
			int multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += separateValues[r][c] * multiplier;
				multiplier *= base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
		*/
	}
	
	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply the first "columbs" columns
	 *  by a different power of the base.
	 */
	public static int[] computeCombinedValues(int separateValues[][], int columns, int base) throws Exception {
		if (columns > separateValues[0].length) {
			throw new Exception("computeCombinedValues: computation request across more columns " +
					columns + " than are available " + separateValues[0].length);
		}
		// Make sure we won't get any overflow here
		if (combinedValuesOverflow(columns, base)) {
			// multiplier has overflown
			throw new Exception("Too many columns " + columns + " for the given base " + base +
					" for this call to computeCombinedValues");
		}
		
		// rows = first index specifies rows
		int rows = separateValues.length;
		int[] combinedValues = new int[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			int combinedRowValue = 0;
			int multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += separateValues[r][c] * multiplier;
				multiplier *= base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
	}

	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply each column
	 *  by a different power of the base.
	 */
	public static long[] computeCombinedValuesLong(int separateValues[][], int base) throws Exception {
		// Number of columns (second index) is sizeof first element
		int columns = separateValues[0].length;
		
		return computeCombinedValuesLong(separateValues, columns, base);
		/*
		// rows = first index specifies rows
		int rows = separateValues.length;
		long[] combinedValues = new long[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			long combinedRowValue = 0;
			long multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += ((long) separateValues[r][c]) * multiplier;
				multiplier *= (long) base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
		*/
	}

	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply the first "columns" columns
	 *  by a different power of the base.
	 */
	public static long[] computeCombinedValuesLong(int separateValues[][], int columns, int base) throws Exception {
		if (columns > separateValues[0].length) {
			throw new Exception("computeCombinedValuesLong: computation request across more columns " +
					columns + " than are available " + separateValues[0].length);
		}
		
		// Make sure we won't get any overflow here
		if (combinedValuesOverflowLong(columns, base)) {
			// multiplier has overflown
			throw new Exception("Too many columns " + columns + " for the given base " + base +
					" for this call to computeCombinedValuesLong");
		}
		
		// rows = first index specifies rows
		int rows = separateValues.length;
		long[] combinedValues = new long[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			long combinedRowValue = 0;
			long multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += ((long) separateValues[r][c]) * multiplier;
				multiplier *= (long) base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
	}

	public static boolean combinedValuesOverflow(int columns, int base) {
		// Make sure we won't get any overflow here
		int multiplier = 1;
		for (int c = columns - 1; c >= 0; c--) {
			if (multiplier < 0) {
				// multiplier has overflown.
				// Technically, it's possible to use one negative value if we were using base-2,
				//  but realistically it's safer if we just call it off now.
				return true;
			}
			multiplier *= (long) base;
		}
		return false;
	}

	public static boolean combinedValuesOverflowLong(int columns, int base) {
		// Make sure we won't get any overflow here
		long multiplier = 1;
		for (int c = columns - 1; c >= 0; c--) {
			if (multiplier < 0) {
				// multiplier has overflown.
				// Technically, it's possible to use one negative value if we were using base-2,
				//  but realistically it's safer if we just call it off now.
				return true;
			}
			multiplier *= (long) base;
		}
		return false;
	}
	
	/**
	 * Select out part of an array.
	 * 
	 * @param data
	 * @param fromIndex
	 * @param length
	 * @return
	 */
	public static double[] select(double[] data, int fromIndex, int length) {
		double[] returnData = new double[length];
		System.arraycopy(data, fromIndex, returnData, 0, length);
		return returnData;
	}
	
	/**
	 * Select out part of an array.
	 * 
	 * @param data
	 * @param fromIndex
	 * @param length
	 * @return
	 */
	public static int[] select(int[] data, int fromIndex, int length) {
		int[] returnData = new int[length];
		System.arraycopy(data, fromIndex, returnData, 0, length);
		return returnData;
	}
	
	public static int[] selectColumn(int matrix[][], int columnNo) {
		int[] column = new int[matrix.length];
		for (int r = 0; r < matrix.length; r++) {
			column[r] = matrix[r][columnNo];
		}
		return column;
	}
	
	public static double[] selectColumn(double matrix[][], int columnNo) {
		double[] column = new double[matrix.length];
		for (int r = 0; r < matrix.length; r++) {
			column[r] = matrix[r][columnNo];
		}
		return column;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param columns
	 * @return
	 */
	public static double[][] selectColumns(double matrix[][], int columns[]) {
		double[][] data = new double[matrix.length][columns.length];
		for (int r = 0; r < matrix.length; r++) {
			for (int cIndex = 0; cIndex < columns.length; cIndex++) {
				data[r][cIndex] = matrix[r][columns[cIndex]];
			}
		}
		return data;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param includeColumnFlags
	 * @return
	 */
	public static double[][] selectColumns(double matrix[][], boolean includeColumnFlags[]) {
		Vector<Integer> v = new Vector<Integer>();
		
		for (int i = 0; i < includeColumnFlags.length; i++) {
			if (includeColumnFlags[i]) {
				v.add(new Integer(i));
			}
		}
		double[][] data = new double[matrix.length][v.size()];
		for (int r = 0; r < matrix.length; r++) {
			for (int outputColumnIndex = 0; outputColumnIndex < v.size(); outputColumnIndex++) {
				int outputColumn = v.get(outputColumnIndex);
				data[r][outputColumnIndex] = matrix[r][outputColumn];
			}
		}
		return data;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param columns
	 * @return
	 */
	public static double[][] selectColumns(double matrix[][], Vector<Integer> columns) {
		double[][] data = new double[matrix.length][columns.size()];
		for (int r = 0; r < matrix.length; r++) {
			for (int cIndex = 0; cIndex < columns.size(); cIndex++) {
				data[r][cIndex] = matrix[r][columns.elementAt(cIndex).intValue()];
			}
		}
		return data;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param columns
	 * @return
	 */
	public static double[][] selectRows(double matrix[][], int fromRow, int rows) {
		double[][] data = new double[rows][];
		for (int rIndex = 0; rIndex < rows; rIndex++) {
			data[rIndex] = matrix[rIndex + fromRow];
		}
		return data;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param columns
	 * @return
	 */
	public static double[][] selectRowsAndColumns(double matrix[][], int rows[], int columns[]) {
		double[][] data = new double[rows.length][columns.length];
		for (int rIndex = 0; rIndex < rows.length; rIndex++) {
			for (int cIndex = 0; cIndex < columns.length; cIndex++) {
				data[rIndex][cIndex] = matrix[rows[rIndex]][columns[cIndex]];
			}
		}
		return data;
	}

	public static double[][] selectFirstTwoDimenions(double[][][][] matrix, int d2, int d3) {
		double[][] newMatrix = new double[matrix.length][];
		for (int i = 0; i < matrix.length; i++) {
			newMatrix[i] = new double[matrix[i].length];
			for (int j = 0; j < matrix[i].length; j++) {
				newMatrix[i][j] = matrix[i][j][d2][d3];
			}
		}
		return newMatrix;
	}
	
	public static double[][] copyMatrixEliminateRowAndColumn(double[][] matrix,
			int rowToEliminate, int colToEliminate) {
		double[][] newMatrix = new double[matrix.length - 1][matrix[0].length - 1];
		for (int r = 0; r < matrix.length; r++) {
			if (r == rowToEliminate) {
				continue;
			}
			for (int c = 0; c < matrix.length; c++) {
				if (c == colToEliminate) {
					continue;
				}
				int newRow = r;
				int newCol = c;
				if (newRow > rowToEliminate) {
					newRow--;
				}
				if (newCol > colToEliminate) {
					newCol--;
				}
				newMatrix[newRow][newCol] = matrix[r][c];
			}
		}
		return newMatrix;
	}
	
	public static double[] extractSelectedTimePoints(double[] data, int[] timePoints) {
		double[] extracted = new double[timePoints.length];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	public static int[] extractSelectedTimePoints(int[] data, int[] timePoints) {
		int[] extracted = new int[timePoints.length];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	public static double[][] extractSelectedTimePoints(double[][] data, int[] timePoints) {
		int columns = data[0].length;
		double[][] extracted = new double[timePoints.length][columns];
		for (int t = 0; t < timePoints.length; t++) {
			System.arraycopy(data[timePoints[t]], 0, extracted[t], 0, columns);
		}
		return extracted;
	}

	/**
	 * Extraxts the double[] vectors at each of the selected time points.
	 * The return double[][] array is an array of points to the existing
	 * double[] vectors.
	 * 
	 * @param data
	 * @param timePoints
	 * @return
	 */
	public static double[][] extractSelectedTimePointsReusingArrays(double[][] data, int[] timePoints) {
		double[][] extracted = new double[timePoints.length][];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	/**
	 * Extraxts the boolean[] vectors at each of the selected time points.
	 * The return boolean[][] array is an array of points to the existing
	 * boolean[] vectors.
	 * 
	 * @param data
	 * @param timePoints
	 * @return
	 */
	public static boolean[][] extractSelectedTimePointsReusingArrays(boolean[][] data, int[] timePoints) {
		boolean[][] extracted = new boolean[timePoints.length][];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	public static double[][] extractSelectedTimePoints(double[][] data, int[][] timePoints,
			int columnInTimePoints) {
		int columns = data[0].length;
		double[][] extracted = new double[timePoints.length][columns];
		for (int t = 0; t < timePoints.length; t++) {
			System.arraycopy(data[timePoints[t][columnInTimePoints]], 0, extracted[t], 0, columns);
		}
		return extracted;
	}

	/**
	 * Extract from data the vectors for rows corresponding to the time values in
	 *  column columnInTimePoints of each row of timePoints.
	 * 
	 * @param data
	 * @param timePoints
	 * @param columnInTimePoints
	 * @param timeOffset
	 * @return a 2D array of doubles, with timePoints.length rows and data[0].length columns
	 */
	public static double[][] extractSelectedTimePoints(double[][] data, int[][] timePoints,
			int columnInTimePoints, int timeOffset) {
		int columns = data[0].length;
		double[][] extracted = new double[timePoints.length][columns];
		for (int t = 0; t < timePoints.length; t++) {
			System.arraycopy(data[timePoints[t][columnInTimePoints] + timeOffset], 0, extracted[t], 0, columns);
		}
		return extracted;
	}
	
	/**
	 * Return the rows of data, where the conditionalData matched the
	 *  conditionalValue for that given row.
	 * Assumes data.length == conditionalData.length.
	 * 
	 * @param data
	 * @param conditionalData
	 * @param conditionalValue
	 * @return a 2D array of doubles where the conditionalData matched the
	 *  conditionalValue for those rows.
	 */
	public static double[][] extractSelectedPointsMatchingCondition(
			double[][] data, int[] conditionalData, int conditionalValue) {
		
		// Count the number of matching points first.
		int numNewRows = 0;
		for (int t = 0; t < data.length; t++) {
			if (conditionalData[t] == conditionalValue) {
				numNewRows++;
			}
		}
		// Create the new extracted data
		return extractSelectedPointsMatchingCondition(data, conditionalData,
				conditionalValue, numNewRows);
	}

	/**
	 * Return the rows of data, where the conditionalData matched the
	 *  conditionalValue for that given row.
	 * Assumes data.length == conditionalData.length.
	 * Here, the caller knows that there will be at minimum knownNumExtractedValues
	 *  values to be extracted, and only wants those values.
	 * 
	 * @param data
	 * @param conditionalData
	 * @param conditionalValue
	 * @param knownNumExtractedValues the known number of matching values
	 * @return a 2D array of doubles where the conditionalData matched the
	 *  conditionalValue for those rows.
	 */
	public static double[][] extractSelectedPointsMatchingCondition(
			double[][] data, int[] conditionalData, int conditionalValue,
			int knownNumExtractedValues) {
		
		// Create the new extracted data
		int columns = data[0].length;
		double[][] extracted = new double[knownNumExtractedValues][columns];
		int rowsCopied = 0;
		if (knownNumExtractedValues == 0) {
			return extracted;
		}
		for (int t = 0; t < data.length; t++) {
			if (conditionalData[t] == conditionalValue) {
				System.arraycopy(data[t], 0, extracted[rowsCopied++],
						0, columns);
			}
			if (rowsCopied == knownNumExtractedValues) {
				// We've extracted enough values
				break;
			}
		}
		return extracted;
	}

	/**
	 * Inserts the given time points (in the order prescribed in timePoints) 
	 * from the vector originalSourceValuesInJoint into the given column in matrix
	 * 
	 * @param inputValues
	 * @param timePoints
	 * @param matrix
	 * @param column
	 */
	public static void reorderVectorIntoMatrix(double[] inputValues, int[] timePoints,
			double[][] matrix, int column) {
		for (int i = 0; i < timePoints.length; i++) {
			int t = timePoints[i];
			matrix[i][column] = inputValues[t];
		}
	}
	
	/**
	 * Return data[x][y]:
	 *  - y==0: inputValues[x][0]
	 *  - y>0:  inputValues[reordering[y-1][x]][y]
	 * 
	 * @param inputValues holds the raw data values
	 * @param reordering outlines how to rearrange the raw data values for each variable or column.
	 *  First index is variable
	 *  or column number. Reorderings may be supplied for all of the columns of the inputValues,
	 *  or for one less than all of the columns, in which case the first column is not
	 *  reordered. Second index is for the row number or time step. The value at that
	 *  point states which row number to pull the data from.
	 * @return
	 */
	public static double[][] reorderDataForVariables(double[][] inputValues, int[][] reordering) {
		int rows = inputValues.length;
		int columns = inputValues[0].length;
		boolean reorderingFirstColumn = (reordering.length == columns);
		double[][] data = new double[rows][columns];
		for (int r = 0; r < rows; r++) {
			int reorderIndex = 0;
			if (reorderingFirstColumn) {
				data[r][0] = inputValues[reordering[reorderIndex++][r]][0];
			} else {
				data[r][0] = inputValues[r][0];
			}
			for (int c = 1; c < columns; c++) {
				data[r][c] = inputValues[reordering[reorderIndex++][r]][c];
			}
		}
		return data;
	}
	
	/**
	 * Reshapes the given single dimensional array into a 2D array of the given
	 *  size
	 * 
	 * @param data
	 * @param rows
	 * @param columns
	 * @return
	 */
	public static double[][] reshape(double[] data, int rows, int columns) {
		double[][] matrix = new double[rows][columns];
		int i = 0;
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				matrix[r][c] = data[i++];
			}
		}
		return matrix;
	}
	
	/**
	 * Constructs all embedding vectors of size k for the data.
	 * Will be data.length - k + 1 of these
	 * 
	 * @param data
	 * @param k
	 * @return
	 */
	public static double[][] makeDelayEmbeddingVector(double[] data, int k) {
		try {
			return makeDelayEmbeddingVector(data, k, k - 1, data.length - k + 1);
		} catch (Exception e) {
			// The above call should not throw an Exception, handle here 
			//  in a RuntimeException so this method doesn't throw one
			throw new RuntimeException(e);
		}
	}

	/**
	 * Constructs numEmbeddingVectors embedding vectors of size k for the data,
	 * with the first embedding vector having it's last time point at t=startKthPoint
	 * 
	 * @param data
	 * @param k
	 * @param startKthPoint
	 * @param numEmbeddingVectors
	 * @return
	 */
	public static double[][] makeDelayEmbeddingVector(double[] data, int k,
			int startKthPoint, int numEmbeddingVectors) throws Exception {
		if (startKthPoint < k - 1) {
			throw new Exception("Start point t=" + startKthPoint + " is too early for a " +
					k + " length embedding vector");
		}
		if (numEmbeddingVectors + startKthPoint > data.length) {
			throw new Exception("Too many embedding vectors " + numEmbeddingVectors +
					" requested for the given startPoint " + startKthPoint +
					" and time series length " + data.length);
		}
		double[][] embeddingVectors = new double[numEmbeddingVectors][k];
		for (int t = startKthPoint; t < numEmbeddingVectors + startKthPoint; t++) {
			for (int i = 0; i < k; i++) {
				embeddingVectors[t - startKthPoint][i] = data[t - i];
			}
		}
		return embeddingVectors;
	}

	/**
	 * Constructs all embedding vectors of size k for the data.
	 * Will be data.length - k + 1 of these
	 * 
	 * @param data
	 * @param k
	 * @return
	 */
	public static double[][] makeDelayEmbeddingVector(double[][] data, int k) {
		try {
			return makeDelayEmbeddingVector(data, k, k - 1, data.length - k + 1);
		} catch (Exception e) {
			// The above call should not throw an Exception, handle here 
			//  in a RuntimeException so this method doesn't throw one
			throw new RuntimeException(e);
		}
	}

	/**
	 * Constructs numEmbeddingVectors embedding vectors of size k for the data,
	 * with the first embedding vector having it's last time point at t=startKthPoint
	 * 
	 * @param data
	 * @param k
	 * @param startKthPoint
	 * @param numEmbeddingVectors
	 * @return
	 */
	public static double[][] makeDelayEmbeddingVector(double[][] data, int k,
			int startKthPoint, int numEmbeddingVectors) throws Exception {
		if (startKthPoint < k - 1) {
			throw new Exception("Start point t=" + startKthPoint + " is too early for a " +
					k + " length embedding vector");
		}
		if (numEmbeddingVectors + startKthPoint > data.length) {
			throw new Exception("Too many embedding vectors " + numEmbeddingVectors +
					" requested for the given startPoint " + startKthPoint +
					" and time series length " + data.length);
		}
		int columns = data[0].length;
		double[][] embeddingVectors = new double[numEmbeddingVectors][k * columns];
		for (int t = startKthPoint; t < numEmbeddingVectors + startKthPoint; t++) {
			for (int i = 0; i < k; i++) {
				for (int c = 0; c < columns; c++) {
					embeddingVectors[t - startKthPoint][i*columns + c] = data[t - i][c];
				}
			}
		}
		return embeddingVectors;
	}

	public static int[] subArray(int[] array, int startIndex, int theLength) {
		int[] sub = new int[theLength];
		for (int r = 0; r < theLength; r++) {
			sub[r] = array[startIndex + r];
		}
		return sub;
	}
	
	public static double stdDev(double[] array) {
		double mean = 0.0;
		double total = 0.0;
		for (int m = 0; m < array.length; m++) {
			total += array[m];
		}
		mean = total / (double) array.length;
		
		return stdDev(array, mean);
	}
	
	public static double stdDev(double[][] matrix, int column) {
		double mean = 0.0;
		double total = 0.0;
		for (int m = 0; m < matrix.length; m++) {
			total += matrix[m][column];
		}
		mean = total / (double) matrix.length;
		
		return stdDev(matrix, column, mean);
	}
	
	/**
	 * Return the standard deviation of all the elements in array
	 * 
	 * @param array
	 * @param mean
	 * @return
	 */
	public static double stdDev(double[] array, double mean) {
		return stdDev(array, mean, array.length);
	}
	
	/**
	 * Standard deviation for the first arrayLength terms of array
	 * 
	 * @param array
	 * @param mean
	 * @param arrayLength
	 * @return
	 */
	public static double stdDev(double[] array, double mean, int arrayLength) {
		if (arrayLength == 0) {
			return 0.0;
		}
		double sumSqs = 0.0;
		for (int m = 0; m < arrayLength; m++) {
			sumSqs += (array[m] - mean) * (array[m] - mean);
		}
		double std = sumSqs / (double) (arrayLength - 1);
		std = Math.sqrt(std);
		return std;
	}
	
	/**
	 * Compute the standard deviation along the given column, with the known
	 *  given mean.
	 * 
	 * @param matrix
	 * @param column
	 * @param mean
	 * @return
	 */
	public static double stdDev(double[][] matrix, int column, double mean) {
		if (matrix.length == 0) {
			return 0.0;
		}
		double sumSqs = 0.0;
		for (int m = 0; m < matrix.length; m++) {
			sumSqs += (matrix[m][column] - mean) * (matrix[m][column] - mean);
		}
		double std = sumSqs / (double) (matrix.length - 1);
		std = Math.sqrt(std);
		return std;
	}
	
	/**
	 * Compute the standard deviation across all values in the 2D matrix
	 * 
	 * @param matrix
	 * @return
	 */
	public static double stdDev(double[][] matrix) {
		double mean = mean(matrix);
		return stdDev(matrix, mean);
	}
	
	/**
	 * Compute the standard deviation across all values in the 2D matrix
	 * 
	 * @param matrix
	 * @param mean
	 * @return
	 */
	public static double stdDev(double[][] matrix, double mean) {
		if (matrix.length == 0) {
			return 0.0;
		}
		double sumSqs = 0.0;
		for (int m = 0; m < matrix.length; m++) {
			for (int c = 0; c < matrix[m].length; c++) {
				sumSqs += (matrix[m][c] - mean) * (matrix[m][c] - mean);
			}
		}
		double std = sumSqs / (double) ((matrix.length * matrix[0].length) - 1);
		std = Math.sqrt(std);
		return std;
	}


	/**
	 * Compute the standard deviations along each column
	 * 
	 * @param matrix
	 * @param means
	 * @return
	 */
	public static double[] stdDevs(double[][] matrix, double[] means) {
		double[] sumSqs = new double[means.length];
		for (int m = 0; m < matrix.length; m++) {
			for (int c = 0; c < matrix[m].length; c++) {
				sumSqs[c] += (matrix[m][c] - means[c]) * (matrix[m][c] - means[c]);
			}
		}
		double[] stds = new double[means.length];
		for (int c = 0; c < stds.length; c++) {
			stds[c] = sumSqs[c] / (double) (matrix.length - 1);
			stds[c] = Math.sqrt(stds[c]);
		}
		return stds;
	}

	/**
	 * Compute the standard deviations along each row
	 * 
	 * @param matrix
	 * @param means
	 * @return
	 */
	public static double[] stdDevsOfRows(double[][] matrix, double[] means) {
		double[] stds = new double[matrix.length];
		for (int r = 0; r < matrix.length; r++) {
			double sumSqs = 0.0;
			for (int c = 0; c < matrix[r].length; c++) {
				sumSqs += (matrix[r][c] - means[r]) * (matrix[r][c] - means[r]);
			}
			stds[r] = sumSqs / (double) (matrix[r].length - 1);
			stds[r] = Math.sqrt(stds[r]);
		}
		return stds;
	}

	public static double max(double[][][] matrix) {
		// double max = 0.0;
		double max = matrix[0][0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				for (int k = 0; k < matrix[i][j].length; k++) {
					if (matrix[i][j][k] > max) {
						max = matrix[i][j][k];
					}
				}
			}
		}
		return max;
	}

	/**
	 * Normalises the elements in the given array
	 * 
	 * @param array
	 */
	public static void normalise(double[] array) {
		double mean = MatrixUtils.mean(array);
		double stdDev = MatrixUtils.stdDev(array, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < array.length; t++) {
				array[t] = (array[t] - mean);
			}			
		} else {
			// stdDev is non zero
			for (int t = 0; t < array.length; t++) {
				array[t] = (array[t] - mean) / stdDev;
			}
		}
	}
	
	/**
	 * Returns a normalised array of the elements in the given array 
	 * 
	 * @param array
	 */
	public static double[] normaliseIntoNewArray(double[] array) {
		double[] newArray = new double[array.length];
		double mean = MatrixUtils.mean(array);
		double stdDev = MatrixUtils.stdDev(array, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < array.length; t++) {
				newArray[t] = (array[t] - mean);
			}
		} else {
			// stdDev is non zero
			for (int t = 0; t < array.length; t++) {
				newArray[t] = (array[t] - mean) / stdDev;
			}
		}
		return newArray;
	}

	/**
	 * Normalises the elements in the given column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 * @param column column number to be normalised
	 */
	public static void normalise(double[][] matrix, int column) {
		double mean = MatrixUtils.mean(matrix, column);
		double stdDev = MatrixUtils.stdDev(matrix, column, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < matrix.length; t++) {
				matrix[t][column] = (matrix[t][column] - mean);
			}
		} else {
			// stdDev is non zero
			for (int t = 0; t < matrix.length; t++) {
				matrix[t][column] = (matrix[t][column] - mean) / stdDev;
			}
		}
	}

	/**
	 * Normalises the elements in the given column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 * @param column column number to be normalised
	 */
	public static double[] normaliseIntoNewArray(double[][] matrix, int column) {
		double[] newArray = new double[matrix.length];
		double mean = MatrixUtils.mean(matrix, column);
		double stdDev = MatrixUtils.stdDev(matrix, column, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < matrix.length; t++) {
				newArray[t] = (matrix[t][column] - mean);
			}
		} else {
			// stdDev is non zero
			for (int t = 0; t < matrix.length; t++) {
				newArray[t] = (matrix[t][column] - mean) / stdDev;
			}
		}
		return newArray;
	}

	/**
	 * Normalises the elements along each column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 */
	public static double[][] normaliseIntoNewArray(double[][] matrix) {
		double[][] newMatrix = new double[matrix.length][matrix[0].length];
		double[] means = means(matrix);
		double[] stds = stdDevs(matrix, means);
		for (int r = 0; r < newMatrix.length; r++) {
			for (int c = 0; c < newMatrix[r].length; c++) {
				newMatrix[r][c] = matrix[r][c] - means[c];
				if (!Double.isInfinite(1.0 /  stds[c])) {
					newMatrix[r][c] /= stds[c];
				} // else we just subtract off the mean
			}
		}
		return newMatrix;		
	}

	public static double max(double[][] matrix) {
		// double max = 0.0;
		double max = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (Double.isNaN(max) || (matrix[i][j] > max)) {
					max = matrix[i][j];
				}
			}
		}
		return max;
	}

	public static int max(int[][] matrix) {
		// int max = 0;
		int max = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (matrix[i][j] > max) {
					max = matrix[i][j];
				}
			}
		}
		return max;
	}

	public static double max(double[] array) {
		return maxStartFromIndex(array, 0);
	}

	public static double maxStartFromIndex(double[] array, int startFromIndex) {
		// double max = 0.0;
		double max = array[startFromIndex];
		for (int i = startFromIndex; i < array.length; i++) {
			if (Double.isNaN(max) || (array[i] > max)) {
				max = array[i];
			}
		}
		return max;
	}

	public static int max(int[] array) {
		// int max = 0;
		int max = array[0];
		for (int i = 0; i < array.length; i++) {
			if (array[i] > max) {
				max = array[i];
			}
		}
		return max;
	}

	/**
	 * Works out the maximum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static double max(double[][] matrix, int column) {
		// double max = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		double max = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (Double.isNaN(max) || (matrix[i][column] > max)) {
				max = matrix[i][column];
			}
		}
		return max;
	}

	/**
	 * Works out the maximum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static int max(int[][] matrix, int column) {
		// double max = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		int max = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (matrix[i][column] > max) {
				max = matrix[i][column];
			}
		}
		return max;
	}

	public static double min(double[][] matrix) {
		// double min = 0.0;
		double min = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (Double.isNaN(min) || (matrix[i][j] < min)) {
					min = matrix[i][j];
				}
			}
		}
		return min;
	}

	public static int min(int[][] matrix) {
		// int min = 0;
		int min = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (matrix[i][j] < min) {
					min = matrix[i][j];
				}
			}
		}
		return min;
	}

	public static double min(double[] array) {
		return minStartFromIndex(array, 0);
	}

	/**
	 * Find the kth minimum value in the array.
	 * 
	 * @param array
	 * @param k
	 * @return
	 * @throws Exception 
	 */
	public static double kthMin(double[] array, int k) throws Exception {
		if (k == 1) {
			return min(array); 
		}
		if (array.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					array.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
		}
		for (int t = 0; t < array.length; t++) {
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (array[t] < mins[k - 1]) {
				mins[k - 1] = array[t];
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (array[t] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						mins[i] = array[t];
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the kth min
		return mins[k-1];
	}

	/**
	 * Find the kth minimum value in the array subject to
	 *  a given condition.
	 * Assumes that the condition is satisfied at least k
	 *  times in the array (this is not checked in here)
	 * 
	 * @param array
	 * @param k
	 * @param extraData
	 * @param extraCondition
	 * @return
	 * @throws Exception 
	 */
	public static double kthMinSubjectTo(double[] array, int k, int[] extraData, int condition) throws Exception {
		// Can't do a quickie for k==1 here since we're subject to 
		//  checking the condition
		if (array.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					array.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
		}
		for (int t = 0; t < array.length; t++) {
			if (extraData[t] != condition) {
				continue;
			}
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (array[t] < mins[k - 1]) {
				mins[k - 1] = array[t];
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (array[t] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						mins[i] = array[t];
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the kth min
		return mins[k-1];
	}

	public static double minIgnoreIndex(double[] array, int indexToIgnore) {
		// double min = 0.0;
		double min;
		if (indexToIgnore != 0) {
			min = array[0];
		} else {
			min = array[1];
		}
		for (int i = 0; i < array.length; i++) {
			if (indexToIgnore == i) {
				continue;
			}
			if (Double.isNaN(min) || (array[i] < min)) {
				min = array[i];
			}
		}
		return min;
	}

	public static double minStartFromIndex(double[] array, int startFromIndex) {
		// double min = 0.0;
		double min = array[startFromIndex];
		for (int i = startFromIndex; i < array.length; i++) {
			if (Double.isNaN(min) || (array[i] < min)) {
				min = array[i];
			}
		}
		return min;
	}

	public static int min(int[] array) {
		// int min = 0;
		int min = array[0];
		for (int i = 0; i < array.length; i++) {
			if (array[i] < min) {
				min = array[i];
			}
		}
		return min;
	}
	
	/**
	 * Works out the minimum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static double min(double[][] matrix, int column) {
		// double min = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		double min = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (Double.isNaN(min) || (matrix[i][column] < min)) {
				min = matrix[i][column];
			}
		}
		return min;
	}

	/**
	 * Works out the minimum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static int min(int[][] matrix, int column) {
		// double min = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		int min = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (matrix[i][column] < min) {
				min = matrix[i][column];
			}
		}
		return min;
	}

	/**
	 * Works out the index of the minimum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static int minIndex(double[][] matrix, int column) {
		// double min = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		double min = matrix[0][column];
		int minIndex = 0;
		for (int i = 1; i < matrix.length; i++) {
			if (Double.isNaN(min) || (matrix[i][column] < min)) {
				min = matrix[i][column];
				minIndex = i;
			}
		}
		return minIndex;
	}

	/**
	 * Works out the index of the k minimum values in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @param k
	 * @return
	 * @throws Exception 
	 */
	public static int[] kMinIndices(double[][] matrix, int column, int k) throws Exception {
		if (matrix.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					matrix.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		int[] minIndices = new int[k];
		if (k == 1) {
			minIndices[0] = minIndex(matrix, column);
			return minIndices;
		}
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
			minIndices[i] = -1;
		}
		for (int t = 0; t < matrix.length; t++) {
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (matrix[t][column] < mins[k - 1]) {
				mins[k - 1] = matrix[t][column];
				minIndices[k-1] = t;
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (matrix[t][column] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						minIndices[i+1] = minIndices[i];
						mins[i] = matrix[t][column];
						minIndices[i] = t;
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the index of the kth min
		return minIndices;
	}

	/**
	 * Works out the index of the k minimum values in the matrix in a given column 
	 *  subject to the extraData matching a given condition.
	 * We do not check whether there are k matches for the extraData to
	 *  the condition here - the caller should check this themselves.
	 * 
	 * @param matrix
	 * @param column
	 * @param k
	 * @param extraData
	 * @param condition
	 * @return
	 * @throws Exception 
	 */
	public static int[] kMinIndicesSubjectTo(double[][] matrix, int column,
			int k, int[] extraData, int condition) throws Exception {
		if (matrix.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					matrix.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		int[] minIndices = new int[k];
		// no quick check for k==1 since we need to check the extra condition
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
			minIndices[i] = -1;
		}
		for (int t = 0; t < matrix.length; t++) {
			if (extraData[t] != condition) {
				continue;
			}
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (matrix[t][column] < mins[k - 1]) {
				mins[k - 1] = matrix[t][column];
				minIndices[k-1] = t;
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (matrix[t][column] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						minIndices[i+1] = minIndices[i];
						mins[i] = matrix[t][column];
						minIndices[i] = t;
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the index of the kth min
		return minIndices;
	}

	/**
	 * Mirrors the matrix in both coordinates
	 * 
	 * @param matrix
	 * @return
	 */
	public static int[][] mirrorMatrixBothCoords(int[][] matrix) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		
		int[][] mirrored = new int[rows][cols];
		
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				mirrored[(rows - 1) - r][(cols - 1) - c] = matrix[r][c];
			}
		}
		return mirrored;
	}
	
	/**
	 * Mirrors the matrix in both coordinates
	 * 
	 * @param matrix
	 * @return
	 */
	public static double[][] mirrorMatrixBothCoords(double[][] matrix) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		
		double[][] mirrored = new double[rows][cols];
		
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				mirrored[(rows - 1) - r][(cols - 1) - c] = matrix[r][c];
			}
		}
		return mirrored;
	}

	/**
	 * Moves the rows of the array up by upBy.
	 * Inserts zeros at the bottom
	 * 
	 * @param matrix
	 * @param upBy
	 */
	public static void moveRowsUp(double[][] matrix, int upBy) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		
		for (int r = 0; r < rows - upBy; r++) {
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = matrix[r + upBy][c];
			}
		}
		for (int r = rows - upBy; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = 0;
			}
		}
	}
	
	/**
	 * <p>Returns the covariance between the two columns of data.</p>
	 * <p>See - <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 * </p>
	 * 
	 * @param data
	 * @return the covariance
	 */
	public static double covariance(double[][] data) {
		double c = 0;
		double mean0 = mean(data, 0);
		double mean1 = mean(data, 1);
		for (int t = 0; t < data.length; t++) {
			c += (data[t][0] - mean0)*(data[t][1]-mean1);
		}
		return c / (double) data.length;
	}
	
	/**
	 * <p>Returns the covariance between the two arrays of data.</p>
	 * <p>See - <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @return the covariance
	 */
	public static double covariance(double[] x, double[] y) {
		double c = 0;
		double meanX = mean(x);
		double meanY = mean(y);
		for (int t = 0; t < x.length; t++) {
			c += (x[t] - meanX)*(y[t]-meanY);
		}
		return c / (double) x.length;
	}

	/**
	 * <p>Returns the correlation between the two arrays of data.</p>
	 * <p>The arrays are asssumed to have the same lengths</p>
	 * <p>See - <a href="http://en.wikipedia.org/wiki/Correlation">Wikipedia</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @return the covariance
	 */
	public static double correlation(double[] x, double[] y) {
		return correlation(x, y, x.length);
	}
	
	/**
	 * <p>Returns the correlation between the two arrays of data.</p>
	 * <p>See - <a href="http://en.wikipedia.org/wiki/Correlation">Wikipedia</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @param dataLength - number of terms in each vector to consider (we look at the first dataLength terms).
	 * 	Precondition: dataLength is less than min(x.length, y.length)
	 * @return the covariance
	 */
	public static double correlation(double[] x, double[] y, int dataLength) {
		// return covariance(x, y) / stdDev(x) / stdDev(y);
		// Save some code time by reusing the code from inside covariance:
		double c = 0;
		double meanX = mean(x, 0, dataLength);
		double meanY = mean(y, 0, dataLength);
		for (int t = 0; t < dataLength; t++) {
			c += (x[t] - meanX)*(y[t]-meanY);
		}
		double covariance = c / (double) (dataLength - 1);
		return covariance / stdDev(x, meanX, dataLength) / stdDev(y, meanY, dataLength);
	}
	
	/**
	 * <p>Returns the correlation between the two arrays of data,
	 *    ignoring any Nan values</p>
	 * <p>See - <a href="http://en.wikipedia.org/wiki/Correlation">Wikipedia</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @param dataLength - number of terms in each vector to consider (we look at the first dataLength terms).
	 * 	Precondition: dataLength is less than min(x.length, y.length)
	 * @return the covariance
	 */
	public static double correlationIgnoreNans(double[] x, double[] y, int dataLength) {
		// return covariance(x, y) / stdDev(x) / stdDev(y);
		// Save some code time by reusing the code from inside covariance:
		double c = 0;
		double meanX = 0;
		double meanY = 0;
		int count = 0;
		for (int i = 0; i < dataLength; i++) {
			if ((!Double.isNaN(x[i])) && (!Double.isNaN(y[i]))) {
				// Only add the values in if they are not NaN
				meanX += x[i];
				meanY += y[i];
				count++;
			}
		}
		// Adjust for the values we've skipped:
		meanX = meanX / count;
		meanY = meanY / count;
		
		for (int t = 0; t < dataLength; t++) {
			if ((!Double.isNaN(x[t])) && (!Double.isNaN(y[t]))) {
				// Only add the product in if it is not NaN
				c += (x[t] - meanX)*(y[t]-meanY);
			}
		}
		double covariance = c / (double) (count - 1);
		
		// Now work out the std devs of each:
		double sumSqsX = 0.0;
		double sumSqsY = 0.0;
		for (int m = 0; m < dataLength; m++) {
			if ((!Double.isNaN(x[m])) && (!Double.isNaN(y[m]))) {
				// Ignore if one is NaN
				sumSqsX += (x[m] - meanX) * (x[m] - meanX);
				sumSqsY += (y[m] - meanY) * (y[m] - meanY);
			}
		}
		double stdX = sumSqsX / (double) (count - 1);
		stdX = Math.sqrt(stdX);
		double stdY = sumSqsY / (double) (count - 1);
		stdY = Math.sqrt(stdY);
		
		return covariance / stdX / stdY;
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			matrix[r] = value;
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 * @param offset where in the array to start from
	 * @param length length in the array to fill
	 */
	public static void fill(int[] matrix, int value, int offset, int length) {
		for (int r = offset; r < offset + length; r++) {
			matrix[r] = value;
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[][] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = value;
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[][][] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					matrix[r][c][h] = value;
				}
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[][][][] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					int depth = matrix[r][c][h].length;
					for (int d = 0; d < depth; d++) {
						matrix[r][c][h][d] = value;
					}
				}
			}
		}
	}
	
	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			matrix[r] = value;
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[][] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = value;
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[][][] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					matrix[r][c][h] = value;
				}
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[][][][] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					int depth = matrix[r][c][h].length;
					for (int d = 0; d < depth; d++) {
						matrix[r][c][h][d] = value;
					}
				}
			}
		}
	}
	
	public static double[][] transpose(double[][] matrix) {
		double[][] newMatrix = new double[matrix[0].length][matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				newMatrix[j][i] = matrix[i][j];
			}
		}
		return newMatrix;
	}
	
	public static int[][] transpose(int[][] matrix) {
		int[][] newMatrix = new int[matrix[0].length][matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				newMatrix[j][i] = matrix[i][j];
			}
		}
		return newMatrix;
	}

	/**
	 * Converts an int array to a double array
	 * 
	 * @param input
	 * @return
	 */
	public static double[][] convertMatrix(int[][] input) {
		double[][] outputArray = new double[input.length][];
		for (int i = 0; i < input.length; i++) {
			outputArray[i] = new double[input[i].length];
			for (int j = 0; j < input[i].length; j++) {
				outputArray[i][j] = input[i][j];
			}
		}
		return outputArray;
	}
	
	/**
	 * Converts a double array to an int array
	 * 
	 * @param input
	 * @param valueOffset value to be subtracted from each value
	 * @return
	 */
	public static int[] convertMatrix(double[] input, int valueOffset) {
		int[] outputArray = new int[input.length];
		for (int i = 0; i < input.length; i++) {
			outputArray[i] = (int) (input[i]) - valueOffset;
		}
		return outputArray;
	}

	/**
	 * Converts a double array to an int array
	 * 
	 * @param input
	 * @return
	 */
	public static int[] convertMatrix(double[] input) {
		return convertMatrix(input, 0);
	}

	/**
	 * <p> Returns the determinant of the input matrix.
	 * </p>
	 * 
	 * @param matrix
	 * @return determinant of matrix
	 */
	public static double determinant(double[][] matrix) throws Exception {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			if (matrix[r].length != rows) {
				throw new Exception("Cannot compute the determinant of a non-square matrix");
			}
		}
		return recursiveDeterminant(matrix);
	}
	
	/**
	 * <p>Private function to compute the determinant recursively.
	 * determinant() calls this after checking the matrix dimensions. <br/>
	 * See - http://mathworld.wolfram.com/Determinant.html
	 * </p>
	 * 
	 * @param matrix
	 * @return
	 */
	private static double recursiveDeterminant(double[][] matrix) {
		int rows = matrix.length;
		double result = 0;
		// Base cases:
		if (rows == 1) {
			return matrix[0][0];
		}
		if (rows == 2) {
			return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
		}
		// Recursive case
		int multiplier = 1;
		for(int col = 0; col < rows; col++) {
			// Construct the next sub-matrix to compute the determinant of
			double minor[][] = copyMatrixEliminateRowAndColumn(matrix, 0, col);
			result += (double) multiplier * matrix[0][col] * recursiveDeterminant(minor);
			multiplier *= -1;
		}

		return result; 
	}
	
	public static void printMatrix(PrintStream out, double[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				out.print(matrix[r][c] + " ");
			}
			out.println();
		}
	}
	
	public static void printMatrix(PrintStream out, int[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				out.print(matrix[r][c] + " ");
			}
			out.println();
		}
	}

	public static void printArray(PrintStream out, double[] array) {
		for (int r = 0; r < array.length; r++) {
				out.print(array[r] + " ");
		}
		out.println();
	}

	public static void printArray(PrintStream out, int[] array) {
		for (int r = 0; r < array.length; r++) {
				out.print(array[r] + " ");
		}
		out.println();
	}
	
	/**
	 * Discretizes using even bin sizes
	 * 
	 * @param data
	 * @param base
	 * @return
	 */
	public static int[] discretise(double data[], int base) {
		int[] discretised = new int[data.length];
		double min = min(data);
		double max = max(data);
		double binInterval = (max - min) / base;
		
		for (int t = 0; t < data.length; t++) {
			discretised[t] = (int) ((data[t] - min) / binInterval);
			if (discretised[t] == base) {
				// This occurs for the maximum value; put it in the largest bin (base - 1)
				discretised[t]--;
			}
		}
		return discretised;
	}
	
	/**
	 * Discretizes each column of the data independently,
	 * using a maximum entropy partitioning
	 * 
	 * @param data
	 * @param base
	 * @return
	 */
	public static int[][] discretiseMaxEntropy(double data[][], int base){
		int lastCol = data[0].length;	
		int lastRow = data.length;
		int[][] newData = new int[lastRow][lastCol];
		for(int j=0;j<lastCol;j++){
			double[] tempData = new double[lastRow];
			for (int i=0;i<lastRow;i++){
				tempData[i] = data[i][j]; 			 
			}
		
			Arrays.sort(tempData);

			int compartmentSize;
			double[] cutOffValues = new double[base]; 
			for(int i=0;i<base;i++){
				compartmentSize = (int)((double)(i+1)*(double)(lastRow)/(double)base)-1;
//				System.out.println(compartmentSize);
				cutOffValues[i]=tempData[compartmentSize];
			}
			
			for (int i=0;i<lastRow;i++){				
				for(int m=0;m<base;m++){
					if (data[i][j] <= cutOffValues[m]){
						newData[i][j] = m;
						m = base;
					}
				}
			}
		}
		return newData;
	}
	
	/**
	 * Take the logical AND of all variables in each row 
	 * 
	 * @param data
	 * @return
	 */
	public static boolean[] andRows(boolean[][] data) {
		boolean[] result = new boolean[data.length];
		
		for (int i = 0; i < data.length; i++) {
			result[i] = true;
			for (int j = 0; j < data[i].length; j++) {
				result[i] &= data[i][j];
			}
		}
		return result;
	}
	
	/**
	 * Take the logical AND of selected variables in each row 
	 * 
	 * @param data
	 * @param columns which variables to take the AND over
	 * @return
	 */
	public static boolean[] andRowsOverSelectedColumns(boolean[][] data, int[] columns) {
		boolean[] result = new boolean[data.length];
		
		for (int i = 0; i < data.length; i++) {
			result[i] = true;
			for (int c = 0; c < columns.length; c++) {
				result[i] &= data[i][columns[c]];
			}
		}
		return result;
	}
	
	/**
	 * Convert a double array to an int array.
	 * This is designed specifically for use of the toolkit in Octave
	 *  where all native arrays are considered as doubles for Java,
	 *  and octave-java can't properly identify valid method signatures
	 *  unless the arrays are converted to int arrays first.
	 * 
	 * @param array
	 * @return
	 */
	public static int[] doubleToIntArray(double[] array) {
		if (array == null) {
			return null;
		}
		int[] intArray = new int[array.length];
		for (int i = 0; i < array.length; i++) {
			intArray[i] = (int) array[i];
		}
		return intArray;
	}

	/**
	 * Convert a 2D double array to an int array.
	 * This is designed specifically for use of the toolkit in Octave
	 *  where all native arrays are considered as doubles for Java,
	 *  and octave-java can't properly identify valid method signatures
	 *  unless the arrays are converted to int arrays first.
	 * 
	 * @param array
	 * @return
	 */
	public static int[][] doubleToIntArray(double[][] array) {
		if (array == null) {
			return null;
		}
		int[][] intArray = new int[array.length][];
		for (int i = 0; i < array.length; i++) {
			if (array[i] == null) {
				intArray[i] = null;
			} else {
				intArray[i] = new int[array[i].length];
				for (int j = 0; j < array[i].length; j++) {
					intArray[i][j] = (int) array[i][j];
				}
			}
		}
		return intArray;
	}
}
