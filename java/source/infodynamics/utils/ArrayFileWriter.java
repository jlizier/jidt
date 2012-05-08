package infodynamics.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ArrayFileWriter {
	
	/**
	 * Outputs different types of 2D Matrix to a file
	 * 
	 * @param outputFileName and matrix
	 */
	
	public static void makeBooleanMatrixFile(boolean matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;
		
		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				if (matrix[i][j]==false){
					out.write("0\t");
				}else{
					out.write("1\t");
				}
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}
	
	public static void makeDoubleMatrixFile(double matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.valueOf(matrix[i][j])+"\t");
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	public static void makeDoubleMatrixFile(double matrix[][],String outputFileName, int decimalPlaces) throws IOException {
		String template = String.format("%%.%df\t", decimalPlaces);
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.format(template, matrix[i][j]));
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	public static void makeIntMatrixFile(int matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.valueOf(matrix[i][j])+"\t");
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	public static void makeMatrixFile(double matrix[],String outputFileName) throws IOException {
    	int rowSize = matrix.length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			out.write(String.valueOf(matrix[i]) + "\n");
		}
		out.close();
	}

	public static void makeMatrixFile(Number matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.valueOf(matrix[i][j])+"\t");
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	public static void makeMatrixFile(Number matrix[],String outputFileName) throws IOException {
    	int rowSize = matrix.length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			out.write(String.valueOf(matrix[i]) + "\n");
		}
		out.close();
	}
	
	private static void createDirectories(String filename) {
		File file = new File(filename);
		File parentDir = file.getParentFile();
		if ((parentDir != null) && !parentDir.isDirectory()) {
			parentDir.mkdirs();
		}
	}
}
