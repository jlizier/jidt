package infodynamics.utils;

/**
 * <p>Class used to sort an array of array of doubles based on the first element in each array.
 * </p> 
 * 
 * <p>Usage: <code>java.util.Arrays.sort(double[][], FirstIndexComparatorDouble.getInstance())</code>;
 * </p>
 * 
 * @author Joseph Lizier
 */
public class FirstIndexComparatorDouble implements java.util.Comparator<double[]> {
	
	private static FirstIndexComparatorDouble instance = null;
	
	private FirstIndexComparatorDouble() {
	}

	public static FirstIndexComparatorDouble getInstance() {
		if (instance == null) {
			instance = new FirstIndexComparatorDouble();
		}
		return instance;
	}
	
	public int compare(double[] a1, double[] a2) {
		if (a1[0] < a2[0]) {
			return -1;
		} else if (a1[0] == a2[0]) {
			return 0;
		} else {
			return 1;
		}
	}
}