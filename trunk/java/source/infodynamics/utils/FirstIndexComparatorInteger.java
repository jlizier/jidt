package infodynamics.utils;

/**
 * <p>Class used to sort an array of array of doubles based on the first element in each array.
 * </p> 
 * 
 * <p>Usage: <code>java.util.Arrays.sort(int[][], FirstIndexComparatorInteger.getInstance())</code>;
 * </p>
 * 
 * @author Joseph Lizier
 */
public class FirstIndexComparatorInteger implements java.util.Comparator<int[]> {
	
	private static FirstIndexComparatorInteger instance = null;
	
	private FirstIndexComparatorInteger() {
	}

	public static FirstIndexComparatorInteger getInstance() {
		if (instance == null) {
			instance = new FirstIndexComparatorInteger();
		}
		return instance;
	}
	
	public int compare(int[] a1, int[] a2) {
		if (a1[0] < a2[0]) {
			return -1;
		} else if (a1[0] == a2[0]) {
			return 0;
		} else {
			return 1;
		}
	}
}