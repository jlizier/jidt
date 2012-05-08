package infodynamics.utils;

import java.util.Arrays;

/**
 * 
 * @author Joseph Lizier
 *
 * LongArrayWrapper is used to wrap a long array so as to be able to give it a hash value
 *  for placement into a hash table.
 */
public class IntArrayWrapper {

	private int[] array;
	private int firstCols = 0;
	
	public IntArrayWrapper() {
		super();
	}

	public IntArrayWrapper(int[] newArray) {
		super();
		setArray(newArray, newArray.length);
	}

	public IntArrayWrapper(int[] newArray, int firstcolumns) {
		super();
		setArray(newArray, firstcolumns);
	}

	public int hashCode() {
		return Arrays.hashCode(array);
	}
	
	public int[] getArray() {
		return array;
	}

	public void setArray(int[] newArray, int firstcolumns) {
		array = newArray;
		firstCols = firstcolumns;
	}
	
	public int getArrayUsedLength() {
		return firstCols;
	}
	
	public boolean equals(Object o) {
		if (o == this) {
			return true;
		}
		
		IntArrayWrapper iaw2 = (IntArrayWrapper) o;
		if (iaw2.getArrayUsedLength() != firstCols) {
			return false;
		}
		
		// check deeply inside arrays
		int[] inArray = iaw2.getArray();
		for (int i = 0; i < firstCols; i++) {
			if (array[i] != inArray[i]) {
				return false;
			}
		}
		return true;
	}
}
