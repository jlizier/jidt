package infodynamics.utils;

public class MathsUtils {

	private static final double EULER_MASCHERONI_CONSTANT = 0.5772156;
	
	private static int highestDigammaArgCalced = 0;
	private static final int NUM_STORED_DIGAMMAS = 10000;
	private static double[] storedDigammas;
	
	public MathsUtils() {
		super();
		// TODO Auto-generated constructor stub
	}

	/**
	 * Returns the integer result of base^power
	 * 
	 * @param base base integer of the operation
	 * @param power power that base is raised to
	 * @return base raised to exponent power (rounded by integer operations)
	 */
	public static int power(int base, int power) {
		int result = 1;
		int absPower = Math.abs(power);
		for (int p = 0; p < absPower; p++) {
			result *= base;
		}
		if (power < 0) {
			// This will be zero for any base except 1 or -1
			result = 1 / result;
		}
		return result;
	}
	
	/**
	 * Returns the integer result of base^power
	 * 
	 * Tested - works.
	 * 
	 * @param base
	 * @param power
	 * @return
	 */
	public static long power(long base, long power) {
		long result = 1;
		long absPower = Math.abs(power);
		for (long p = 0; p < absPower; p++) {
			result *= base;
		}
		if (power < 0) {
			// This will be zero for any base except 1 or -1
			result = 1 / result;
		}
		return result;
	}
	
	public static long factorial(int n) {
		long result = 1;
		for (int i = 1; i <= n; i++) {
			result *= (long) i;
		}
		return result;
	}
	
	public static int factorialCheckBounds(int n) throws Exception {
		long result = 1;
		for (int i = 1; i <= n; i++) {
			result *= (long) i;
			if (result > Integer.MAX_VALUE) {
				throw new Exception("n! causes integer overflow");
			}
		}
		return (int) result;
	}

	public static double factorialAsDouble(int n) {
		double result = 1;
		for (int i = 1; i <= n; i++) {
			result *= (double) i;
		}
		return result;		
	}
	
	/**
	 * Computes n! as a double (to provide extended range over a long).
	 * 
	 * We include the divisor here since n! hits with n at about 340.
	 * So if there is anything the result would have been divided by, we include it here,
	 *  thus extending the range of n the function is suitable for. 
	 * 
	 * @param n
	 * @param divisor
	 * @return
	 */
	public static double factorialAsDoubleIncludeDivisor(int n, double divisor) {
		double result = 1.0 / divisor;
		for (int i = 1; i <= n; i++) {
			result *= (double) i;
		}
		return result;		
	}

	/**
	 * n!!
	 * see http://en.wikipedia.org/wiki/Double_factorial#Double_factorial
	 * 
	 * @param n
	 * @return
	 */
	public static long doubleFactorial(int n) {
		long result = 1;
		int startValue;
		if (n % 2 == 0) {
			// n even
			startValue = 2;
		} else {
			// n odd
			startValue = 3;
		}
		for (int i = startValue; i <= n; i += 2) {
			result *= (long) i;
		}
		return result;
	}

	/**
	 * n!!
	 * see http://en.wikipedia.org/wiki/Double_factorial#Double_factorial
	 * 
	 * @param n
	 * @return
	 */
	public static double doubleFactorialAsDouble(int n) {
		double result = 1.0;
		int startValue;
		if (n % 2 == 0) {
			// n even
			startValue = 2;
		} else {
			// n odd
			startValue = 3;
		}
		for (int i = startValue; i <= n; i += 2) {
			result *= (double) i;
		}
		return result;
	}
	
	/**
	 * n!!
	 * see http://en.wikipedia.org/wiki/Double_factorial#Double_factorial
	 * 
	 * We include the divisor here since the gamma(d/2+1) hits Inf with d just over 300.
	 * So if there is anything the result would have been divided by, we include it here,
	 *  thus extending the range of d the function is suitable for. 
	 * 
	 * @param n
	 * @param divisor
	 * @return
	 */
	public static double doubleFactorialAsDoublewithDivisor(int n, double divisor) {
		double result = 1.0 / divisor;
		int startValue;
		if (n % 2 == 0) {
			// n even
			startValue = 2;
		} else {
			// n odd
			startValue = 3;
		}
		for (int i = startValue; i <= n; i += 2) {
			result *= (double) i;
		}
		return result;
	}

	/**
	 * Computes gamma(d/2 + 1)
	 * See http://en.wikipedia.org/wiki/Gamma_function
	 *  for description of the analytical result for d odd.
	 * For d even, we have gamma of an integer, which is equal to
	 *  (d/2)!
	 * 
	 * @param d
	 * @return
	 */
	public static double gammaOfArgOn2Plus1(int d) {
		if (d % 2 == 0) {
			// d even
			return factorialAsDouble(d/2);
		} else {
			// d odd
			return Math.sqrt(Math.PI) * (double) doubleFactorialAsDouble(d) / 
				(double) Math.pow(2, ((double) (d + 1)) / 2.0);
		}
	}

	/**
	 * Computes gamma(d/2 + 1)/divisor
	 * See http://en.wikipedia.org/wiki/Gamma_function
	 *  for description of the analytical result for d odd.
	 * For d even, we have gamma of an integer, which is equal to
	 *  (d/2)!
	 * 
	 * We include the divisor here since the gamma(d/2+1) hits Inf with d just over 300.
	 * So if there is anything the result would have been divided by, we include it here,
	 *  thus extending the range of d the function is suitable for. 
	 * 
	 * @param d
	 * @param divisor
	 * @return
	 */
	public static double gammaOfArgOn2Plus1IncludeDivisor(int d, double divisor) {
		if (d % 2 == 0) {
			// d even
			return factorialAsDoubleIncludeDivisor(d/2, divisor);
		} else {
			// d odd
			return doubleFactorialAsDoublewithDivisor(d,
					divisor * Math.pow(2, ((double) (d + 1)) / 2.0) / Math.sqrt(Math.PI));
		}
	}

	/**
	 * Compute digamma(d).
	 * 
	 * Stores previous calculations to speed up computation here, though some precision may
	 *  be lost because we're adding in larger numbers first.
	 * 
	 * @param d
	 * @return
	 * @throws Exception
	 */
	public static double digamma(int d) throws Exception {
		if (d < 1) {
			return Double.NaN;
		}
		if (storedDigammas == null) {
			// allocate space to store our results
			storedDigammas = new double[NUM_STORED_DIGAMMAS];
			storedDigammas[0] = Double.NaN;
			storedDigammas[1] = -EULER_MASCHERONI_CONSTANT;
			highestDigammaArgCalced = 1;
		}
		if (d <= highestDigammaArgCalced) {
			// We've already calculated this one
			return storedDigammas[d];
		}
		// else need to calculate it
		double result = storedDigammas[highestDigammaArgCalced];
		for (int n = highestDigammaArgCalced + 1; n <= d; n++) {
			result += 1.0 / (double) (n-1);
			if (d < NUM_STORED_DIGAMMAS) {
				storedDigammas[n] = result;
			}
		}
		if (d < NUM_STORED_DIGAMMAS) {
			highestDigammaArgCalced = d;
		} else {
			highestDigammaArgCalced = NUM_STORED_DIGAMMAS - 1;
		}
		return result;		
	}
	
	/**
	 * Compute the digamma function from first principles
	 * 
	 * @param d
	 * @return
	 * @throws Exception
	 */
	public static double digammaByDefinition(int d) throws Exception {
		if (d < 1) {
			return Double.NaN;
		}
		double result = 0;
		for (int n = d; n > 1; n--) {
			result += 1.0 / (double) (n-1);
		}
		// Now add in result for n == 1
		result += -EULER_MASCHERONI_CONSTANT;
		return result;
	}

	/**
	 * Return the number of possible combinations of p from n (i.e. n choose p)
	 * 
	 * @param n
	 * @param p
	 * @return
	 * @throws Exception if the number would be greater than Integer.MAX_INT
	 */
	public static int numOfSets(int n, int p) throws Exception {
		// Compute how many sets there will be
		long counter = n;
		long numSets = 1;
		for (int x = 1; x <= p; x++) {
			numSets *= counter;
			numSets /= x;
			if (numSets > Integer.MAX_VALUE) {
				throw new Exception("nCp causes integer overflow");
			}
			counter--;
		}
		// numSets counts the number of permutations of n.
		// Need to get rid of repeats to make is combinations:
		return (int) numSets;
	}

	/**
	 * Return an array of all possible combinations of p from n 
	 * 
	 * @param n
	 * @param p
	 * @return
	 * @throws Exception when the number of sets is greaterr than Integer.MAX_INT
	 */
	public static int[][] generateAllSets(int n, int p) throws Exception {
		int numOfSets = numOfSets(n,p);
		int[][] sets = new int[numOfSets][p];
		int[] currentSet = new int[p];
		writeSetsIn(n, p, 0, 0, currentSet, sets, 0);
		return sets;
	}

	/**
	 * Recursive call used by generateAllSets.
	 * 
	 * @param n
	 * @param p
	 * @param currentIndexInSet current index in currentSet that we are writing into
	 * @param currentSet current set containing indices already written into the upper parts
	 * @param sets array to write generated sets into
	 * @param upToSetNum
	 * @return new value of upToSetNum
	 */
	private static int writeSetsIn(int n, int p, int currentIndexInSet,
			int firstCandidate, int[] currentSet, int[][] sets, int upToSetNum) {
		/*
		String indent = "";
		for (int i = 0; i < currentIndexInSet; i++) {
			indent += " ";
		}
		System.out.println(indent + String.format("currentIndex=%d", currentIndexInSet));
		*/
		// Put every candidate into this position:
		for (int candidate = firstCandidate; candidate < n - p + currentIndexInSet + 1; candidate++) {
			// System.out.println(indent + candidate);
			currentSet[currentIndexInSet] = candidate;
			if (currentIndexInSet == p - 1) {
				// We just wrote the last index, so copy this one in and return
				// System.out.println(indent + "writing into line " + upToSetNum);
				System.arraycopy(currentSet, 0, sets[upToSetNum++], 0, p);
			} else {
				// There are more indices to be written in, so make a recursive call to write the 
				//  next ones in
				upToSetNum = writeSetsIn(n, p, currentIndexInSet + 1, candidate + 1, currentSet, sets, upToSetNum);
			}
		}
		return upToSetNum;
	}

	public static void main(String args[]) throws Exception {
		/*
		System.out.println(numOfSets(158,4));
		System.out.println(numOfSets(158,3));
		System.out.println(numOfSets(158,2));
		*/
		// int[][] sets = generateAllSets(6,4);
		// MatrixUtils.printMatrix(System.out, sets);
		
		System.out.printf("digamma()  digammaOld()\n");
		for (int n = 0; n < 100; n++) {
			System.out.printf("%d  %.3f  %.3f\n", n, MathsUtils.digamma(n), MathsUtils.digammaByDefinition(n));
		}
		for (int n = 0; n < 101; n++) {
			System.out.printf("%d  %.3f  %.3f\n", n, MathsUtils.digamma(n), MathsUtils.digammaByDefinition(n));
		}
	}
}
