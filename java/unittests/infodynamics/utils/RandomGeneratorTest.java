package infodynamics.utils;

import junit.framework.TestCase;

public class RandomGeneratorTest extends TestCase {

	public void testGenerateRandomPerturbations() {
		checkGenerateRandomPerturbations(100, 20);
		// Now check if we ask for more perturbations than exist
		// (this should be ok)
		checkGenerateRandomPerturbations(5, 200);
		// Now test perturbations of a long set of digits
		//  to check timing (should be ~ 1 sec)
		checkGenerateRandomPerturbations(60000, 1000);
	}
	
	public void checkGenerateRandomPerturbations(int N, int p) {
		// Generate some perturbations:
		RandomGenerator rg = new RandomGenerator();
		int[][] perturbations = rg.generateRandomPerturbations(N, p);
		if (perturbations.length != p) {
			fail("Number of returned perturbations " + perturbations.length
					+ " does not match " + p);
		}
		// Make sure each one uses all of the available integers:
		for (int ip = 0; ip < p; ip++) {
			// Check that all integers < N used here:
			if (perturbations[ip].length != N) {
				fail("Length of perturbation " + ip + " is " +
						perturbations[ip].length + " instead of " + N);
			}
			boolean[] used = new boolean[N];
			for (int i = 0; i < N; i++) {
				if (used[perturbations[ip][i]]) {
					fail("Duplicate value " + perturbations[ip][i] +
							" in " + ip + "th array");
				}
				used[perturbations[ip][i]] = true;
			}
			// Postcondition: perturbation ip has no duplicates.
		}
	}

}
