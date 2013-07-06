package infodynamics.utils;

/**
 * An exception class to signal that an input matrix was not positive definite
 * 
 * @author Joseph Lizier
 *
 */
public class NonPositiveDefiniteMatrixException extends Exception {

	/**
	 * Default serialVersionUID
	 */
	private static final long serialVersionUID = 1L;

	public NonPositiveDefiniteMatrixException(String message) {
		super(message);
	}
}
