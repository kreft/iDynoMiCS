

/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

/**
 * @since Feb 2004
 * @version 1.0
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 */

package exceptions;

public class ModelRuntimeException extends RuntimeException {

	// recommended for classes that implement serializable
	private static final long serialVersionUID = 666L;

	/**
	 *  
	 */
	public ModelRuntimeException() {
		super();
	}
	/**
	 * @param message
	 */
	public ModelRuntimeException(String message) {
		super(message);
	}
	/**
	 * Create an exception to be thrown when a given number is not valid in
	 * runtime
	 * 
	 * @param f
	 */
	public ModelRuntimeException(double f) {
		super(f + " is not a valid number");
	}
}
