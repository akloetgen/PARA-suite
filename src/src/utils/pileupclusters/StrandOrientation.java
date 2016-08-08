package utils.pileupclusters;

/**
 * Container for strand orientation.
 * 
 * @author akloetgen
 * 
 */
public class StrandOrientation {

	/**
	 * shows whether a read is set on the reverse strand of a reference
	 */
	private Boolean isReverse;

	/**
	 * 
	 * @return true if read is mapping on the reverse strand of a reference
	 *         sequence
	 */
	public boolean isReverse() {
		return isReverse;
	}

	/**
	 * 
	 * @param isReverse
	 *            sets local variable
	 */
	public void setReverse(Boolean isReverse) {
		this.isReverse = isReverse;
	}

	/**
	 * 
	 * @return true if read maps to reverse strand of a reference sequence
	 */
	public Boolean getReverse() {
		return isReverse;
	}

	/**
	 * returns string represantation of the orientation: +, - and +/- for null
	 * 
	 * @return + for forward strand, - for negative strand, +/- for variable not
	 *         set
	 */
	public String getStrandOrientation() {
		if (isReverse == null) {
			return "+/-";
		} else if (isReverse) {
			return "-";
		} else {
			return "+";
		}
	}
}
