package utils.postprocessing;

import htsjdk.samtools.CigarOperator;

public class IndelTupel {

	private int indelPosition;
	private CigarOperator indelType;

	public IndelTupel(int indelPosition, CigarOperator indelType) {
		this.setIndelPosition(indelPosition);
		this.setIndelType(indelType);
	}

	public int getIndelPosition() {
		return indelPosition;
	}

	public void setIndelPosition(int indelPosition) {
		this.indelPosition = indelPosition;
	}

	public CigarOperator getIndelType() {
		return indelType;
	}

	public void setIndelType(CigarOperator indelType) {
		this.indelType = indelType;
	}

}
