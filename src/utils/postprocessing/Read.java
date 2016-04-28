package utils.postprocessing;

import htsjdk.samtools.SAMRecord;

import java.util.LinkedList;
import java.util.List;

/**
 * Container for a single sequencing read. Saves additional information
 * necessary for calculating the mapping position in the reference genome for
 * reads aligned against a transcriptome sequence.
 * 
 * @author akloetgen
 * 
 */
public class Read {

	private String readName;
	private int primaryIndex;
	private List<String> genesHitted;
	private List<Integer> newGenomicPositions;
	private List<String> newGenomicCigars;
	private List<SAMRecord> records;

	public Read(String readName) {
		this.readName = readName;
		genesHitted = new LinkedList<String>();
		newGenomicPositions = new LinkedList<Integer>();
		newGenomicCigars = new LinkedList<String>();
		records = new LinkedList<SAMRecord>();
	}

	public String getReadName() {
		return readName;
	}

	public void setReadName(String readName) {
		this.readName = readName;
	}

	public List<String> getGenesHitted() {
		return genesHitted;
	}

	public void addGeneHit(String gene) {
		genesHitted.add(gene);
	}

	public List<Integer> getNewGenomicPositions() {
		return newGenomicPositions;
	}

	public void addHitToGenomicPosition(Integer genomicPosition) {
		newGenomicPositions.add(genomicPosition);
	}

	public List<SAMRecord> getRecords() {
		return records;
	}

	public void addRecord(SAMRecord newRecord) {
		records.add(newRecord);
	}

	public int getPrimaryIndex() {
		return primaryIndex;
	}

	public void setPrimaryIndex(int primaryIndex) {
		this.primaryIndex = primaryIndex;
	}

	public List<String> getNewGenomicCigars() {
		return newGenomicCigars;
	}

	public void addHitToGenomicCigars(String newGenomicCigar) {
		newGenomicCigars.add(newGenomicCigar);
	}

}
