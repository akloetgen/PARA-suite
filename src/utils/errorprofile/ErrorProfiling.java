package utils.errorprofile;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFrame;

import main.MappingLogger;

import org.math.plot.Plot2DPanel;

/**
 * 
 * Represents the error profile for particular a sequencing run
 * 
 * @author akloetgen
 * 
 */
public class ErrorProfiling {

	private String mappingFileName;
	private String referenceFileName;
	private int maxReadLength;
	// ?????????????????

	// if int gets too small, use long (double precision as int)
	private int[][][] positionConversions;
	private int[][] totalBaseCountsPerPos;
	private double[][][] conversionsPercentagePerPos;
	private double[][] totalErrorCounts;
	private List<Integer>[] baseQualitiesPerPos;
	private double[] baseQualitiesPerPosMean;
	private double[] baseQualitiesPerPosSd;
	private double[] insertionsPerPos;
	private double[] deletionsPerPos;
	private int[] totalCountsPerPos;
	private String[] basePositions;
	private int totalBasesChecked;
	private int totalBasesMissed;

	@SuppressWarnings("unchecked")
	public ErrorProfiling(String mappingFileName, String referenceFileName,
			int maxReadLength) {
		this.mappingFileName = mappingFileName;
		this.referenceFileName = referenceFileName;
		this.maxReadLength = maxReadLength;

		positionConversions = new int[this.maxReadLength][4][4];
		totalBaseCountsPerPos = new int[this.maxReadLength][4];
		conversionsPercentagePerPos = new double[this.maxReadLength][4][4];
		totalErrorCounts = new double[4][4];
		baseQualitiesPerPosMean = new double[this.maxReadLength];
		baseQualitiesPerPosSd = new double[this.maxReadLength];
		baseQualitiesPerPos = new List[this.maxReadLength];
		for (int i = 0; i < baseQualitiesPerPos.length; i++) {
			List<Integer> newList = new LinkedList<Integer>();
			baseQualitiesPerPos[i] = newList;
		}
		insertionsPerPos = new double[this.maxReadLength];
		deletionsPerPos = new double[this.maxReadLength];
		totalCountsPerPos = new int[this.maxReadLength];
		basePositions = new String[4];
		basePositions[0] = "A";
		basePositions[1] = "C";
		basePositions[2] = "G";
		basePositions[3] = "T";
		totalBasesChecked = 0;
		totalBasesMissed = 0;
	}

	/**
	 * Infers the error profile of a given sequencing run.
	 * 
	 * @param isInferQualities
	 *            whether base calling qualities should be measured
	 * @param isShowErrorPlot
	 *            whether error profile should be plotted
	 */
	public void inferErrorProfile(boolean isInferQualities,
			boolean isShowErrorPlot) {
		try {
			File mappingFile = new File(mappingFileName);
			final SamReaderFactory factory = SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
					.validationStringency(ValidationStringency.LENIENT);
			SamReader samFileReader = factory.open(mappingFile);

			IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(
					new File(referenceFileName));
			FileWriter fileWriterProfile = new FileWriter(mappingFile
					+ ".errorprofile");
			FileWriter fileWriterVCFProfile = new FileWriter(mappingFile
					+ ".errorprofile.vcf");
			FileWriter fileWriterIndels = new FileWriter(mappingFile
					+ ".indels");
			FileWriter fileWriterIndelsMean = new FileWriter(mappingFile
					+ ".indelprofile");
			FileWriter fileWriterQualities = new FileWriter(mappingFile
					+ ".qualities");

			if (!samFileReader.getFileHeader().getSortOrder()
					.equals(SAMFileHeader.SortOrder.coordinate)) {
				MappingLogger.getLogger().debug(
						samFileReader.getFileHeader().getSortOrder());
				MappingLogger.getLogger().error(
						"BAM file is not sorted. Please provide "
								+ "a sorted BAM-file as input alignment file.");
				System.exit(0);
			}

			int numReadsProcessed = 0;
			int unmapped = 0;
			int failedVendor = 0;
			int duplicates = 0;
			int startZero = 0;
			int indelRead = 0;
			int skippedReads = 0;
			int longerIndels = 0;

			boolean isFilter = false;

			// multi-threading missing. Still very slow...
			for (SAMRecord readHit : samFileReader) {
				// StringBuffer readSeq = new StringBuffer();
				// for (int basePosition = 0; basePosition < readHit
				// .getReadBases().length; basePosition++) {
				// readSeq.append(Character.toUpperCase((char) readHit
				// .getReadBases()[basePosition]));
				// }
				// System.out.println("readHit=" + readHit.getReadName()
				// + " with seq=" + readSeq);
				if (readHit.getReadUnmappedFlag()) {
					unmapped++;
					continue;
				}
				if (readHit.getDuplicateReadFlag()) {
					duplicates++;
					continue;
				}
				if (readHit.getAlignmentStart() == 0) {
					startZero++;
					continue;
				}

				byte readSequence[] = readHit.getReadBases();
				byte refSequenceForRead[] = reference.getSubsequenceAt(
						readHit.getReferenceName(),
						readHit.getAlignmentStart(), readHit.getAlignmentEnd())
						.getBases();

				numReadsProcessed++;
				if (numReadsProcessed % 250000 == 0) {
					MappingLogger.getLogger().info(
							numReadsProcessed + " number reads processed");
				}

				if (refSequenceForRead[0] == 0) {
					continue;
				}

				// introduce indels/exon-exon junctions; can be done before
				// rev-compl, because no mismatch information is needed and
				// cigar string is for the forward strand
				boolean skip = false;
				int mappingLength;
				if (readSequence.length > refSequenceForRead.length) {
					mappingLength = readSequence.length;
				} else {
					mappingLength = refSequenceForRead.length;
				}
				if (readSequence.length != refSequenceForRead.length) {
					byte refSequenceForReadTemp[] = new byte[mappingLength];
					byte readSequenceTemp[] = new byte[mappingLength];

					// TODO
					// USE ALIGNMENT_BLOCKS, TAKE LENGHT + REF_START TO
					// RECONSTRUCT BOTH BYTE_ARRAYS AND USE -1/-2 FOR INDELS

					int passed = 0;
					int passedRef = 0;
					int passedRead = 0;
					int passedMatches = 0;
					for (CigarElement elem : readHit.getCigar()
							.getCigarElements()) {

						// IF IT IS AN I OR D, THE CORRESPONDING BASE IS
						// MIISSING AND THIS IS SHIT WHEN CALCULATING MISMATCH
						// NUMBERS; HOW TO IGNORE THEMN LATER??!?! BY INDICATING
						// THAT WITH A NEW BYTE in the byte-array?!?!?!?!?

						if (elem.getOperator().equals(CigarOperator.M)
								|| elem.getOperator().equals(CigarOperator.X)
								|| elem.getOperator().equals(CigarOperator.EQ)) {
							for (int z = 0; z < elem.getLength(); z++) {
								try {
									refSequenceForReadTemp[z + passedMatches] = refSequenceForRead[z
											+ passedRef];
									// if (readHit.getCigar().getCigarElements()
									// .contains(CigarOperator.D)) {
									// System.out.println("start="
									// + readHit.getAlignmentStart()
									// + " and end="
									// + readHit.getAlignmentEnd());
									readSequenceTemp[z + passedMatches] = readSequence[z
											+ passedRead];
									// }
								} catch (ArrayIndexOutOfBoundsException e) {
									// SKIP!!! UNRESOLVED ERROR!!!!!!!!!!!!!!
									// i think this error is resolved, but
									// should be tested...
									MappingLogger.getLogger().error(
											"random error. sorry.");
									// e.printStackTrace();
									skip = true;
								}
							}
							passedMatches += elem.getLength();
							passedRef += elem.getLength();
							passedRead += elem.getLength();
						} else if (elem.getOperator().toString()
								.equals(CigarOperator.N.toString())) {
							// passed += elem.getLength();
							passedRef += elem.getLength();
							passedRead += elem.getLength();
						} else if (elem.getOperator().toString()
								.equals(CigarOperator.I.toString())) {
							// read is LONGER than reference

							for (int z = 0; z < elem.getLength(); z++) {
								refSequenceForReadTemp[passedMatches + z] = 45;
								// readSequenceTemp[passedMatches + z] =
								// readSequence[passed
								// + z];
							}
							passedMatches += elem.getLength();
							passedRead += elem.getLength();

							for (int q = 1; q <= elem.getLength(); q++) {
								insertionsPerPos[passedMatches + q] += 1.0;
							}
							if (elem.getLength() > 1) {
								longerIndels++;
							}

						} else if (elem.getOperator().toString()
								.equals(CigarOperator.D.toString())) {

							// read is SHORTER than reference
							for (int z = 0; z < elem.getLength(); z++) {
								readSequenceTemp[passedMatches + z] = 45;
								// refSequenceForReadTemp[passedMatches + z] =
								// refSequenceForRead[passed
								// + z];
							}
							passedMatches += elem.getLength();
							passedRef += elem.getLength();
							// passed += elem.getLength();

							for (int q = 1; q <= elem.getLength(); q++) {
								deletionsPerPos[passedMatches + q] += 1.0;
							}

							if (elem.getLength() > 1) {
								longerIndels++;
							}

						}
					}
					indelRead++;
					refSequenceForRead = refSequenceForReadTemp;
					readSequence = readSequenceTemp;
				}

				byte[] readQualities = readHit.getBaseQualities();
				int error = 0;
				if (skip) {
					skippedReads++;
					continue;
				}
				skip = false;

				// create rev-compl if negative strand flag is set to calculate
				// correct mismatches
				if (readHit.getReadNegativeStrandFlag()) {
					htsjdk.samtools.util.SequenceUtil
							.reverseComplement(readSequence);
					htsjdk.samtools.util.SequenceUtil
							.reverseComplement(refSequenceForRead);
				}

				for (int i = 0; i < readSequence.length; i++) {
					int arrayPosRef = calculateArrayPos(refSequenceForRead[i]);
					int arrayPosRead = calculateArrayPos(readSequence[i]);

					if (arrayPosRef != arrayPosRead) {
						error++;
					}
					// in average it should be 75% errors if e.g. an exon was
					// wrongly annotated. as this can be the second half of the
					// read, thereby reducing this chance to 37.5% i would
					// postulate an error threshold of around 30% of the read
					// length.
					if (error >= (readSequence.length * 30 / 100) && isFilter) {
						// skip these reads, because they originate from spliced
						// alignments that are 1 or 2 nucs shifted because of
						// wrong exon annotations. maybe could be solved smarter
						// OR can be corrected?!?!??!

						// i found a few examples in the ensembl genes v75 db
						// where wrong annotations led to those errors. Maybe
						// resolved on hg38 and Ensembl genes v78?
						skippedReads++;
						skip = true;
						break;
					}
				}
				if (skip) {
					continue;
				}
				for (int i = 0; i < readSequence.length; i++) {
					int arrayPosRef = calculateArrayPos(refSequenceForRead[i]);
					int arrayPosRead = calculateArrayPos(readSequence[i]);

					int refPos = -1;
					if (readHit.getReadNegativeStrandFlag()) {
						refPos = i + readHit.getAlignmentStart();
					} else {
						refPos = readHit.getAlignmentEnd() - i;
					}

					StringBuffer readSeq = new StringBuffer();
					for (int basePosition = 0; basePosition < readSequence.length; basePosition++) {
						readSeq.append(Character
								.toUpperCase((char) readSequence[basePosition]));
					}
					StringBuffer refSeq = new StringBuffer();
					for (int basePosition = 0; basePosition < refSequenceForRead.length; basePosition++) {
						refSeq.append(Character
								.toUpperCase((char) refSequenceForRead[basePosition]));
					}
//					if (arrayPosRef == 0 && arrayPosRead == 2) {
//						System.out.println("A to G hit in read="
//								+ readHit.getReadName() + " with seq="
//								+ readSeq + " and ref-seq=" + refSeq);
//					}

					if (arrayPosRef >= 0 && arrayPosRead >= 0) {
						positionConversions[i][arrayPosRef][arrayPosRead]++;
						totalBasesChecked++;
					} else {
						// System.out.println("wrong bases: REF=" + arrayPosRef
						// + " (" + refSequenceForRead[i]
						// + " at position " + refPos + ")" + "\tREAD="
						// + arrayPosRead + " (" + readSequence[i]
						// + " at position " + i + ")" + " for read "
						// + readHit.getReadName() + " orient="
						// + readHit.getReadNegativeStrandFlag());
						// totalBasesMissed++;
					}
					if (isInferQualities) {
						// um hier aus der liste ein array zu machen müsste die
						// gesamte Read-Anzahl vorher bekannt sein...
						baseQualitiesPerPos[i].add((int) readQualities[i]);
						baseQualitiesPerPosMean[i] += readQualities[i];
					}
				}
			}
			double percentIndelReads = (double) indelRead / numReadsProcessed;
			MappingLogger.getLogger().debug(
					"Stats: unmapped: " + unmapped + "; failedVendor: "
							+ failedVendor + "; duplicate: " + duplicates
							+ "; startZero: " + startZero + "; %indelReads: "
							+ percentIndelReads
							+ "; skippedReads with >5 errors: " + skippedReads
							+ "; longer Indels: " + longerIndels
							+ "; totalBasesChecked=" + totalBasesChecked
							+ "; totalBasesMissed=" + totalBasesMissed);

			if (isInferQualities) {
				for (int i = 0; i < baseQualitiesPerPos.length; i++) {
					baseQualitiesPerPosMean[i] = (double) baseQualitiesPerPosMean[i]
							/ baseQualitiesPerPos[i].size();
					double tempSdValue = 0.0;
					for (Integer val : baseQualitiesPerPos[i]) {
						tempSdValue += Math.pow(
								(val - baseQualitiesPerPosMean[i]), 2);
					}
					baseQualitiesPerPosSd[i] = Math.sqrt((double) tempSdValue
							/ baseQualitiesPerPos[i].size());

					fileWriterQualities.write(baseQualitiesPerPosMean[i] + "\t"
							+ baseQualitiesPerPosSd[i]
							+ System.getProperty("line.separator"));
				}
			}

			double[] totalBaseCounts = new double[4];
			// double[][] totalBaseCountsPerPos = new double[maxReadLength][4];
			for (int i = 0; i < maxReadLength; i++) {
				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 4; k++) {
						totalErrorCounts[j][k] += positionConversions[i][j][k];
						totalBaseCounts[j] += positionConversions[i][j][k];
						totalCountsPerPos[i] += positionConversions[i][j][k];
						totalBaseCountsPerPos[i][j] += positionConversions[i][j][k];
					}
				}
			}

			double[] mutationsT2C = new double[maxReadLength];
			double[] averageMut = new double[maxReadLength];
			double[] averageMutExceptT2C = new double[maxReadLength];
			for (int i = 0; i < maxReadLength; i++) {
				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 4; k++) {
						// ???
						// if (totalCountsPerPosPerBase[i][j] == 0) {
						// continue;
						// }
						conversionsPercentagePerPos[i][j][k] = (double) positionConversions[i][j][k]
								/ numReadsProcessed * 100;

					}
				}

				mutationsT2C[i] = (double) conversionsPercentagePerPos[i][3][1];
				averageMut[i] = (double) (conversionsPercentagePerPos[i][0][1]
						+ conversionsPercentagePerPos[i][0][2]
						+ conversionsPercentagePerPos[i][0][3]
						+ conversionsPercentagePerPos[i][1][0]
						+ conversionsPercentagePerPos[i][1][2]
						+ conversionsPercentagePerPos[i][1][3]
						+ conversionsPercentagePerPos[i][2][0]
						+ conversionsPercentagePerPos[i][2][1]
						+ conversionsPercentagePerPos[i][2][3]
						+ conversionsPercentagePerPos[i][3][0]
						+ conversionsPercentagePerPos[i][3][1] + conversionsPercentagePerPos[i][3][2]);
				averageMutExceptT2C[i] = (double) (conversionsPercentagePerPos[i][0][1]
						+ conversionsPercentagePerPos[i][0][2]
						+ conversionsPercentagePerPos[i][0][3]
						+ conversionsPercentagePerPos[i][1][0]
						+ conversionsPercentagePerPos[i][1][2]
						+ conversionsPercentagePerPos[i][1][3]
						+ conversionsPercentagePerPos[i][2][0]
						+ conversionsPercentagePerPos[i][2][1]
						+ conversionsPercentagePerPos[i][2][3]
						+ conversionsPercentagePerPos[i][3][0] + conversionsPercentagePerPos[i][3][2]);

				// calculate error-profile for each base (16 values: A-A, A-C,
				// A-G, A-T, C-A, C-C, ..., T-T)

			}

			for (int j = 0; j < 4; j++) {
				String errorCounts = "";
				for (int k = 0; k < 4; k++) {
					fileWriterVCFProfile.write(basePositions[j] + "\t"
							+ basePositions[k] + "\t" + totalErrorCounts[j][k]);
					fileWriterVCFProfile.write(System
							.getProperty("line.separator"));
					totalErrorCounts[j][k] = (double) totalErrorCounts[j][k]
							/ totalBaseCounts[j];
					errorCounts += totalErrorCounts[j][k] + "\t";
					fileWriterProfile.write(totalErrorCounts[j][k] + "\t");

					// totalErrorCountsPerPosition[i][j][k] = (double)
					// positionConversions[i][j][k]
					// / totalBaseCountsPerPos[i][j];
					// fileWriterProfile
					// .write(totalErrorCountsPerPosition[i][j][k]
					// + "\t");
				}
				MappingLogger.getLogger().debug(errorCounts);
				fileWriterProfile.write(System.getProperty("line.separator"));
				fileWriterVCFProfile
						.write(System.getProperty("line.separator"));
			}
			double averageT2CEPR = 0.0;
			for (int j = 0; j < mutationsT2C.length; j++) {
				if (mutationsT2C[j] > 0) {
					averageT2CEPR += mutationsT2C[j];
				} else {
					averageT2CEPR = (double) averageT2CEPR / (j + 1);
					break;
				}
			}
			MappingLogger.getLogger().info(
					"Averaged T2C = " + averageT2CEPR + " EPR");
			// }

			double insertionsOverall = 0.0;
			double deletionsOverall = 0.0;
			int insertionsZero = 0;
			int deletionsZero = 0;

			// to avoid NaN values if no insertion/deletion was set during
			// alignment process, just set indel probabilites to 0.

			for (int i = 0; i < maxReadLength; i++) {
				if (totalCountsPerPos[i] == 0.0) {
					insertionsPerPos[i] = 0.0;
					deletionsPerPos[i] = 0.0;
					insertionsZero++;
					deletionsZero++;
				} else {
					insertionsPerPos[i] = (double) insertionsPerPos[i]
							/ totalCountsPerPos[i];
					if (insertionsPerPos[i] > 0) {
						insertionsOverall += insertionsPerPos[i];
					} else {
						insertionsZero++;
					}
					deletionsPerPos[i] = (double) deletionsPerPos[i]
							/ totalCountsPerPos[i];
					if (deletionsPerPos[i] > 0) {
						deletionsOverall += deletionsPerPos[i];
					} else {
						deletionsZero++;
					}
				}
				fileWriterIndels.write(insertionsPerPos[i] + "");
				fileWriterIndels.write("\t");
				fileWriterIndels.write(deletionsPerPos[i] + "");
				fileWriterIndels.write(System.getProperty("line.separator"));
			}
			if (maxReadLength == insertionsZero
					&& maxReadLength == deletionsZero) {
				insertionsOverall = 0.0;
				deletionsOverall = 0.0;
			} else {
				insertionsOverall = (double) insertionsOverall
						/ (maxReadLength - insertionsZero);
				deletionsOverall = (double) deletionsOverall
						/ (maxReadLength - deletionsZero);
			}
			fileWriterIndelsMean.write(insertionsOverall + "\t"
					+ deletionsOverall);

			if (isShowErrorPlot) {
				Plot2DPanel plot = new Plot2DPanel();

				plot.addLinePlot("T2C errors", Color.BLUE, mutationsT2C);
				plot.addLinePlot("All errors", Color.RED, averageMut);
				plot.addLinePlot("All errors except T2C", Color.GREEN,
						averageMutExceptT2C);
				plot.setAxisLabel(0, "position within read");
				plot.setAxisLabel(1, "Errors Per Reads (EPR) (* 100)");
				// plot.setFixedBounds(1, 0, 0.2);
				plot.setFixedBounds(0, 0, maxReadLength);
				plot.addLegend("WEST");

				JFrame frame = new JFrame("Error profile for file: "
						+ mappingFileName);
				frame.setSize(1000, 450);
				frame.setLocation(100, 100);
				frame.setContentPane(plot);
				frame.setVisible(true);
			}

			samFileReader.close();
			reference.close();
			fileWriterProfile.close();
			fileWriterVCFProfile.close();
			fileWriterIndels.close();
			fileWriterQualities.close();
			fileWriterIndelsMean.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			MappingLogger.getLogger().error(
					"file not found; probably index for genome-fasta missing?");
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private byte calculateArrayPos(byte formerByte) {
		byte arrayPosForByte = -1;

		switch (formerByte) {
		case 65:
			arrayPosForByte = 0;
			break;
		case 67:
			arrayPosForByte = 1;
			break;
		case 71:
			arrayPosForByte = 2;
			break;
		case 84:
			arrayPosForByte = 3;
			break;
		case 97:
			arrayPosForByte = 0;
			break;
		case 99:
			arrayPosForByte = 1;
			break;
		case 103:
			arrayPosForByte = 2;
			break;
		case 116:
			arrayPosForByte = 3;
			break;
		}

		return arrayPosForByte;
	}
}
