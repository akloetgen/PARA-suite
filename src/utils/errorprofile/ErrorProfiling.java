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
 * represents the profile HMM to store data on error profile for particular
 * sequencing run
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
	}

	public void inferErrorProfile(boolean isInferQualities,
			boolean isShowErrorPlot) {
		try {
			File mappingFile = new File(mappingFileName);
			final SamReaderFactory factory = SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
					.validationStringency(ValidationStringency.LENIENT);
			SamReader samFileReader = factory.open(mappingFile);
			// File mappingFile = new File(mappingFileName);
			// SAMFileReader samFileReader = new SAMFileReader(mappingFile);

			IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(
					new File(referenceFileName));
			FileWriter fileWriterProfile = new FileWriter(mappingFile
					+ ".errorprofile");
			FileWriter fileWriterIndels = new FileWriter(mappingFile
					+ ".indels");
			FileWriter fileWriterIndelsMean = new FileWriter(mappingFile
					+ ".indelprofile");
			FileWriter fileWriterQualities = new FileWriter(mappingFile
					+ ".qualities");

			// try to use the iterator for establishing multi threading

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

			// MappingLogger.getLogger().debug("counting reads started");
			// int mappedReads = 0;
			// for (SAMRecord readHit : samFileReader) {
			// mappedReads++;
			// }
			// MappingLogger.getLogger().debug(
			// "counting reads finished; " + mappedReads
			// + " were counted.");

			// MappingLogger.getLogger().debug("separating reads started");
			//
			// // for (int z = 0; z < threads; z++) {
			// for (SAMRecord readHit : samFileReader) {
			//
			// for (int r = 0; r < (int) threads / mappedReads; r++) {
			// // will create 20 mio. SAMRecord objects -> not practicable
			// }
			//
			// }
			// MappingLogger.getLogger().debug("separating reads finished");

			// MAYBE RUN EACH READ ANALYSIS ON A SINGLE THREAD?!?! BUT PROBABLY
			// INVOKING THREADS COSTS MORE TIME THAN JUST LET IT RUN!??! MAYBE
			// TRUE FOR ONLY 2-3 THREADS?!?!

			// SO, CREATE 10 THREADS, GO THROUGH SAMFILEREADER AND START THE
			// NEXT FREE THREAD WITH THIS SAMRECROD, GO TO NEXT AND WAIT FOR THE
			// NEXT FREE THREAD AND SO ON; IN THE END, COMBINE ALL 10
			// TEMP_STATISTICS AS THEY ARE NOT WRITE-PROTECTED WHILE USING 10
			// THREADS!?!?!!!!!

			// QueryInterval[] intervals = new QueryInterval[threads];
			// SAMRecordIterator it = samFileReader.query(intervals, true);
			// while(it.hasNext()) {
			for (SAMRecord readHit : samFileReader) {
				// SAMRecord readHit = it.next();
				// for (samFileReader.iterator()) {

				// HERE IS AN IDEA!!!!
				// int useThreadIndexForTempStatistics = numReadsProcessed
				// % threads;

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
					MappingLogger.getLogger().debug(
							numReadsProcessed + " number reads processed");
					// break;
				}

				if (refSequenceForRead[0] == 0) {
					continue;
				}

				// introduce indels/exon-exon junctions; can be done before
				// rev-compl, because no mismatch information is needed and
				// cigar string is for the forward strand
				boolean skip = false;
				if (readSequence.length != refSequenceForRead.length) {
					// refSequenceForRead = reference.getSubsequenceAt(
					// readHit.getReferenceName(),
					// readHit.getAlignmentStart(),
					// readHit.getAlignmentEnd()).getBases();
					byte refSequenceForReadTemp[] = new byte[readSequence.length];
					// MappingLogger.getLogger().debug(
					// "NEXT: stradn: "
					// + readHit.getReadNegativeStrandFlag());
					// TODO
					// USE ALIGNMENT_BLOCKS, TAKE LENGHT + REF_START TO
					// RECONSTRUCT BOTH BYTE_ARRAYS AND USE -1/-2 FOR INDELS

					// System.out.println("read sequence vs. ref sequence");
					// for (int i = 0; i < readSequence.length; i++) {
					// System.out.print(readSequence[i]);
					// }
					// System.out.println();
					// for (int i = 0; i < refSequenceForRead.length; i++) {
					// System.out.print(refSequenceForRead[i]);
					// }
					// System.exit(0);
					int passed = 0;
					int passedMatches = 0;
					for (CigarElement elem : readHit.getCigar()
							.getCigarElements()) {

						// IF IT IS AN I OR D, THE CORRESPONDING BASE IS
						// MIISSING AND THIS IS SHIT WHEN CALCULATING MISMATCH
						// NUMBERS; HOW TO IGNORE THEMN LATER??!?! BY INDICATING
						// THAT WITH A NEW BYTE in the byte-array?!?!?!?!?

						if (elem.getOperator().toString()
								.equals(CigarOperator.M.toString())) {
							// MappingLogger.getLogger()
							// .debug("matches: "
							// + elem.getOperator().toString());
							// MappingLogger
							// .getLogger()
							// .debug("passedMatches: "
							// + passedMatches
							// + "; passed: "
							// + passed
							// + "; length elem: "
							// + elem.getLength()
							// + "; byte array size: "
							// + refSequenceForReadTemp.length
							// + "; CIGAR: "
							// + readHit.getCigarString()
							// + "; strand: "
							// + readHit
							// .getReadNegativeStrandFlag());

							for (int z = 0; z < elem.getLength(); z++) {
								try {
									refSequenceForReadTemp[z + passedMatches] = refSequenceForRead[z
											+ passed];
								} catch (ArrayIndexOutOfBoundsException e) {
									// SKIP!!! UNRESOLVED ERROR!!!!!!!!!!!!!!
									skip = true;
									// MappingLogger.getLogger()
									// .debug("CIGAR="
									// + readHit.getCigarString());
								}
								// MappingLogger.getLogger().debug(
								// "read: "
								// + readSequence[passedMatches
								// + z]
								// + "; refNew: "
								// + refSequenceForReadTemp[z
								// + passedMatches]);
							}
							passedMatches += elem.getLength();
							passed += elem.getLength();

							// } else {
							// MappingLogger.getLogger().debug(
							// "no matches: "
							// + readHit.getCigar()
							// .getCigarElement(z)
							// .getLength());
							// passed += elem.getLength();
						} else if (elem.getOperator().toString()
								.equals(CigarOperator.N.toString())) {
							passed += elem.getLength();
						} else if (elem.getOperator().toString()
								.equals(CigarOperator.I.toString())) {
							// read is LONGER than reference
							// MappingLogger.getLogger().debug(
							// "read length: " + readSequence.length
							// + "; ref length: "
							// + refSequenceForRead.length
							// + "; ref temp length: "
							// + refSequenceForReadTemp.length
							// + "; cigar length: "
							// + elem.getLength());
							// MappingLogger
							// .getLogger()
							// .debug("I match: passedMatches: "
							// + passedMatches
							// + "; passed: "
							// + passed
							// + "; length elem: "
							// + elem.getLength()
							// + "; byte array size: "
							// + refSequenceForReadTemp.length
							// + "; CIGAR: "
							// + readHit.getCigarString()
							// + "; strand: "
							// + readHit
							// .getReadNegativeStrandFlag());
							for (int z = 0; z < elem.getLength(); z++) {
								refSequenceForReadTemp[passedMatches + z] = -1;
							}
							passedMatches += elem.getLength();
							for (int q = 1; q <= elem.getLength(); q++) {
								insertionsPerPos[passedMatches + q] += 1.0;
							}
							if (elem.getLength() > 1) {
								longerIndels++;
							}

						} else if (elem.getOperator().toString()
								.equals(CigarOperator.D.toString())) {
							for (int q = 1; q <= elem.getLength(); q++) {
								deletionsPerPos[passedMatches + q] += 1.0;
							}

							if (elem.getLength() > 1) {
								longerIndels++;
							}
							// read is SHORTER than reference
							// MappingLogger.getLogger().debug(
							// "read length: " + readSequence.length
							// + "; ref length: "
							// + refSequenceForRead.length
							// + "; ref temp length: "
							// + refSequenceForReadTemp.length
							// + "; cigar length: "
							// + elem.getLength());
							// MappingLogger
							// .getLogger()
							// .debug("D match: passedMatches: "
							// + passedMatches
							// + "; passed: "
							// + passed
							// + "; length elem: "
							// + elem.getLength()
							// + "; byte array size: "
							// + refSequenceForReadTemp.length
							// + "; CIGAR: "
							// + readHit.getCigarString()
							// + "; strand: "
							// + readHit
							// .getReadNegativeStrandFlag());
						}

						// z++;
						// readHit.getCigar().getCigarElement(0);
					}
					// MappingLogger.getLogger().debug(
					// "read length: " + readSequence.length
					// + "; ref length: "
					// + refSequenceForRead.length + "; cigar: "
					// + readHit.getCigarString());
					// for (AlignmentBlock block : readHit.getAlignmentBlocks())
					// {
					// MappingLogger.getLogger().debug(
					// "block length: " + block.getLength());
					//
					// }

					indelRead++;
					refSequenceForRead = refSequenceForReadTemp;
					// System.exit(0);
					// skippedReads++;
					// continue;
				}

				// MappingLogger.getLogger().debug(
				// "readHit: " + readHit.getReadName());

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
					// rev compl vom Read wird im BAM-File gespeichert wenn
					// rev-strand!! Getestet an MSI1 PARCLIP Datensatz!!! :)

					htsjdk.samtools.util.SequenceUtil
							.reverseComplement(readSequence);
					htsjdk.samtools.util.SequenceUtil
							.reverseComplement(refSequenceForRead);

				}

				for (int i = 0; i < readSequence.length; i++) {
					int arrayPosRef = calculateArrayPos(refSequenceForRead[i]);
					int arrayPosRead = calculateArrayPos(readSequence[i]);

					// MappingLogger.getLogger().debug(
					// "read: " + readSequence[i] + "; ref: "
					// + refSequenceForRead[i]);

					if (arrayPosRef != arrayPosRead) {
						error++;
					}
					// in average it should be 75% errors if e.g. an exon was
					// wrongly annotated. as this can be the second half of the
					// read, thereby reducing this chance to 37.5% i would
					// postulate an error threshold of around 30% of the read
					// length.
					if (error >= (readSequence.length * 30 / 100)) {
						// skip these reads, because they originate from spliced
						// alignments that are 1 or 2 nucs shifted because of
						// wrong
						// exon annotations. maybe could be solver smarter OR
						// can be corrected?!?!??!
						skippedReads++;
						skip = true;
						// System.exit(0);
						break;
					}
				}
				if (skip) {
					continue;
				}
				// if (readHit.getReadName().equals(
				// "SEQ_ID:>ENST00000492267|3|82855666|82857172-2:6")) {
				// MappingLogger.getLogger().debug(
				// "FOUND: with errors: " + error);
				// }

				for (int i = 0; i < readSequence.length; i++) {
					int arrayPosRef = calculateArrayPos(refSequenceForRead[i]);
					int arrayPosRead = calculateArrayPos(readSequence[i]);

					if (arrayPosRef >= 0 && arrayPosRead >= 0) {
						positionConversions[i][arrayPosRef][arrayPosRead]++;
					}
					if (isInferQualities) {
						// um hier aus der liste ein array zu machen müsste die
						// gesamte Read-Anzahl vorher bekannt sein...
						baseQualitiesPerPos[i].add((int) readQualities[i]);
						baseQualitiesPerPosMean[i] += readQualities[i];
					}
				}

				// calculate SD for base qualities
				// if (isInferQualities) {
				// for (int i = 0; i < readSequence.length; i++) {
				// baseQualitiesPerPosSd[i]
				// }
				// }

			}
			double percentIndelReads = (double) indelRead / numReadsProcessed;
			MappingLogger.getLogger().debug(
					"Stats: unmapped: " + unmapped + "; failedVendor: "
							+ failedVendor + "; duplicate: " + duplicates
							+ "; startZero: " + startZero + "; %indelReads: "
							+ percentIndelReads
							+ "; skippedReads with >5 errors: " + skippedReads
							+ "; longer Indels: " + longerIndels);

			if (isInferQualities) {
				// fileWriterQualities.write("mean\tsd"
				// + System.getProperty("line.separator"));
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
						// wieso stand vorher total[i][j]???????
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
						// wieso war das drin??
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
			// double[][][] totalErrorCountsPerPosition = new
			// double[maxReadLength][4][4];
			// for (int i = 0; i < maxReadLength; i++) {
			// fileWriterProfile.write("@error rates"
			// + System.getProperty("line.separator"));
			for (int j = 0; j < 4; j++) {
				String errorCounts = "";
				for (int k = 0; k < 4; k++) {
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
			MappingLogger.getLogger().debug(
					"Averaged T2C = " + averageT2CEPR + " EPR");
			// }

			double insertionsOverall = 0.0;
			double deletionsOverall = 0.0;
			int insertionsZero = 0;
			int deletionsZero = 0;

			// fileWriterProfile.write("@indel rates"
			// + System.getProperty("line.separator"));
			for (int i = 0; i < maxReadLength; i++) {
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
				fileWriterIndels.write(insertionsPerPos[i] + "");
				fileWriterIndels.write("\t");
				fileWriterIndels.write(deletionsPerPos[i] + "");
				fileWriterIndels.write(System.getProperty("line.separator"));
			}
			insertionsOverall = (double) insertionsOverall
					/ (maxReadLength - insertionsZero);
			deletionsOverall = (double) deletionsOverall
					/ (maxReadLength - deletionsZero);
			fileWriterIndelsMean.write(insertionsOverall + "\t"
					+ deletionsOverall);

			/*
			 * Double insertionsRel = (double) insertions / numReadsProcessed;
			 * Double deletionsRel = (double) deletions / numReadsProcessed;
			 * fileWriterIndels.write(insertionsRel.toString());
			 * fileWriterIndels.write(System.getProperty("line.separator"));
			 * fileWriterIndels.write(deletionsRel.toString());
			 * fileWriterIndels.write(System.getProperty("line.separator"));
			 * MappingLogger.getLogger().debug("insertion rate: " +
			 * insertionsRel); MappingLogger.getLogger().debug("deletion rate: "
			 * + deletionsRel);
			 */

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

	// private String byteToNucleotide(int arrayPosition) {
	// String nucleotide = "";
	//
	// switch (arrayPosition) {
	// case 0:
	// nucleotide = "A";
	// break;
	// case 1:
	// nucleotide = "C";
	// break;
	// case 2:
	// nucleotide = "G";
	// break;
	// case 3:
	// nucleotide = "T";
	// break;
	// }
	//
	// return nucleotide;
	// }

}
