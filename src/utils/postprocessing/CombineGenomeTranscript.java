package utils.postprocessing;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import main.MappingLogger;

public class CombineGenomeTranscript {

	private SAMFileHeader genomeSamFileHeader;
	private SAMFileWriter combinedResults;

	private int mappedReads;
	private int splicedReads;

	public CombineGenomeTranscript() {
		mappedReads = 0;
		splicedReads = 0;
	}

	public void combine(String genomeMappingFileName,
			String transcriptMappingFileName, String combinedFileName) {
		try {
			File genomeMappingFile = new File(genomeMappingFileName);
			final SamReaderFactory factory = SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
					.validationStringency(ValidationStringency.LENIENT);
			SamReader genomeSamFileReader = factory.open(genomeMappingFile);
			// File genomeMappingFile = new File(genomeMappingFileName);
			// SAMFileReader genomeSamFileReader = new
			// SAMFileReader(genomeMappingFile);
			SAMFileWriterFactory fac = new SAMFileWriterFactory();
			combinedResults = fac.makeBAMWriter(genomeSamFileReader
					.getFileHeader(), false, new File(combinedFileName));
			int missedTranscriptAlignments = 0;

			for (SAMRecord readHit : genomeSamFileReader) {
				combinedResults.addAlignment(readHit);
				// if (readHit.getReadName().equals(
				// "HWI-ST737:305:C232HACXX:6:2106:2427:18829")) {
				// MappingLogger.getLogger().debug("übereinstimmung genomic.");
				// }
				mappedReads++;
			}
			genomeSamFileHeader = genomeSamFileReader.getFileHeader();

			genomeSamFileReader.close();

			File transcriptMappingFile = new File(transcriptMappingFileName);
			// final SamReaderFactory factoryTranscript =
			// SamReaderFactory.makeDefault()
			// .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
			// .validationStringency(ValidationStringency.LENIENT);
			SamReader transcriptSamFileReader = factory
					.open(transcriptMappingFile);
			// File transcriptMappingFile = new File(transcriptMappingFileName);
			// SAMFileReader transcriptSamFileReader = new SAMFileReader(
			// transcriptMappingFile);

			// save info on: read; read-hit; gene; genomic relocation: after
			// all,
			// check whether one read has several hits and in that case check
			// whether it is always the same gene and (nearly) same genomic
			// position. accept in that case.

			// create new read-hit for genomic stuff and add it to genomic
			// mapping
			// file. re-sort.

			if (!transcriptSamFileReader.getFileHeader().getSortOrder()
					.equals(SAMFileHeader.SortOrder.queryname)) {
				MappingLogger.getLogger().debug(
						transcriptSamFileReader.getFileHeader().getSortOrder());
				MappingLogger.getLogger().error(
						"BAM file " + transcriptMappingFileName
								+ " is not sorted. Please provide "
								+ "a sorted BAM-file as input alignment file.");
				System.exit(0);
			}

			String readNameTemp = "";
			HashMap<String, Read> readIDToRead = new HashMap<String, Read>();

			// CHECK FOR NAME-SORTED, THEN GO THROUGH SAM-RECORDS AND PRINT IF
			// NEXT
			// READ-NAME, DONT FOREGT LAST GROUP OF READ-HITS. WILL OVERCOME OUT
			// OF
			// MEMORY ERROR!!!
			for (SAMRecord readHit : transcriptSamFileReader) {
				if (readHit.getReferenceName().equals("*")) {
					continue;
				}

				List<IndelTupel> indelTupelList = new LinkedList<IndelTupel>();
				int cigarLengths = 0;
				for (CigarElement cigarElem : readHit.getCigar()
						.getCigarElements()) {
					if (cigarElem.getOperator().equals(CigarOperator.D)
							|| cigarElem.getOperator().equals(
									CigarOperator.DELETION)
							|| cigarElem.getOperator().equals(CigarOperator.I)
							|| cigarElem.getOperator().equals(
									CigarOperator.INSERTION)) {
						// position maybe +1??
						indelTupelList.add(new IndelTupel(cigarLengths + 1,
								cigarElem.getOperator()));
					}
					cigarLengths += cigarElem.getLength();
				}
				for (IndelTupel indel : indelTupelList) {
//					MappingLogger.getLogger().debug(
//							"position: " + indel.getIndelPosition()
//									+ " and type: "
//									+ indel.getIndelType().toString());
				}

				// if (readHit.getReadName().equals(
				// "HWI-ST737:305:C232HACXX:6:2106:2427:18829")) {
				// MappingLogger.getLogger().debug("übereinstimmung transcript");
				// }

				if (!readNameTemp.equals(readHit.getReadName())) {
					// check for multiple re-locations, otherwise add to
					// combine.bam
					// file
					printReadsToBamFile(readIDToRead);
					readIDToRead = new HashMap<String, Read>();
					readNameTemp = readHit.getReadName();
				}

				String refName = readHit.getReferenceName();
				// MappingLogger.getLogger().debug(refName);
				String[] readNameSplitted = refName.split("\\|");

				// MappingLogger.getLogger().debug(
				// "align-start: " + read.getAlignmentStart()
				// + "; align-end: " + read.getAlignmentEnd());

				// check second best mappings for same gene as primary mapping?
				// in
				// that case, accept primary mapping as good although it could
				// have
				// low MAPQ.

				// reallocate genomic position of transcript hits and add them
				// to
				// genomic mapping-BAM file
				String exonStartsString = readNameSplitted[3];
				// MappingLogger.getLogger().debug(exonStartsString);
				String exonEndsString = readNameSplitted[4];
				// MappingLogger.getLogger().debug(exonEndsString);
				String transcriptStrand = readNameSplitted[5];

				// exon positions of ensembl seem to be 1-based
				String[] exonStarts = exonStartsString.split(";");
				String[] exonEnds = exonEndsString.split(";");
				Arrays.sort(exonStarts);
				Arrays.sort(exonEnds);
				Integer newGenomicPositionStart = -1;
				String newGenomicCigar = "";

				// MappingLogger.getLogger().debug(
				// readHit.getAlignmentStart() + "-"
				// + readHit.getAlignmentEnd() + "; le: "
				// + readHit.getReadLength());
				//

				// if (readHit.getReadName().equals(
				// "SEQ_ID:>ENST00000488297|5|87677774|87678542-3:9")) {
				// if (transcriptStrand.equals("-1")) {
				// for (int i = 0; i < exonStarts.length; i++) {
				// MappingLogger.getLogger().debug(
				// "exon " + i + ": " + exonStarts[i] + " - "
				// + exonEnds[i]);
				// }
				// }
				int lengthPassed = 0;
				int matchingPositionsPassed = 0;
				boolean spliceJunctionCovered = false;
				// int convertedBases = 0;
				// int currentCigarElemIdx = 0;
				// int passedCigarElemLength = 0;
				// List<CigarElement> cigarElems = readHit.getCigar()
				// .getCigarElements();

				// CHECK WHETHER CIGAR CONTAINS "I"/"D" AND USE ANOTHER CIGAR
				// STRING
				// BUILDING FUNCTION!!!!!!!!!!!!!!!!

				// PRINT OUT 2-3 "I"/"D"s THAT ARE SPANNING (INDICATED BY "N",
				// TAKE
				// LOCATION AND PREVIOUS CIGAR, THEN TRY TO OPTIMIZE WITH THOSE
				// TEST-CASES)

				// works also on unstranded protocols
				if (transcriptStrand.equals("1")) {
					// foward

					for (int i = 0; i < exonStarts.length; i++) {
						int lengthPassedTemp = lengthPassed;
						lengthPassed += Integer.parseInt(exonEnds[i])
								- Integer.parseInt(exonStarts[i]) + 1;

						if (readHit.getAlignmentStart() <= lengthPassed) {
							if (newGenomicPositionStart == -1) {
								newGenomicPositionStart = Integer
										.parseInt(exonStarts[i])
										+ (readHit.getAlignmentStart() - lengthPassedTemp)
										- 1;
								// BAM IS 1-BASED!!!!
								// -1;
							}
						}

						if (readHit.getAlignmentEnd() <= lengthPassed) {
							if (newGenomicPositionStart >= Integer
									.parseInt(exonStarts[i])) {
								// read maps to 1 exon continuously --- MAYBE
								// EXCHANGE THIS WITH THE CIGAR OF THE MAPPING
								// ITSELF TO AUTOMATICALLY INCLUDE
								// INDELS?!?!?!?!?!?!
								newGenomicCigar = readHit.getCigarString();
								// newGenomicCigar = (readHit.getReadLength()) +
								// "M";
							} else {
								// adding Ms from exon-intron boundary to
								// alignment
								// end ---- USE THE CIGAR OF THE MAPPING ITSELF
								// TO
								// INCLUDE INDELS?!?!?!
								// if (true) {
								// continue;
								// }
								// if (readHit.getCigarString().contains(
								// CigarOperator.DELETION.toString())
								// || readHit.getCigarString().contains(
								// CigarOperator.INSERTION
								// .toString())) {
								// if (!spliceJunctionCovered) {
								// // just copy transcriptomic cigar until
								// it is possible to the current position,
								// inculding the indel
								// int cigarElemLenghts = 0;
								// for (CigarElement elem :
								// readHit.getCigar().getCigarElements()) {
								// cigarElemLenghts += elem.getLength();
								// if (cigarElemLenghts >
								// matchingPositionsPassed) {
								//
								// } else {
								// newGenomicCigar += elem.toString();
								// }
								// }
								// } else {
								// // recalculate the actual indel position
								// AFTER the exon-exon junction!
								// }
								// missedTranscriptAlignments++;
								// break;
								// }

								// newGenomicCigar += readHit.getCigarString()
								// .substring(convertedBases);

//								MappingLogger.getLogger().debug(
//										"matching positions passed before: "
//												+ matchingPositionsPassed);
								newGenomicCigar += (readHit.getAlignmentEnd() - lengthPassedTemp)
										+ "M";
								matchingPositionsPassed += (readHit
										.getAlignmentEnd() - lengthPassedTemp);
//								MappingLogger.getLogger().debug(
//										"matching positions passed after: "
//												+ matchingPositionsPassed);
							}
							break;
						} else if (newGenomicPositionStart != -1) {
							// if (true) {
							// continue;
							// }
							if (readHit.getCigarString().contains(
									CigarOperator.DELETION.toString())
									|| readHit.getCigarString().contains(
											CigarOperator.INSERTION.toString())) {
								missedTranscriptAlignments++;
								break;
							}
							if (newGenomicPositionStart >= Integer
									.parseInt(exonStarts[i])) {

								// from alignment start to exon-intron boundary
								// convertedBases +=
								// Integer.parseInt(exonEnds[i])
								// - newGenomicPositionStart + 1;

								// while (true) { // ÄNDERN AUF ABBRUCH!
								// if (passedCigarElemLength
								// + cigarElems.get(currentCigarElemIdx)
								// .getLength() < convertedBases) {
								// newGenomicCigar += cigarElems
								// .get(currentCigarElemIdx);
								// passedCigarElemLength += cigarElems.get(
								// currentCigarElemIdx).getLength();
								// currentCigarElemIdx++;
								// } else {
								// newGenomicCigar += (convertedBases -
								// passedCigarElemLength)
								// + cigarElems
								// .get(currentCigarElemIdx)
								// .getOperator().toString();
								// break;
								// }
								// }

//								MappingLogger.getLogger().debug(
//										"matching positions passed before: "
//												+ matchingPositionsPassed);
								matchingPositionsPassed += (Integer
										.parseInt(exonEnds[i])
										- newGenomicPositionStart + 1);
								newGenomicCigar += (Integer
										.parseInt(exonEnds[i])
										- newGenomicPositionStart + 1) + "M";

//								MappingLogger.getLogger().debug(
//										"matching positions passed after: "
//												+ matchingPositionsPassed);
							} else {
								// entire exon covered by read:
								// KOMPLIZIERT?!?! VLLT. INT MITLAUFEN LASSEN UM
								// START POS IM READ CIGAR ZU
								// FINDEN?!?!?!?!??!?!?!

//								MappingLogger.getLogger().debug(
//										"matching positions passed before: "
//												+ matchingPositionsPassed);
								matchingPositionsPassed += (Integer
										.parseInt(exonEnds[i])
										- Integer.parseInt(exonStarts[i]) + 1);
								newGenomicCigar += (Integer
										.parseInt(exonEnds[i])
										- Integer.parseInt(exonStarts[i]) + 1)
										+ "M";

//								MappingLogger.getLogger().debug(
//										"matching positions passed after: "
//												+ matchingPositionsPassed);
								// newGenomicCigar += readHit
								// .getCigarString()
								// .substring(
								// convertedBases,
								// convertedBases
								// + Integer
								// .parseInt(exonEnds[i])
								// - Integer
								// .parseInt(exonStarts[i])
								// + 1);
								// convertedBases +=
								// Integer.parseInt(exonEnds[i])
								// - Integer.parseInt(exonStarts[i]) + 1;
							}
							if (i < exonStarts.length - 1) {
								int intronLength = (Integer
										.parseInt(exonStarts[i + 1])
										- Integer.parseInt(exonEnds[i]) - 1);
								if (intronLength <= 0) {
									break;
								}
								newGenomicCigar += intronLength + "N";
								spliceJunctionCovered = true;
							} else {
								break;
							}
						}
					}
				} else if (transcriptStrand.equals("-1")) {
					// if (readHit.getReadName().equals(
					// "SEQ_ID:>ENST00000488297|5|87677774|87678542-3:9")) {
					// MappingLogger.getLogger().debug(
					// "read: " + readHit.getReadString());
					// MappingLogger.getLogger().debug(
					// "legnth: " + readHit.getReadLength() + "; start: "
					// + readHit.getAlignmentStart() + " to end "
					// + readHit.getAlignmentEnd());
					// }
					int newGenomicPositionEnd = -1;
					for (int i = exonStarts.length - 1; i >= 0; i--) {
						int lengthPassedTemp = lengthPassed;
						lengthPassed += Integer.parseInt(exonEnds[i])
								- Integer.parseInt(exonStarts[i]) + 1;
						// if (readHit.getReadName().equals(
						// "SEQ_ID:>ENST00000488297|5|87677774|87678542-3:9")) {
						// MappingLogger.getLogger().debug(
						// "lp: " + lengthPassed + "; lpT = "
						// + lengthPassedTemp);
						// }

						if (readHit.getAlignmentStart() <= lengthPassed) {
							if (newGenomicPositionEnd == -1) {
								newGenomicPositionEnd = Integer
										.parseInt(exonEnds[i])
										- (readHit.getAlignmentStart() - lengthPassedTemp)
										+ 1;
								// if (readHit
								// .getReadName()
								// .equals("SEQ_ID:>ENST00000488297|5|87677774|87678542-3:9"))
								// {
								// MappingLogger.getLogger().debug(
								// "newEndPos: " + newGenomicPositionEnd);
								// }

							}
						}
						if (readHit.getAlignmentEnd() <= (lengthPassed)) {
							if (newGenomicPositionEnd <= Integer
									.parseInt(exonEnds[i])) {
								// read maps to 1 exon continuously --- MAYBE
								// EXCHANGE THIS WITH THE CIGAR OF THE MAPPING
								// ITSELF TO AUTOMATICALLY INCLUDE
								// INDELS?!?!?!?!?!?!
								// newGenomicCigar = (readHit.getReadLength()) +
								// "M";
								newGenomicCigar = readHit.getCigarString();
								newGenomicPositionStart = newGenomicPositionEnd
										- readHit.getReadLength() + 1;
							} else {
								if (readHit.getCigarString().contains(
										CigarOperator.DELETION.toString())
										|| readHit.getCigarString().contains(
												CigarOperator.INSERTION
														.toString())) {
									missedTranscriptAlignments++;
									break;
								}
								// fill in Ms from exon-intron boundary to read
								// end
								newGenomicCigar = (readHit.getAlignmentEnd() - lengthPassedTemp)
										+ "M" + newGenomicCigar;
								newGenomicPositionStart = Integer
										.parseInt(exonEnds[i])
										- (readHit.getAlignmentEnd() - lengthPassedTemp)
										+ 1;
							}
							// if (readHit
							// .getReadName()
							// .equals("SEQ_ID:>ENST00000488297|5|87677774|87678542-3:9"))
							// {
							// MappingLogger.getLogger().debug(
							// "newStartPos: " + newGenomicPositionStart);
							// }
							break;
						} else if (newGenomicPositionEnd != -1) {
							if (readHit.getCigarString().contains(
									CigarOperator.DELETION.toString())
									|| readHit.getCigarString().contains(
											CigarOperator.INSERTION.toString())) {
								missedTranscriptAlignments++;
								break;
							}
							// if (true) {
							// continue;
							// }
							if (newGenomicPositionEnd < Integer
									.parseInt(exonEnds[i])) {
								// fill in Ms for alignment start to exon-intron
								// boundary
								newGenomicCigar = (newGenomicPositionEnd
										- Integer.parseInt(exonStarts[i]) + 1)
										+ "M" + newGenomicCigar;
							} else {
								// fill in entire exon
								newGenomicCigar = (Integer
										.parseInt(exonEnds[i])
										- Integer.parseInt(exonStarts[i]) + 1)
										+ "M" + newGenomicCigar;
							}
							if (i >= 1) {
								int intronLength = (Integer
										.parseInt(exonStarts[i])
										- Integer.parseInt(exonEnds[i - 1]) - 1);
								if (intronLength <= 0) {
									// SOMETHING IS WRONG HERE: MAYBE THE
									// SORTING
									// ALGORITHM ISNT CORRENT WORKING ON THAT
									// LARGE
									// NUMBERS??? OR WHAT THE HELL???
									// OUTPUT ALSO ALL EXON VALUES FROM
									// ARRAYS!!!!
									// ALSO PAIRS!!!
									// MappingLogger.getLogger().error(
									// read.getReferenceName());
									break;
								}
								newGenomicCigar = intronLength + "N"
										+ newGenomicCigar;
							} else {
								// MappingLogger.getLogger().error(
								// "match exceeds last exon on rev strand");
								break;
							}
						}
					}
				}

				// SEQ_ID:>ENST00000496359|10|81564241|81586109-2:3
				// if (readHit.getReadName().equals(
				// "SEQ_ID:>ENST00000476173|10|81839051|81851108-2:11")
				// || readHit.getReadName().equals(
				// "SEQ_ID:>ENST00000488297|5|87677774|87678542-3:9")) {
				// MappingLogger.getLogger().debug(
				// "found with cigar: " + newGenomicCigar + " and new start "
				// + newGenomicPositionStart + "; ref map "
				// + readHit.getReferenceName());
				// // System.exit(0);
				// }

				// MappingLogger.getLogger().debug("newCigar: " +
				// newGenomicCigar);

				// IF NO GENOMIC POSITION COULD BE LOCALIZED!!! CONTINUE!!!
				if (newGenomicPositionStart == (-1)) {
					continue;
				}
				if (newGenomicCigar.contains("N")) {
					splicedReads++;
				}

				// HASHMAP RESULTS IN "OUT OF MEMORY" ERROR FOR LARGE
				// DATASETS!!!s
				if (!readIDToRead.containsKey(readHit.getReadName())) {
					Read newReadObj = new Read(readHit.getReadName());
					if (!readHit.getNotPrimaryAlignmentFlag()) {
						newReadObj.setPrimaryIndex(newReadObj.getGenesHitted()
								.size());
					}
					newReadObj.addGeneHit(readNameSplitted[0]);
					newReadObj.addHitToGenomicPosition(newGenomicPositionStart);
					newReadObj.addRecord(readHit);
					newReadObj.addHitToGenomicCigars(newGenomicCigar);

					readIDToRead.put(readHit.getReadName(), newReadObj);

				} else {
					Read oldReadObj = readIDToRead.get(readHit.getReadName());
					if (!readHit.getNotPrimaryAlignmentFlag()) {
						oldReadObj.setPrimaryIndex(oldReadObj.getGenesHitted()
								.size());
					}
					oldReadObj.addGeneHit(readNameSplitted[0]);
					oldReadObj.addHitToGenomicPosition(newGenomicPositionStart);
					oldReadObj.addRecord(readHit);
					oldReadObj.addHitToGenomicCigars(newGenomicCigar);
				}
				// END OF FOR OVER ALL READ-HITS
			}
			// save hit for the last read
			printReadsToBamFile(readIDToRead);

			MappingLogger.getLogger().debug(
					"an overall of " + mappedReads
							+ " reads have high mapping potential.");
			MappingLogger
					.getLogger()
					.info("an overall of " + splicedReads
							+ " reads are spanning at least 1 splice junction!");
			MappingLogger
					.getLogger()
					.info(missedTranscriptAlignments
							+ " transcript alignments were skipped due to indels + splicing");

			transcriptSamFileReader.close();
		} catch (FileNotFoundException e) {
			MappingLogger.getLogger()
					.error("File not found: " + e.getMessage());
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		MappingLogger.getLogger().debug("Saving BAM file...");
		combinedResults.close();
	}

	private void printReadsToBamFile(HashMap<String, Read> readIDToRead) {
		for (String readID : readIDToRead.keySet()) {
			Read currentRead = readIDToRead.get(readID);
			boolean add = true;

			int newGenomicLocationTemp = currentRead.getNewGenomicPositions()
					.get(0);
			for (int newGenomicLocation : currentRead.getNewGenomicPositions()) {
				if (newGenomicLocation != newGenomicLocationTemp) {
					add = false;
					break;
				}
			}

			// add hit to genomic mapping BAM-file
			if (add) {
				SAMRecord newRecord = currentRead.getRecords().get(
						currentRead.getPrimaryIndex());
				// MappingLogger.getLogger().debug("add readhit");

				// reset some values
				String[] readRecordNameSplitted = newRecord.getReferenceName()
						.split("\\|");
				if (genomeSamFileHeader.getSequenceIndex("chr"
						+ readRecordNameSplitted[2]) == -1) {
					continue;
				}
				newRecord.setHeader(genomeSamFileHeader);
				if (readRecordNameSplitted[2].equals("MT")) {
					readRecordNameSplitted[2] = "M";
				}

				newRecord.setReferenceIndex(genomeSamFileHeader
						.getSequenceIndex("chr" + readRecordNameSplitted[2]));
				newRecord.setAlignmentStart(currentRead
						.getNewGenomicPositions().get(
								currentRead.getPrimaryIndex()));

				newRecord.setCigarString(currentRead.getNewGenomicCigars().get(
						currentRead.getPrimaryIndex()));
				newRecord.setMappingQuality(10);
				// newRecord.set

				if (readRecordNameSplitted[5].equals("-1")) {
					// rev compl
					if (newRecord.getReadNegativeStrandFlag()) {
						newRecord.setFlags(newRecord.getFlags() - 16);
					} else {
						newRecord.setFlags(newRecord.getFlags() + 16);
					}

					byte[] revComplSeq = newRecord.getReadBases();
					htsjdk.samtools.util.SequenceUtil
							.reverseComplement(revComplSeq);
					newRecord.setReadBases(revComplSeq);
				}
				// for (CigarElement elem : newRecord.getCigar()
				// .getCigarElements()) {
				// MappingLogger.getLogger().debug(
				// "operator=" + elem.getOperator().toString());
				// }

				// add modified record!!!
				combinedResults.addAlignment(newRecord);

				mappedReads++;
			}
		}
	}
}
