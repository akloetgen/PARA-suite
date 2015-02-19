package utils.pileupclusters;

import htsjdk.samtools.AlignmentBlock;
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
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFrame;

import main.MappingLogger;

import org.math.plot.Plot2DPanel;

/**
 * 
 * @author akloetgen
 * 
 */
public class PileupClusters {

	private IndexedFastaSequenceFile reference;
	private boolean isDerivingErrorProfile = false;
	private int skippedDueIndel = 0;
	private SNPCalling snpCalling;
	private int snpHit;

	public void calculateReadPileups(String alignmentFile,
			String referenceFile, String outputFile, String snpVcfFile,
			int minReadCoverage) {
		try {
			reference = new IndexedFastaSequenceFile(new File(referenceFile));
			File mappingFile = new File(alignmentFile);
			final SamReaderFactory factory = SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
					.validationStringency(ValidationStringency.LENIENT);
			SamReader samFileReader = factory.open(mappingFile);

			FileWriter fileWriterPileup = new FileWriter(outputFile);
			FileWriter fileWriterErrorProfile = new FileWriter(mappingFile
					+ ".sitefrequency.tsv");
			FileWriter fileWriterCCRsFASTA = new FileWriter(outputFile
					+ ".ccr.fasta");
			FileWriter fileWriterCCRsInfo = new FileWriter(outputFile
					+ ".ccr.tsv");
			FileWriter fileWriterAllelePositions = new FileWriter(mappingFile
					+ ".sitepositions.tsv");
			snpCalling = new SNPCalling(snpVcfFile);

			if (!samFileReader.getFileHeader().getSortOrder()
					.equals(SAMFileHeader.SortOrder.coordinate)) {
				MappingLogger
						.getLogger()
						.error("BAM file is not sorted. Please provide"
								+ " a sorted BAM-file as input alignment file.");
				System.exit(0);
			}

			// prepare output-file with header
			String lineOutput = "ClusterID\tChr\tStart\tEnd\tStrand\t#reads\t#T2C\t#T2C sites\tT2C Fraction\tSeqenece\tCombStrand\tSeqLength"
					+ System.getProperty("line.separator");
			fileWriterPileup.write(lineOutput);

			String lineOutputCCR = "Protein_Group\tCluster ID\tStrand\tChromosome"
					+ "\tCluster_Begin\tCluster_End\tAnchor_FlankSeq_Begin\tAnchor_FlankSeq_End"
					+ "\tAnchor_FlankSeq\tAnchor_Position\tCluster_Clone_Count\tNumber_of_T2C_Positions"
					+ "\tT2C_Freq_at_Anchor_Position\tT2C_Fract_at_Anchor_Position\tT2C_Freq_Whole_Cluster"
					+ "\tT2C_Fract_Whole_Cluster"
					+ System.getProperty("line.separator");
			fileWriterCCRsInfo.write(lineOutputCCR);

			// general information
			int numReadsProcessed = 0;
			int numCrosslinkedClusters = 0;
			int numAllelePositions = 0;

			// allele frequency information
			List<Double> alleleFrequencyInformation = new LinkedList<Double>();
			int[] allelePositions = new int[51];
			boolean[] alleleFrequencyPositionsTemp = new boolean[51];

			// cluster information
			int tempClusterStart = 0;
			int tempClusterEnd = 0;
			String tempClusterChr = "";
			byte[] tempClusterBytes = new byte[0];
			int numReadsPerCluster = 0;
			int numT2CMutationPerCluster = 0;
			int numT2CSitesPerCluster = 0;
			int doubleStranded = 0;
			HashMap<Integer, Integer> mutationMap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> baseCoveredMap = new HashMap<Integer, Integer>();
			StrandOrientation isReverse = new StrandOrientation();
			isReverse.setReverse(false);
			boolean tempIsReverse = false;
			List<Double> sortedT2CAmounts = new LinkedList<Double>();
			String clusterID = "";
			int runningID = 1;

			int SNPs = 0;

			for (SAMRecord readHit : samFileReader) {
				numReadsProcessed++;
				if (numReadsProcessed % 250000 == 0) {
					MappingLogger.getLogger().debug(
							numReadsProcessed + " number reads processed.");
					// break;
				}

				// MORE CHECKS????????????????
				if (readHit.getReadUnmappedFlag()) {
					continue;
				}

				// INSERTION = READ IS LONGER
				// DELETION = REFERENCE IS LONGER
				if (readHit.getCigarString().contains("I")
						|| readHit.getCigarString().contains("D")) {
					// TODO so far indel reads are not handled, so continue!
					// continue;
					// skippedDueIndel++;
					// MappingLogger.getLogger().debug(
					// "CIGAR: " + readHit.getCigarString());
					// byte readSequence[] = readHit.getReadBases();
					// byte refSequenceForRead[] = reference.getSubsequenceAt(
					// readHit.getReferenceName(),
					// readHit.getAlignmentStart(),
					// readHit.getAlignmentEnd()).getBases();
					// MappingLogger.getLogger().debug(
					// "Read length = " + readSequence.length);
					// MappingLogger.getLogger().debug(
					// "Ref length = " + refSequenceForRead.length);
					// MappingLogger.getLogger().debug(
					// "number align-blocks = "
					// + readHit.getAlignmentBlocks().size());
				}

				// CHECK FOR TOO MANY ERRORS!! THIS CAN BE DONE SMARTER IN THE
				// LOOP AT THE BOTTOM!!! this is not necessary for correct
				// alignment, isnt it?
				/*
				 * int overallMutations = 0; byte readSequence[] =
				 * readHit.getReadBases(); byte refSequenceForRead[] =
				 * reference.getSubsequenceAt( readHit.getReferenceName(),
				 * readHit.getAlignmentStart(), readHit.getAlignmentEnd())
				 * .getBases(); if (readHit.getReadNegativeStrandFlag()) {
				 * net.sf.samtools.util.SequenceUtil
				 * .reverseComplement(readSequence);
				 * net.sf.samtools.util.SequenceUtil
				 * .reverseComplement(refSequenceForRead); } for (int i = 0; i <
				 * readSequence.length; i++) { if
				 * (calculateArrayPos(readSequence[i]) !=
				 * calculateArrayPos(refSequenceForRead[i])) {
				 * overallMutations++; } } if (overallMutations >
				 * overallMutationsThreshold) {
				 * 
				 * StringBuffer tempReadSequence = new StringBuffer(); for (int
				 * o = 0; o < readSequence.length; o++) {
				 * tempReadSequence.append((char) readSequence[o]); }
				 * StringBuffer tempRefSequence = new StringBuffer(); for (int o
				 * = 0; o < refSequenceForRead.length; o++) {
				 * tempRefSequence.append((char) refSequenceForRead[o]); }
				 * MappingLogger.getLogger().debug( "read id: " +
				 * readHit.getReadName() + "\nread seq: " + tempReadSequence +
				 * "\nref seq: " + tempRefSequence);
				 * 
				 * continue; }
				 */

				if ((tempClusterEnd - readHit.getAlignmentStart()) < 5
						|| !readHit.getReferenceName().equals(tempClusterChr)) {
					// NEW CLUSTER FOUND!
					boolean isSave = true;
					// MappingLogger.getLogger().debug(
					// "tempClusterEnd="
					// + tempClusterEnd
					// + "; readAlignmentStart="
					// + readHit.getAlignmentStart()
					// + "; result="
					// + (tempClusterEnd - readHit
					// .getAlignmentStart()));
					double fractionT2CMutationPerCluster = 0.0;
					if (numReadsPerCluster >= minReadCoverage) {
						numT2CSitesPerCluster = mutationMap.size();
						int tempBestMutationPos = -1;
						double tempBestMutationValue = 0.0;

						// EXCLUDE SNPs FROM FOUND T2C MUTATIONS WITHIN A
						// CLUSTER
						HashMap<Integer, Integer> mutationMapTemp = new HashMap<Integer, Integer>();
						mutationMapTemp.putAll(mutationMap);
						for (int mutationKey : mutationMap.keySet()) {
							// MappingLogger.getLogger().debug(
							// "test for SNP at " + tempClusterChr + ":"
							// + mutationKey);
							if (snpCalling.querySNP(tempClusterChr,
									mutationKey, "T", "C")) {
								mutationMapTemp.remove(mutationKey);
								snpHit++;
							}
							// System.out.println(mutationKey);
						}
						mutationMap.clear();
						mutationMap.putAll(mutationMapTemp);

						if (numT2CSitesPerCluster > 0) {
							sortedT2CAmounts.clear();

							for (int mutationKey : mutationMap.keySet()) {
								// TEST FOR SNP?!?!?!?!?!
								if (((double) mutationMap.get(mutationKey) / baseCoveredMap
										.get(mutationKey)) == 1.0) {
									// MappingLogger
									// .getLogger()
									// .debug("found 100% T-C, skip as this seems to be a SNP");
									// AT THIS POINT, ALL INFOs HAVE TO BE
									// DOWNCALCULATED!!!! SO FAR, SKIP THIS
									// UNTIL REAL SNP FINDING IS IMPLEMENTED

									// SNPs++;
									// continue;
								}
								// GO ON AS USUAL
								double mutationValue = (double) mutationMap
										.get(mutationKey)
										/ baseCoveredMap.get(mutationKey);
								// MappingLogger.getLogger().debug(
								// "mutVal=" + mutationValue
								// + "; tempBestVal="
								// + tempBestMutationValue);
								if (mutationValue >= tempBestMutationValue) {
									// MappingLogger.getLogger().debug(
									// "CHANGE: mutKey=" + mutationKey);
									tempBestMutationValue = mutationValue;
									tempBestMutationPos = mutationKey;
									// MappingLogger.getLogger().debug(
									// "CHANGE: tempBestMutationPos="
									// + tempBestMutationPos);
								}
								sortedT2CAmounts.add(mutationValue);
							}
							Collections.sort(sortedT2CAmounts,
									Collections.reverseOrder());

							// calculate frequency only if it is sure that this
							// is a T2C allele and not a sequencing error
							// allele!! So check for high T2C, e.g. above 20%.
							if (sortedT2CAmounts.size() == 0) {
								// if (tempClusterStart == 20891481) {
								// MappingLogger.getLogger().error(
								// "Found 20891481. DO NOT SAVE");
								// }
								// isSave = false;
							} else if (sumUpList(sortedT2CAmounts) >= 0.2) {
								for (int k = 0; k < sortedT2CAmounts.size(); k++) {
									if (alleleFrequencyInformation.size() > k) {
										alleleFrequencyInformation
												.set(k,
														(double) (alleleFrequencyInformation
																.get(k) + sortedT2CAmounts
																.get(k)));
									} else if (alleleFrequencyInformation
											.isEmpty()) {
										alleleFrequencyInformation
												.addAll(sortedT2CAmounts);
									} else {
										alleleFrequencyInformation
												.add(sortedT2CAmounts.get(k));
									}
								}
								numCrosslinkedClusters++;
								for (int j = 0; j < alleleFrequencyPositionsTemp.length; j++) {
									if (alleleFrequencyPositionsTemp[j]) {
										allelePositions[j]++;
										numAllelePositions++;
									}
								}
							}

							for (Double t2cAmount : sortedT2CAmounts) {
								// if (tempClusterStart == 17469751) {
								// MappingLogger.getLogger().debug(
								// "17469751: t2c-amount: "
								// + t2cAmount);
								//
								// }
								fractionT2CMutationPerCluster += t2cAmount;
							}
							// PRINT OUT CCR FILE FOR OTHER PURPOSES.
							if (tempBestMutationPos > 0) {
								// IS NOT ADJUSTED FOR SPLICED READS!!!
								byte[] ccrSequenceBytes = reference
										.getSubsequenceAt(tempClusterChr,
												tempBestMutationPos - 20,
												tempBestMutationPos + 20)
										.getBases();
								if (isReverse.getStrandOrientation()
										.equals("-")) {
									htsjdk.samtools.util.SequenceUtil
											.reverseComplement(ccrSequenceBytes);
								}
								StringBuffer ccrSequence = new StringBuffer();
								for (int basePosition = 0; basePosition < ccrSequenceBytes.length; basePosition++) {
									ccrSequence
											.append(Character
													.toUpperCase((char) ccrSequenceBytes[basePosition]));
								}

								fileWriterCCRsFASTA.write(">" + clusterID
										+ " 20-anchor-20 " + tempClusterChr
										+ ":"
										+ isReverse.getStrandOrientation()
										+ ":" + (tempBestMutationPos - 20)
										+ "-" + (tempBestMutationPos + 20)
										+ System.getProperty("line.separator"));
								fileWriterCCRsFASTA.write(ccrSequence
										+ System.getProperty("line.separator"));

								lineOutputCCR = "Gene\t" + clusterID + "\t"
										+ isReverse.getStrandOrientation()
										+ "\t" + tempClusterChr + "\t"
										+ tempClusterStart + "\t"
										+ tempClusterEnd + "\t"
										+ (tempBestMutationPos - 20) + "\t"
										+ (tempBestMutationPos + 20) + "\t"
										+ ccrSequence + "\t"
										+ tempBestMutationPos + "\t"
										+ numReadsPerCluster + "\t"
										+ numT2CSitesPerCluster + "\t"
										+ mutationMap.get(tempBestMutationPos)
										+ "\t" + tempBestMutationValue + "\t"
										+ numT2CMutationPerCluster + "\t"
										+ fractionT2CMutationPerCluster;
								fileWriterCCRsInfo.write(lineOutputCCR
										+ System.getProperty("line.separator"));

							}
						}
						if (isSave) {

							// tempClusterBytes = reference.getSubsequenceAt(
							// tempClusterChr, tempClusterStart,
							// tempClusterEnd).getBases();
							if (tempIsReverse) {
								htsjdk.samtools.util.SequenceUtil
										.reverseComplement(tempClusterBytes);
							}
							StringBuffer tempClusterSequence = new StringBuffer();
							for (int i = 0; i < tempClusterBytes.length; i++) {
								tempClusterSequence
										.append((char) tempClusterBytes[i]);
							}

							// ADD OUTPUT FOR SPLICED CLUSTER!!!

							lineOutput = clusterID + "\t" + tempClusterChr
									+ "\t" + tempClusterStart + "\t"
									+ tempClusterEnd + "\t"
									+ orientationToString(tempIsReverse) + "\t"
									+ numReadsPerCluster + "\t"
									+ numT2CMutationPerCluster + "\t"
									+ numT2CSitesPerCluster + "\t"
									+ fractionT2CMutationPerCluster + "\t"
									+ tempClusterSequence + "\t"
									+ isReverse.getStrandOrientation() + "\t"
									+ tempClusterSequence.length()
									+ System.getProperty("line.separator");
							fileWriterPileup.write(lineOutput);
						}
					}

					tempClusterStart = readHit.getAlignmentStart();
					tempClusterEnd = readHit.getAlignmentEnd();
					tempClusterChr = readHit.getReferenceName();
					numReadsPerCluster = 1;
					numT2CMutationPerCluster = 0;
					numT2CSitesPerCluster = 0;
					isReverse.setReverse(false);
					mutationMap.clear();
					baseCoveredMap.clear();
					runningID++;
					clusterID = "cl_" + runningID + "_" + tempClusterChr;
					alleleFrequencyPositionsTemp = new boolean[51];
					// numA2GMutationPerCluster = 0;

					// calculate everything
					numT2CMutationPerCluster = calculateClusterInformation(
							readHit, isReverse, mutationMap, baseCoveredMap,
							alleleFrequencyPositionsTemp,
							numT2CMutationPerCluster);
					tempIsReverse = isReverse.getReverse();

					// INITIAL CLUSTER SEQUENCE SET
					tempClusterBytes = new byte[0];
					for (AlignmentBlock block : readHit.getAlignmentBlocks()) {
						byte[] tempReferenceSequence = reference
								.getSubsequenceAt(
										readHit.getReferenceName(),
										block.getReferenceStart(),
										block.getReferenceStart()
												+ block.getLength() - 1)
								.getBases();
						tempClusterBytes = mergeByteArrays(tempClusterBytes,
								tempReferenceSequence);
					}

				} else {
					// MEMBER FOR PREVIOUS CLUSTER FOUND!
					// calculate everything
					tempClusterChr = readHit.getReferenceName();
					// if (readHit.getAlignmentStart() < tempClusterStart) {

					// byte[] additionalNucs = reference.getSubsequenceAt(
					// readHit.getReferenceName(),
					// readHit.getAlignmentStart(), tempClusterStart)
					// .getBases();
					// tempClusterStart = readHit.getAlignmentStart();

					// UPDATE CLUSTER SEQ

					// }
					if (readHit.getAlignmentEnd() > tempClusterEnd) {
						// UPDATE CLUSTER SEQ
						// if (tempClusterStart == 71146850) {
						// MappingLogger.getLogger().debug(
						// "Found 71146850. Process additional read");
						// }
						for (AlignmentBlock block : readHit
								.getAlignmentBlocks()) {
							// MappingLogger.getLogger().debug(
							// "block start in ref: "
							// + block.getReferenceStart());
							if ((block.getReferenceStart() + block.getLength() - 1) < tempClusterEnd) {
								continue;
							}

							byte[] additionalNucs = reference.getSubsequenceAt(
									readHit.getReferenceName(),
									block.getReferenceStart(),
									block.getReferenceStart()
											+ block.getLength() - 1).getBases();
							int overhangPosition = tempClusterEnd
									- block.getReferenceStart() - 1;
							// if (tempClusterStart == 71146850) {
							// MappingLogger
							// .getLogger()
							// .debug("Found 71146850. adds nucs to the end with overhangPos = "
							// + overhangPosition
							// + " and length "
							// + (block.getReferenceStart()
							// + block.getLength() - 1
							// - tempClusterEnd - overhangPosition));
							// }
							if (overhangPosition > 0) {
								// if (tempClusterStart == 71146850) {
								// MappingLogger.getLogger().debug(
								// "Found 71146850. merge with 0 - "
								// + tempClusterBytes.length
								// + " and "
								// + overhangPosition + "-"
								// + block.getLength());
								// }
								tempClusterBytes = mergeByteSubArrays(
										tempClusterBytes, 0,
										tempClusterBytes.length,
										additionalNucs, overhangPosition,
										block.getLength());
							} else {
								tempClusterBytes = mergeByteArrays(
										additionalNucs, tempClusterBytes);
							}

						}

						tempClusterEnd = readHit.getAlignmentEnd();
					}
					numReadsPerCluster++;

					numT2CMutationPerCluster = calculateClusterInformation(
							readHit, isReverse, mutationMap, baseCoveredMap,
							alleleFrequencyPositionsTemp,
							numT2CMutationPerCluster);
					// MappingLogger.getLogger().debug("temP: " +
					// tempIsReverse);
					// if (isReverse.getReverse() == null) {
					// MappingLogger.getLogger().debug("isRev = null");
					// } else {
					// MappingLogger.getLogger().debug(
					// "isRev = " + isReverse.getReverse());
					//
					// }
					if (isReverse.getReverse() != null
							&& tempIsReverse != isReverse.getReverse()) {
						doubleStranded++;
						// MappingLogger.getLogger().debug(
						// "double stranded found.");
						isReverse.setReverse(null);
					}
				}
			}

			MappingLogger.getLogger().debug(
					"Double stranded clusters found: " + doubleStranded);
			MappingLogger.getLogger()
					.debug("Loci found that are SNPs: " + SNPs);
			MappingLogger.getLogger().debug(
					skippedDueIndel + "insertion or deletion skipped");
			MappingLogger.getLogger().debug(
					"T-C mutations identified as SNPs: " + snpHit);

			// LETZTER CLUSTER WIRD NOCH NICHT EINGETRAGEN!!! NACHHOLEN!!! EVTL
			// WRITING METHODE AUSLAGERN UND 2 MAL AUFRUFEN?!??!!

			double[] alleleFrequencyInformationArray = new double[alleleFrequencyInformation
					.size()];
			for (int i = 0; i < alleleFrequencyInformationArray.length; i++) {
				alleleFrequencyInformationArray[i] = (double) alleleFrequencyInformation
						.get(i).doubleValue() / numCrosslinkedClusters;
				fileWriterErrorProfile.write(alleleFrequencyInformationArray[i]
						+ System.getProperty("line.separator"));
			}
			double[] alleleFrequencyPositions = new double[allelePositions.length];
			for (int i = 0; i < allelePositions.length; i++) {
				alleleFrequencyPositions[i] = (double) allelePositions[i]
						/ numAllelePositions;
				fileWriterAllelePositions.write(alleleFrequencyPositions[i]
						+ System.getProperty("line.separator"));
			}

			if (isDerivingErrorProfile) {
				Plot2DPanel plot = new Plot2DPanel();

				plot.addBarPlot("T2C allele frequencies", Color.BLUE,
						alleleFrequencyInformationArray);

				plot.getAxis(0).setLabelText(
						"sorted T2C alleles within clusters");
				plot.getAxis(1).setLabelText("% of T2C rate for given allele");

				plot.addLegend("WEST");

				// plot.setFixedBounds(1, 0, 1);
				plot.setFixedBounds(0, 0, 35);

				JFrame frame = new JFrame(
						"Allele frequencies for PAR-CLIP data");
				frame.setSize(750, 400);
				frame.setLocation(100, 100);
				frame.setContentPane(plot);
				frame.setVisible(true);
			}
			samFileReader.close();
			reference.close();
			fileWriterPileup.close();
			fileWriterErrorProfile.close();
			fileWriterCCRsFASTA.close();
			fileWriterCCRsInfo.close();
			fileWriterAllelePositions.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private int calculateClusterInformation(SAMRecord readHit,
			StrandOrientation isReverse, HashMap<Integer, Integer> mutationMap,
			HashMap<Integer, Integer> baseCoveredMap,
			boolean[] mutationMapInRead, int numT2CMutationPerCluster) {
		byte tempReadSequence[] = readHit.getReadBases();
		byte readSequence[] = new byte[0];
		byte refSequenceForRead[] = new byte[0];
		// for (CigarElement elem : readHit.getCigar().getCigarElements()) {
		// if (elem.equals(CigarOperator.M.toString())) {
		//
		// }
		// }
		for (AlignmentBlock block : readHit.getAlignmentBlocks()) {
			// int currentReadSequenceLength = 0;
			// if (readSequence.length > 0) {
			// currentReadSequenceLength = readSequence.length + 1;
			// }
			readSequence = mergeByteSubArrays(readSequence, 0,
					readSequence.length, tempReadSequence,
					block.getReadStart() - 1,
					block.getReadStart() - 1 + block.getLength());
			// if (readHit.getCigarString().contains("D")) {
			//
			// MappingLogger.getLogger().debug(
			// "Block-start: " + block.getReadStart()
			// + "; block-end: " + blockEnd);
			// }
			byte[] tempReferenceSequence = reference.getSubsequenceAt(
					readHit.getReferenceName(), block.getReferenceStart(),
					block.getReferenceStart() + block.getLength() - 1)
					.getBases();
			refSequenceForRead = mergeByteArrays(refSequenceForRead,
					tempReferenceSequence);
		}
		// byte refSequenceForRead[] = reference.getSubsequenceAt(
		// readHit.getReferenceName(), readHit.getAlignmentStart(),
		// readHit.getAlignmentEnd()).getBases();
		// int startPosition;
		if (readHit.getReadNegativeStrandFlag()) {
			htsjdk.samtools.util.SequenceUtil.reverseComplement(readSequence);
			htsjdk.samtools.util.SequenceUtil
					.reverseComplement(refSequenceForRead);
			if (isReverse != null) {
				isReverse.setReverse(true);
			}

			// StringBuffer tempReadSequence = new StringBuffer();
			// for (int o = 0; o < readSequence.length; o++) {
			// tempReadSequence.append((char) readSequence[o]);
			// }
			// MappingLogger.getLogger().debug(
			// "read id: " + read.getReadName() + "read seq: "
			// + tempReadSequence);

			// startPosition = read.getAlignmentEnd();
			// } else {
			// startPosition = read.getAlignmentStart();
		}

		// if (readSequence.length != refSequenceForRead.length) {
		// if (readHit.getCigarString().contains("I")) {
		// // TODO error handling needed: this routine is calculating
		// // information for a read, but what if it has an indel??!?! But it
		// // is already included in the cluster?!?!?!?
		// // MappingLogger.getLogger().debug(
		// // "reflength=" + refSequenceForRead.length
		// // + " and readlenght=" + readSequence.length);
		// MappingLogger.getLogger().debug(
		// "Cigar: " + readHit.getCigarString());
		// StringBuffer tempRefSequence = new StringBuffer();
		// StringBuffer tempReadSequenceBuffer = new StringBuffer();
		// for (int i = 0; i < refSequenceForRead.length; i++) {
		// tempRefSequence.append((char) refSequenceForRead[i]);
		// tempReadSequenceBuffer.append((char) readSequence[i]);
		// }
		//
		// MappingLogger.getLogger()
		// .debug("ref sequence:\t" + tempRefSequence);
		// MappingLogger.getLogger().debug(
		// "readsequence:\t" + tempReadSequenceBuffer);
		//
		// return 0;
		// }

		for (int i = 0; i < readSequence.length; i++) {
			int checkPosition;
			if (readHit.getReadNegativeStrandFlag()) {
				checkPosition = readHit.getAlignmentEnd() - i;
			} else {
				checkPosition = readHit.getAlignmentStart() + i;
			}
			// CHECK FOR ALL MUTATIONS! IN CASE THAT THERE IS A SHIFTED
			// ALIGNMENT VS.
			// TRANSCRIPTOME, THIS WILL PRODUCE AN ERROR! DO NOT OUTPUT THOSE
			// ALIGNMENTS/CLUSTERS!
			// if (readSequence[i] != refSequenceForRead[i]) {
			// overallMutations++;
			// }

			if (calculateArrayPos(refSequenceForRead[i]) == 3
					&& calculateArrayPos(readSequence[i]) == 1) {
				// if (!snpCalling.querySNP(readHit.getReferenceName(),
				// checkPosition, "T", "C")) {
				// snpHit++;
				// } else {
				numT2CMutationPerCluster++;
				mutationMapInRead[i] = true;
				if (mutationMap.containsKey(checkPosition)) {
					mutationMap.put(checkPosition,
							mutationMap.get(checkPosition) + 1);
				} else {
					mutationMap.put(checkPosition, 1);
				}
				// }
			}/*
			 * else if (isAllowForAG && ((refSequenceForRead[i] == 65 &&
			 * readSequence[i] == 71) || (refSequenceForRead[i] == 97 &&
			 * readSequence[i] == 103))) { numA2GMutationPerCluster++; }
			 */
			if (baseCoveredMap.containsKey(checkPosition)) {
				baseCoveredMap.put(checkPosition,
						baseCoveredMap.get(checkPosition) + 1);
			} else {
				baseCoveredMap.put(checkPosition, 1);
			}
		}
		// if (overallMutations > overallMutationsThreshold) {
		// return -1;
		// }
		return numT2CMutationPerCluster;
	}

	private String orientationToString(boolean tempIsReverse) {
		if (tempIsReverse) {
			return "-";
		} else {
			return "+";
		}
	}

	/**
	 * Converts bytes representing A/C/G/T ASCII characters to index values.
	 * 
	 * @param formerByte
	 *            byte representation of A/C/G/T
	 * @return index value for A/C/G/T
	 */
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

	private byte[] mergeByteArrays(byte[] array1, byte[] array2) {
		byte[] tempByteArray = new byte[array1.length + array2.length];
		for (int i = 0; i < array1.length; i++) {
			tempByteArray[i] = array1[i];
		}
		for (int j = 0; j < array2.length; j++) {
			tempByteArray[j + array1.length] = array2[j];
		}

		return tempByteArray;
	}

	private byte[] mergeByteSubArrays(byte[] array1, int array1Start,
			int array1End, byte[] array2, int array2Start, int array2End) {
		byte[] tempByteArray = new byte[(array1End - array1Start)
				+ (array2End - array2Start)];

		for (int i = array1Start; i < array1End; i++) {
			tempByteArray[i] = array1[i];
		}
		for (int j = 0; j < (array2End - array2Start); j++) {
			tempByteArray[j + (array1End - array1Start)] = array2[j
					+ array2Start];
		}

		return tempByteArray;
	}

	private double sumUpList(List<Double> list) {
		double sum = 0.0;
		for (Double val : list) {
			sum += val;
		}
		return sum;
	}

}
