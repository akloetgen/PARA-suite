package utils.pileupclusters;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMException;
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
 * Hierarchical clustering algorithm to pile up aligned reads into clusters
 * representing RBP-bound regions identified by CLIP.
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
	private int highFrequentError;

	/**
	 * Apply hierarchical clustering and exclude T-C SNPs annotated in the SNP
	 * VCF file.
	 * 
	 * @param alignmentFile
	 *            BAM file of the read alignment
	 * @param referenceFile
	 *            indexed reference genome sequence file
	 * @param outputFile
	 *            file-prefix where all outputs are saved
	 * @param snpVcfFile
	 *            SNP db in a VCF format, requires tabix index!
	 * @param minReadCoverage
	 *            minimal coverage of reads per cluster, e.g. 5
	 */
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
			FileWriter fileWriterLog = new FileWriter(outputFile + ".report");
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
					MappingLogger.getLogger().info(
							numReadsProcessed + " number reads processed.");
					// break;
				}

				// MORE CHECKS????????????????
				if (readHit.getReadUnmappedFlag()) {
					continue;
				}

				// INSERTION = READ IS LONGER
				// DELETION = REFERENCE IS LONGER
				if ((readHit.getCigarString().contains("I") || readHit
						.getCigarString().contains("D"))
						&& readHit.getCigarString().contains("N")) {
					// TODO so far indel reads are not handled, so continue!
					skippedDueIndel++;
					continue;

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

				if ((tempClusterEnd - readHit.getAlignmentStart()) < 5
						|| !readHit.getReferenceName().equals(tempClusterChr)) {
					// NEW CLUSTER FOUND!
					boolean isSave = true;
					double fractionT2CMutationPerCluster = 0.0;
					if (numReadsPerCluster >= minReadCoverage) {
						numT2CSitesPerCluster = mutationMap.size();
						int tempBestMutationPos = -1;
						double tempBestMutationValue = 0.0;

						// exclude T-C SNPs and 100% T-C sites (most likely
						// SNVs)
						HashMap<Integer, Integer> mutationMapTemp = new HashMap<Integer, Integer>();
						mutationMapTemp.putAll(mutationMap);
						for (int mutationKey : mutationMap.keySet()) {
							if (snpCalling.querySNP(tempClusterChr,
									mutationKey, "T", "C")) {
								mutationMapTemp.remove(mutationKey);
								snpHit++;
							}
							if (mutationMap.get(mutationKey) == 1.0) {
								highFrequentError++;
								continue;
							}
						}
						mutationMap.clear();
						mutationMap.putAll(mutationMapTemp);

						if (numT2CSitesPerCluster > 0) {
							sortedT2CAmounts.clear();

							for (int mutationKey : mutationMap.keySet()) {
								// This could be another check for unusual T-C
								// conversion sites. So far excluded.
								// if (((double) mutationMap.get(mutationKey) /
								// baseCoveredMap
								// .get(mutationKey)) == 1.0) {
								// }
								double mutationValue = (double) mutationMap
										.get(mutationKey)
										/ baseCoveredMap.get(mutationKey);

								if (mutationValue >= tempBestMutationValue) {
									tempBestMutationValue = mutationValue;
									tempBestMutationPos = mutationKey;
								}
								sortedT2CAmounts.add(mutationValue);
							}
							Collections.sort(sortedT2CAmounts,
									Collections.reverseOrder());

							// calculate frequency only if it is sure that this
							// is a T2C allele and not a sequencing error
							// allele!! So check for high T2C, e.g. above 20%.
							if (sortedT2CAmounts.size() == 0) {
								// ???
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
								fractionT2CMutationPerCluster += t2cAmount;
							}
							// PRINT OUT CCR FILE FOR OTHER PURPOSES.
							if (tempBestMutationPos > 0) {
								// IS NOT ADJUSTED FOR SPLICED READS!!!
								byte[] ccrSequenceBytes;
								try {
									ccrSequenceBytes = reference
											.getSubsequenceAt(tempClusterChr,
													tempBestMutationPos - 20,
													tempBestMutationPos + 20)
											.getBases();
									if (isReverse.getStrandOrientation()
											.equals("-")) {
										htsjdk.samtools.util.SequenceUtil
												.reverseComplement(ccrSequenceBytes);
									}
								} catch (SAMException e) {
									ccrSequenceBytes = new byte[0];
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

					// calculate everything
					numT2CMutationPerCluster = calculateClusterInformation(
							readHit, isReverse, mutationMap, baseCoveredMap,
							alleleFrequencyPositionsTemp,
							numT2CMutationPerCluster);
					tempIsReverse = isReverse.getReverse();

					// INITIAL CLUSTER SEQUENCE SET
					tempClusterBytes = new byte[0];
					// readHit.getCigar().getCigarElement(0).getLength()

					int currentBlockStartPosition = readHit.getAlignmentStart();
					// for (AlignmentBlock block : readHit.getAlignmentBlocks())
					// {
					for (CigarElement cigarElem : readHit.getCigar()
							.getCigarElements()) {
						byte[] tempReferenceSequence;
						// if (!readHit.getReadNegativeStrandFlag()) {

						// tempReferenceSequence = reference.getSubsequenceAt(
						// readHit.getReferenceName(),
						// block.getReferenceStart(),
						// block.getReferenceStart() + block.getLength()
						// - 1).getBases();

						if (cigarElem.getOperator().equals(CigarOperator.D)
								|| cigarElem.getOperator().equals(
										CigarOperator.M)) {
							tempReferenceSequence = reference.getSubsequenceAt(
									readHit.getReferenceName(),
									currentBlockStartPosition,
									currentBlockStartPosition
											+ cigarElem.getLength() - 1)
									.getBases();
							tempClusterBytes = mergeByteArrays(
									tempClusterBytes, tempReferenceSequence);

							// currentBlockStartPosition +=
							// cigarElem.getLength();
							// } else if (cigarElem.equals(CigarOperator.I)
							// || cigarElem.equals(CigarOperator.N)) {
							// // TODO nothing todo here? delete?!
						}
						if (!cigarElem.getOperator().equals(CigarOperator.I)) {
							currentBlockStartPosition += cigarElem.getLength();
						}
						// } else {
						// tempReferenceSequence = reference.getSubsequenceAt(
						// readHit.getReferenceName(),
						// block.getReferenceStart() + 1,
						// block.getReferenceStart()
						// + block.getLength()).getBases();
						// }

						// blockNumber++;
					}

				} else {
					// MEMBER FOR PREVIOUS CLUSTER FOUND!
					// calculate everything
					tempClusterChr = readHit.getReferenceName();

					if (readHit.getAlignmentEnd() > tempClusterEnd) {
						// UPDATE CLUSTER SEQ
						// for (AlignmentBlock block : readHit
						// .getAlignmentBlocks()) {
						int currentBlockStartPosition = readHit
								.getAlignmentStart();
						for (CigarElement cigarElem : readHit.getCigar()
								.getCigarElements()) {
							if ((currentBlockStartPosition
									+ cigarElem.getLength() - 1) < tempClusterEnd) {
								if (!cigarElem.getOperator().equals(
										CigarOperator.I)) {
									currentBlockStartPosition += cigarElem
											.getLength();
								}
								continue;
							}
							if (cigarElem.getOperator().equals(CigarOperator.D)
									|| cigarElem.getOperator().equals(
											CigarOperator.M)) {
								byte[] additionalNucs = reference
										.getSubsequenceAt(
												readHit.getReferenceName(),
												currentBlockStartPosition,
												currentBlockStartPosition
														+ cigarElem.getLength()
														- 1).getBases();

								int overhangPosition = tempClusterEnd
										- currentBlockStartPosition + 1;
								if (overhangPosition > 0) {

									tempClusterBytes = mergeByteSubArrays(
											tempClusterBytes, 0,
											tempClusterBytes.length,
											additionalNucs, overhangPosition,
											cigarElem.getLength());

								} else {
									// MappingLogger
									// .getLogger()
									// .debug("overhangpositions < 0 ---> WHY?!?!?!?!?");
									tempClusterBytes = mergeByteArrays(
											additionalNucs, tempClusterBytes);
								}
								// if (clusterID
								// .equals("cl_37_chr6:36673421-36675326")) {
								// StringBuffer tempBuffer = new StringBuffer();
								// for (int i = 0; i < tempClusterBytes.length;
								// i++)
								// {
								// tempBuffer
								// .append((char) tempClusterBytes[i]);
								// }
								// MappingLogger.getLogger().debug(
								// "tempCLsuterBytes after concat: "
								// + tempBuffer);
								// }
							}
							tempClusterEnd = readHit.getAlignmentEnd();
							if (!cigarElem.getOperator()
									.equals(CigarOperator.I)) {
								currentBlockStartPosition += cigarElem
										.getLength();
							}
						}
					}
					numReadsPerCluster++;

					numT2CMutationPerCluster = calculateClusterInformation(
							readHit, isReverse, mutationMap, baseCoveredMap,
							alleleFrequencyPositionsTemp,
							numT2CMutationPerCluster);
					if (isReverse.getReverse() != null
							&& tempIsReverse != isReverse.getReverse()) {
						doubleStranded++;
						isReverse.setReverse(null);
					}
				}
			}

			fileWriterLog.write("Double stranded clusters found: "
					+ doubleStranded + System.getProperty("line.separator"));
			fileWriterLog.write("Loci found that are SNPs: " + SNPs
					+ System.getProperty("line.separator"));
			fileWriterLog.write(skippedDueIndel
					+ " insertion or deletion skipped"
					+ System.getProperty("line.separator"));
			fileWriterLog.write("T-C mutations identified as SNPs: " + snpHit
					+ System.getProperty("line.separator"));
			fileWriterLog
					.write("T-C mutations identified as SNVs (100% T-C in 1 site): "
							+ highFrequentError
							+ System.getProperty("line.separator"));

			MappingLogger.getLogger().debug(
					"Double stranded clusters found: " + doubleStranded);
			MappingLogger.getLogger()
					.debug("Loci found that are SNPs: " + SNPs);
			MappingLogger.getLogger().debug(
					skippedDueIndel + "insertion or deletion skipped");
			MappingLogger.getLogger().debug(
					"T-C mutations identified as SNPs: " + snpHit);
			MappingLogger.getLogger().debug(
					"T-C mutations identified as SNVs (100% T-C in 1 site): "
							+ highFrequentError);

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
			fileWriterLog.close();

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

		for (AlignmentBlock block : readHit.getAlignmentBlocks()) {
			readSequence = mergeByteSubArrays(readSequence, 0,
					readSequence.length, tempReadSequence,
					block.getReadStart() - 1,
					block.getReadStart() - 1 + block.getLength());
			byte[] tempReferenceSequence = reference.getSubsequenceAt(
					readHit.getReferenceName(), block.getReferenceStart(),
					block.getReferenceStart() + block.getLength() - 1)
					.getBases();
			refSequenceForRead = mergeByteArrays(refSequenceForRead,
					tempReferenceSequence);
		}

		if (readHit.getReadNegativeStrandFlag()) {
			htsjdk.samtools.util.SequenceUtil.reverseComplement(readSequence);
			htsjdk.samtools.util.SequenceUtil
					.reverseComplement(refSequenceForRead);
			if (isReverse != null) {
				isReverse.setReverse(true);
			}
		}

		// if (readSequence.length != refSequenceForRead.length) {
		// if (readHit.getCigarString().contains("I")) {
		// // TODO error handling needed: this routine is calculating
		// // information for a read, but what if it has an indel??!?! But it
		// // is already included in the cluster?!?!?!?
		// // MappingLogger.getLogger().debug(
		// // "reflength=" + refSequenceForRead.length
		// // + " and readlenght=" + readSequence.length);

		// MappingLogger.getLogger().debug("Cigar: " +
		// readHit.getCigarString());
		// StringBuffer tempRefSequence = new StringBuffer();
		// StringBuffer tempReadSequenceBuffer = new StringBuffer();
		// for (int i = 0; i < refSequenceForRead.length; i++) {
		// tempRefSequence.append((char) refSequenceForRead[i]);
		// tempReadSequenceBuffer.append((char) readSequence[i]);
		// }
		//
		// MappingLogger.getLogger().debug("ref sequence:\t" + tempRefSequence);
		// MappingLogger.getLogger().debug(
		// "readsequence:\t" + tempReadSequenceBuffer);

		for (int i = 0; i < readSequence.length; i++) {
			int checkPosition;
			if (readHit.getReadNegativeStrandFlag()) {
				checkPosition = readHit.getAlignmentEnd() - i;
			} else {
				checkPosition = readHit.getAlignmentStart() + i;
			}
			// CHECK FOR ALL MUTATIONS! IN CASE THAT THERE IS A SHIFTED
			// ALIGNMENT VS. TRANSCRIPTOME, THIS WILL PRODUCE AN ERROR! DO NOT
			// OUTPUT THOSE ALIGNMENTS/CLUSTERS!
			// if (readSequence[i] != refSequenceForRead[i]) {
			// overallMutations++;
			// }

			if (calculateArrayPos(refSequenceForRead[i]) == 3
					&& calculateArrayPos(readSequence[i]) == 1) {
				numT2CMutationPerCluster++;
				mutationMapInRead[i] = true;
				if (mutationMap.containsKey(checkPosition)) {
					mutationMap.put(checkPosition,
							mutationMap.get(checkPosition) + 1);
				} else {
					mutationMap.put(checkPosition, 1);
				}
			}
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
