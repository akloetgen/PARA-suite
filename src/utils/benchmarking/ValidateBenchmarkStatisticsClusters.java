package utils.benchmarking;

import htsjdk.samtools.SAMFormatException;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.HashMap;
import java.util.LinkedList;

import main.MappingLogger;

/**
 * Evaluates the alignment accuracy of a respective read aligner performed on
 * simulated PAR-CLIP read dataset.
 * 
 * @author akloetgen
 * 
 */
public class ValidateBenchmarkStatisticsClusters {

	/**
	 * Calculates alignment accuracy values, including recall, precision and F1
	 * score for a given alignment of simulated PAR-CLIP reads
	 * 
	 * @param mappingFileName
	 *            filename of the BAM mapping file
	 * @param outStatistics
	 *            filename where the statistics should be written to
	 * @param readsFile
	 *            filename of the raw simulated PAR-CLIP reads file
	 * @param onlyBoundClusters
	 *            whether only bound PAR-CLIP clusters should be evaluated
	 */
	public void calculateBenchmarkStatistics(String clustersFile,
			String outStatistics, String clustersReferenceFile) {

		// final SamReaderFactory factory = SamReaderFactory.makeDefault()
		// .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
		// .validationStringency(ValidationStringency.LENIENT);
		// SamReader samFileReader = factory.open(new File(mappingFileName));
		// File clustersFileOpen = new File(clustersFile);
		// READ CLUSETERS FILE

		File outFile = new File(outStatistics);
		// int matched_correctly = 0;
		// int wrongly_assigned = 0;
		// int clusters_found_overall = 0;
		// int clusters_found_bound = 0;
		int allclusters = 0;
		// int allclusters_bound = 0;

		// equals simulated cluster marked as bound and found as binding-site by
		// the postprocessing algorithm
		int truePositives = 0;
		// equals simulated cluster marked as not bound and found by clustering
		// but didn't passed filtering to be marked as binding-site
		int trueNegatives = 0;
		// equals simulated cluster marked as bound but either not found by
		// alignment or not found by postprocessing algorithm
		int falsePositives = 0;
		// equals simulated cluster marked as not bound and either not found by
		// alignment or not found by postprocessing algorithm
		int falseNegatives = 0;
		int falsePositivesTotal = 0;
		int falseNegativesTotal = 0;

		int processedClusters = 0;

		String clusterStartString = "";
		String clusterEndString = "";
		String clusterBound = "";
		String compl = "";

		int positives = 0;
		int negatives = 0;

		boolean startsWithChr = false;
		// HashMap<String, String> clustersReferenceMap = new
		// HashMap<String,String>();
		HashMap<String, LinkedList<String[]>> clustersReferenceMap = new HashMap<String, LinkedList<String[]>>();

		// either clustering, bmix or paralyzer
		String clusteringTool = "clustering";

		try {

			BufferedReader inClustersReference = new BufferedReader(
					new FileReader(clustersReferenceFile));
			String line = null;
			while ((line = inClustersReference.readLine()) != null) {
				allclusters++;
				// if (line.startsWith("@SEQ_ID")) {
				String[] clusterSplit = line.split("\t");
				if (clustersReferenceMap.get(clusterSplit[1]) == null) {
					// System.out.println("create list for hash-key="
					// + clusterSplit[1]);
					clustersReferenceMap.put(clusterSplit[1],
							new LinkedList<String[]>());
				}
				clustersReferenceMap.get(clusterSplit[1]).add(clusterSplit);
				if (clusterSplit[1].startsWith("chr")) {
					startsWithChr = true;
				}

				if (clusterSplit[4].equals("1")) {
					positives++;
				} else if (clusterSplit[4].equals("0")) {
					negatives++;
				}

				// clusters_found_bound++;
			}
			inClustersReference.close();

			// check whether chromosomes in cluster references start with "chr"

			// if (clustersReferenceMap.get(arg0)
			// .get(0)[1].startsWith("chr")) {
			// startsWithChr = true;
			// }

			BufferedReader inClusters = new BufferedReader(new FileReader(
					clustersFile));
			inClusters.readLine();
			while ((line = inClusters.readLine()) != null) {
				// MappingLogger.getLogger().debug("check line: " + line);
				// for (SAMRecord readHit : samFileReader) {
				if (processedClusters % 5000 == 0) {
					MappingLogger.getLogger().debug(
							"clusters checked: " + processedClusters);
				}
				// fastq-header = information about origin in simulated data???
				// if yes: parse header, compare, output, ready

				// String[] splittedReadName =
				// readHit.getReadName().split("\\|");

				String[] splittedCluster = line.split("\t");
				// compl = readHit.getReadName();
				String clusterChr = null;
				String clusterPassesFilterString = "";
				// clusterStartString = null;
				// clusterEndString = null;

				// CHECK WHICH CLUSTERING PROGRAM WAS PERFORMED TO LOAD THE DATA
				// PROPERLY:
				// if (clusteringTool.equals("clustering")) {

				clusterChr = splittedCluster[1];
				clusterStartString = splittedCluster[2];
				clusterEndString = splittedCluster[3];
				clusterPassesFilterString = splittedCluster[14];

				// } else if (clusteringTool.equals("bmix")) {
				//
				// } else if (clusteringTool.equals("paralyzer")) {

				// }

				int clusterStart = Integer.parseInt(clusterStartString);
				int clusterEnd = Integer.parseInt(clusterEndString);
				int clusterPassesFilter = Integer
						.parseInt(clusterPassesFilterString);

				if (startsWithChr && !clusterChr.startsWith("chr")) {
					clusterChr = "chr" + clusterChr;
				} else if (!startsWithChr && clusterChr.startsWith("chr")) {
					clusterChr = clusterChr.substring(3);
				}

				if (clusterChr.equals("chrM")) {
					clusterChr = "chrMT";
				}

				// does postprocessing algorithm provides me only with the list
				// of clusters not passed filtering?!?! otherwise, is there a
				// column in the all-output showing whether the cluster passes
				// filters???

				try {
					boolean found = false;
					for (String[] clusterReference : clustersReferenceMap
							.get(clusterChr)) {
						// check for TP
						if (clusterChr.equals(clusterReference[1])
								&& (clusterStart - 5) <= Integer
										.parseInt(clusterReference[2])
								&& (clusterEnd + 5) >= Integer
										.parseInt(clusterReference[3])
								&& clusterReference[4].equals("1")
								&& clusterPassesFilter == 1) {
							// && (!onlyBoundClusters || clusterReference[4]
							// .equals("1"))
							truePositives++;
							found = true;
							break;
						} else if (clusterChr.equals(clusterReference[1])
								&& (clusterStart - 5) <= Integer
										.parseInt(clusterReference[2])
								&& (clusterEnd + 5) >= Integer
										.parseInt(clusterReference[3])
								&& clusterReference[4].equals("0")
								&& clusterPassesFilter == 0) {
							// && (onlyBoundClusters && clusterReference[4]
							// .equals("0"))

							// checks for TN
							trueNegatives++;
							found = true;
							break;
						} else if (clusterChr.equals(clusterReference[1])
								&& (clusterStart - 5) <= Integer
										.parseInt(clusterReference[2])
								&& (clusterEnd + 5) >= Integer
										.parseInt(clusterReference[3])
								&& clusterReference[4].equals("0")
								&& clusterPassesFilter == 1) {
							falsePositives++;
							found = true;
							break;

						} else if (clusterChr.equals(clusterReference[1])
								&& (clusterStart - 5) <= Integer
										.parseInt(clusterReference[2])
								&& (clusterEnd + 5) >= Integer
										.parseInt(clusterReference[3])
								&& clusterReference[4].equals("1")
								&& clusterPassesFilter == 0) {
							falseNegatives++;
							found = true;
							break;
							// } else if (clusterReference[4].equals("1")
							// && clusterPassesFilter == 1) {
							// falsePositives++;
							// } else if (clusterReference[4].equals("1")
							// && clusterPassesFilter == 0) {
							// falseNegatives++;
						}
					}
					if (!found) {
						if (clusterPassesFilter == 1) {
							falsePositives++;
						} else {
							falseNegatives++;
						}
					}
				} catch (NullPointerException e) {
					System.out
							.println("chromosome doesnt exist in reference cluster file: "
									+ clusterChr);
				}

				// calculate falsePositives and falseNegatives
				falsePositivesTotal = positives - truePositives;
				falseNegativesTotal = negatives - trueNegatives;

				processedClusters++;
			}

			inClusters.close();
			MappingLogger.getLogger().info("Calculating statistics finished.");
		} catch (IOException e) {
			MappingLogger.getLogger().error(
					"Abort due to error: IO Exception. ");
			e.printStackTrace();
		} catch (NumberFormatException e2) {
			e2.printStackTrace();
			MappingLogger.getLogger().error(
					"Abort due to error. Debug Info: start: "
							+ clusterStartString + " end: " + clusterEndString
							+ " compl: " + compl);
		} catch (SAMFormatException e3) {
			MappingLogger
					.getLogger()
					.debug("SAM Format error found. Skipping wrongly formatted read...");
		}

		MappingLogger.getLogger().debug("after cathcing exception");

		Writer output;
		try {
			output = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outFile)));

			// float precision = (float) matched_correctly / mapped_overall;
			// float sensitivity = (float) matched_correctly / allreads;
			float precision;
			float recall;
			float accuracy;
			float precisionTotal;
			float recallTotal;
			float accuracyTotal;
			// if (onlyBoundClusters) {
			// precision = (float) matched_correctly / clusters_found_bound;
			// sensitivity = (float) matched_correctly / allclusters_bound;
			// } else {
			precision = (float) truePositives
					/ (truePositives + falsePositives);
			recall = (float) truePositives / (truePositives + falseNegatives);
			accuracy = (float) (truePositives + trueNegatives)
					/ (truePositives + trueNegatives + falsePositives + falseNegatives);

			precisionTotal = (float) truePositives
					/ (truePositives + falsePositivesTotal);
			recallTotal = (float) truePositives
					/ (truePositives + falseNegativesTotal);
			accuracyTotal = (float) (truePositives + trueNegatives)
					/ (positives + negatives);

			// }
			// float mapped_perc = (float) clusters_found_overall / allclusters;
			// float f1_score = (float) 2 * (precision * sensitivity)
			// / (precision + sensitivity);

			MappingLogger.getLogger().info("Writing statistics to file.");
			output.write("\ntruePositives:\t" + truePositives
					+ "\ntrueNegatives:\t" + trueNegatives
					+ "\nfalsePositives:\t" + falsePositives
					+ "\nfalseNegatives:\t" + falseNegatives
					+ "\nfalsePositivesTotal:\t" + falsePositivesTotal
					+ "\nfalseNegativesTotal\t" + falseNegativesTotal
					+ "\nrecall:\t\t" + recall + "\nprecision:\t" + precision
					+ "\naccuracy:\t" + accuracy + "\nrecallTotal:\t"
					+ recallTotal + "\nprecisionTotal\t" + precisionTotal
					+ "\naccuracyTotal:\t" + accuracyTotal);
			// MappingLogger.getLogger().info(
			// "Precision:\t\t\t" + precision + "\nrecall:\t\t\t\t"
			// + recall + "\naccuracy:\t\t\t" + accuracy);
			MappingLogger.getLogger().info(
					"\ntruePositives:\t" + truePositives + "\ntrueNegatives:\t"
							+ trueNegatives + "\nfalsePositives:\t"
							+ falsePositives + "\nfalseNegatives:\t"
							+ falseNegatives + "\nfalsePositivesTotal:\t"
							+ falsePositivesTotal + "\nfalseNegativesTotal\t"
							+ falseNegativesTotal + "\nrecall:\t\t" + recall
							+ "\nprecision:\t" + precision + "\naccuracy:\t"
							+ accuracy + "\nrecallTotal:\t" + recallTotal
							+ "\nprecisionTotal\t" + precisionTotal
							+ "\naccuracyTotal:\t" + accuracyTotal);

			output.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
