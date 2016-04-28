package utils.benchmarking;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;

import main.MappingLogger;

/**
 * Evaluates the alignment accuracy of a respective read aligner performed on
 * simulated PAR-CLIP read dataset.
 * 
 * @author akloetgen
 * 
 */
public class ValidateBenchmarkStatisticsPARCLIP {

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
	public void calculateBenchmarkStatistics(String mappingFileName,
			String outStatistics, String readsFile, boolean onlyBoundClusters) {

		final SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
				.validationStringency(ValidationStringency.LENIENT);
		SamReader samFileReader = factory.open(new File(mappingFileName));
		File outFile = new File(outStatistics);
		int readsProcessed = 0;
		int allreads = 0;

		// reads marked as bound and correctly aligned?
		int truePositives = 0;
		// reads marked as NOT bound and correctly aligned?
		int trueNegatives = 0;
		// reads marked as bound and falsely aligned
		int falsePositives = 0;
		// reads marked as not bound and falsely aligned
		int falseNegatives = 0;
		// reads marked as bound and falsely/not aligned
		// int falsePositivesTotal = 0;
		// reads marked as not bound and falsely/not aligned
		// int falseNegativesTotal = 0;
		int positives = 0;
		int negatives = 0;

		String readStartString = "";
		String readEndString = "";
		String clusterBound = "";
		String compl = "";

		boolean startsWithChr = false;

		try {

			BufferedReader in = new BufferedReader(new FileReader(readsFile));
			String line = null;
			while ((line = in.readLine()) != null) {
				allreads++;
				if (line.startsWith("@SEQ_ID")) {
					String[] splittedHeader = line.split("\\|");
					if (splittedHeader[5].split("-")[0].equals("1")) {
						// allreads_bound++;
						positives++;
					} else if (splittedHeader[5].split("-")[0].equals("0")) {
						negatives++;
					}
					// if (splittedHeader[2].startsWith("chr")) {
					//
					// }

				}

			}
			in.close();
			if (allreads % 4 != 0) {
				MappingLogger.getLogger().error(
						"FASTQ cannot be divided by 4! Exiting...");
				System.exit(0);
			}
			allreads /= 4;

			for (SAMRecord readHit : samFileReader) {
				if (readsProcessed % 100000 == 0) {
					MappingLogger.getLogger().debug(
							"readsProcessed: " + readsProcessed);
				}
				// fastq-header = information about origin in simulated data???
				// if yes: parse header, compare, output, ready

				String[] splittedReadName = readHit.getReadName().split("\\|");
				compl = readHit.getReadName();
				String readChr = splittedReadName[2];
				readStartString = splittedReadName[3];
				readEndString = splittedReadName[4];
				clusterBound = splittedReadName[5].split("-")[0];

				// if (onlyBoundClusters) {
				// clusterBound = splittedReadName[5].split("-")[0];
				// if (clusterBound.equals("1")) {
				// mapped_overall_bound++;
				// }
				// }

				int readStart = Integer.parseInt(readStartString);
				int readEnd = Integer.parseInt(readEndString);

				String chr = readHit.getReferenceName();
				if (chr.startsWith("chr")) {
					startsWithChr = true;
				}
				if (startsWithChr && !readChr.startsWith("chr")) {
					readChr = "chr" + readChr;
				} else if (!startsWithChr && readChr.startsWith("chr")) {
					readChr = readChr.substring(3);
				}
				// System.out.println("ref=" + chr + "; readChr=" + readChr);

				if (readChr.equals("chrM")) {
					readChr = "chrMT";
				}

				if (readChr.equals(chr)
						&& (readStart - 5) <= readHit.getAlignmentStart()
						&& (readEnd + 5) >= readHit.getAlignmentEnd()
						&& clusterBound.equals("1")) {
					// CEHCK THIS AGAIN ON MAPSPLICE; WHETHER READ AND NOT
					// READHIT SHOULD BE CHECKED AGAIN!=?!==!=!=!?!?!??!?!F
					truePositives++;
					// continue;
				} else if (readChr.equals(chr)
						&& (readStart - 5) <= readHit.getAlignmentStart()
						&& (readEnd + 5) >= readHit.getAlignmentEnd()
						&& clusterBound.equals("0")) {
					trueNegatives++;
					// break;
				}

				readsProcessed++;

			}

			falsePositives = positives - truePositives;
			falseNegatives = negatives - trueNegatives;

			System.out.println("TP=" + truePositives + "; TN=" + trueNegatives
					+ "; FP=" + falsePositives + "; FN=" + falseNegatives);

			samFileReader.close();
			MappingLogger.getLogger().info("Calculating statistics finished.");
		} catch (IOException e) {
			MappingLogger.getLogger().error(
					"Abort due to error: IO Exception. ");
			e.printStackTrace();
		} catch (NumberFormatException e2) {
			e2.printStackTrace();
			MappingLogger.getLogger().error(
					"Abort due to error. Debug Info: start: " + readStartString
							+ " end: " + readEndString + " compl: " + compl);
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
			// float precision;
			// float sensitivity;

			float precision = (float) truePositives
					/ (truePositives + falsePositives);
			float recall = (float) truePositives
					/ (truePositives + falseNegatives);
			float accuracy = (float) (truePositives + trueNegatives)
					/ (positives + negatives);

			// if (onlyBoundClusters) {
			// precision = (float) matched_correctly / mapped_overall_bound;
			// sensitivity = (float) matched_correctly / allreads_bound;
			// } else {
			// precision = (float) matched_correctly / mapped_overall;
			// sensitivity = (float) matched_correctly / allreads;
			// }
			// float mapped_perc = (float) mapped_overall / allreads;
			// float f1_score = (float) 2 * (precision * sensitivity)
			// / (precision + sensitivity);

			int alignedCorrectly = truePositives + trueNegatives;
			// float alignedOverall =

			MappingLogger.getLogger().info("Writing statistics to file.");
			output.write("matched correctly:\t" + alignedCorrectly
					+ "\nreadsProcessed:\t" + readsProcessed + "\nall reads:\t"
					+ allreads + "\nprecision:\t" + precision + "\nrecall:\t"
					+ recall + "\naccuracy:\t" + accuracy);
			MappingLogger.getLogger().info(
					"Precision:\t\t\t" + precision + "\nRecall:\t\t\t" + recall
							+ "\nAccuracy:\t\t\t" + accuracy);
			MappingLogger.getLogger().debug(
					"matched correctly:\t" + alignedCorrectly
							+ "\nreadsProcessed:\t" + readsProcessed
							+ "\nall reads:\t" + allreads + "\nprecision:\t"
							+ precision + "\nrecall:\t" + recall
							+ "\naccuracy:\t" + accuracy);

			output.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
