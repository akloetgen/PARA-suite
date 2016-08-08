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
 * 
 * @author akloetgen
 * @deprecated
 */
public class ValidateBenchmarkStatistics {

	public static void main(String[] args) {
		MappingLogger.getLogger().info(
				"Calculating statistics for mapped "
						+ "reads against simulated data-set");
		String mappingFileName = null;
		String outStatistics = null;
		String readsFile = null;

		try {
			if (args[0].equals("-h") || args[0].equals("--help")) {
				MappingLogger
						.getLogger()
						.error("Usage:\njava -jar evaluate_benchmark.jar "
								+ "alignmentFileName outputFileStatistics rawReaedsFile");
				System.exit(0);
			}
			mappingFileName = args[0];
			outStatistics = args[1];
			readsFile = args[2];
		} catch (ArrayIndexOutOfBoundsException e) {
			MappingLogger
					.getLogger()
					.error("Usage:\njava -jar evaluate_benchmark.jar "
							+ "alignmentFileName outputFileStatistics rawReaedsFile");
			System.exit(0);
		}

		// EVTL ZUSÄTZLICHE OPTIONEN EINBAUEN! KÖNNTE SCHNELLER WERDEN!
		final SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
				.validationStringency(ValidationStringency.LENIENT);
		SamReader samFileReader = factory.open(new File(mappingFileName));

		// File mappingFile = new File(mappingFileName);
		// SamReader samFileReader = new SAMFileReader(mappingFile);
		File outFile = new File(outStatistics);
		int matched_correctly = 0;
		int mapped_overall = 0;
		int allreads = 0;

		String readStartString = "";
		String readEndString = "";
		String compl = "";

		try {

			BufferedReader in = new BufferedReader(new FileReader(readsFile));

			while (in.readLine() != null) {
				allreads++;
			}
			in.close();
			if (allreads % 4 != 0) {
				MappingLogger.getLogger().error(
						"FASTQ cannot be divided by 4! Exiting...");
				System.exit(0);
			}
			allreads /= 4;

			// int rounds = 10;

			for (SAMRecord readHit : samFileReader) {
				if (mapped_overall % 100000 == 0) {
					MappingLogger.getLogger().debug(
							"mapped_overall: " + mapped_overall);
				}
				// fastq-header = information about origin in simulated data???
				// if yes: parse header, compare, output, ready
				// System.out.println("read-name: " + read.getReadName());

				// chr:start-end
				// System.out.println("ref-name: " + read.getReferenceName());
				// System.out.println("alstart-name: " +
				// read.getAlignmentStart());
				// System.out.println("alstop-name: " + read.getAlignmentEnd());

				String[] splittedReadName = readHit.getReadName().split("\\|");
				compl = readHit.getReadName();
				String readChr = splittedReadName[1];
				readStartString = splittedReadName[2];

				// following splitting is for ART simulated data
				readEndString = splittedReadName[3].split("-")[0];

				// FILTER FOR RBP-BOUND CLUSTERS! ONLY THOSE ARE OF INTEREST AND
				// ARE DENOTED BY A 1 IN THE NEXT TO LAST SPLIT

				// following splitting is for own simulated data
				// readEndString = splittedReadName[3].split(";")[0];
				int readStart = Integer.parseInt(readStartString);
				int readEnd = Integer.parseInt(readEndString);

				String chr = readHit.getReferenceName();
				if (chr.startsWith("chr")) {
					chr = chr.substring(3);
				}

				if (readChr.equals(chr)
						&& (readStart - 50) <= readHit.getAlignmentStart()
						&& (readEnd + 50) >= readHit.getAlignmentEnd()) {

					// CEHCK THIS AGAIN ON MAPSPLICE; WHETHER READ AND NOT
					// READHIT SHOULD BE CHECKED AGAIN!=?!==!=!=!?!?!??!?!F

					// output.write("MATCH\n");
					// if
					// (readsAlreadyMappedCorrectly.get(readHit.getReadName())
					// == null) {
					matched_correctly++;
					// readsAlreadyMappedCorrectly.put(readHit.getReadName(),
					// readHit);
					//
					// } else {
					// MappingLogger.getLogger().debug(
					// "prev: start: "
					// + readsAlreadyMappedCorrectly.get(
					// readHit.getReadName())
					// .getAlignmentStart()
					// + "; end: "
					// + readsAlreadyMappedCorrectly.get(
					// readHit.getReadName())
					// .getAlignmentEnd()
					// + "; chr: "
					// + readsAlreadyMappedCorrectly.get(
					// readHit.getReadName())
					// .getReferenceName()
					// + "; readName: "
					// + readsAlreadyMappedCorrectly.get(
					// readHit.getReadName())
					// .getReadName()
					// + "; MAPQ: "
					// + readsAlreadyMappedCorrectly.get(
					// readHit.getReadName())
					// .getMappingQuality());
					//
					// MappingLogger.getLogger().debug(
					// "current: start: "
					// + readHit.getAlignmentStart()
					// + "; end: " + readHit.getAlignmentEnd()
					// + "; chr: "
					// + readHit.getReferenceName()
					// + "; readName: "
					// + readHit.getReadName() + "; MAPQ: "
					// + readHit.getMappingQuality());
					//
					// }
					// } else {
					// output.write("NOT MATCHED\n");
					// } else {
					// MappingLogger.getLogger().debug(
					// readHit.getReadName() + " : "
					// + readHit.getReferenceName() + ":"
					// + readHit.getAlignmentStart() + "-"
					// + readHit.getAlignmentEnd());
					// counter++;
					// } else {
					// MappingLogger.getLogger().debug(
					// "wrong: MAPQ: " + readHit.getMappingQuality()
					// + readHit.getReadName() + "; "
					// + readHit.getReferenceName() + ": "
					// + readHit.getAlignmentStart() + "-"
					// + readHit.getAlignmentEnd());
					// if (readHit.getReadFailsVendorQualityCheckFlag()) {
					// MappingLogger.getLogger().debug(
					// "wrong: " + readHit.getReadName() + "; "
					// + readHit.getReferenceName() + ": "
					// + readHit.getAlignmentStart() + "-"
					// + readHit.getAlignmentEnd());
					// }
				}
				// if (counter == 20) {
				// System.exit(0);
				// }

				// if (readsAlreadyMapped.get(readHit.getReadName()) == null) {
				mapped_overall++;
				// TOOOOOOOOOOOOOOOOO
				// SLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// readsAlreadyMapped.put(readHit.getReadName(), readHit);
				// }

				// if (rounds == 0) {
				// break;
				// }
				// rounds--;

			}

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
			// e3.printStackTrace();
			MappingLogger
					.getLogger()
					.debug("SAM Format error found. Skipping wrongly formatted read...");
		}

		MappingLogger.getLogger().debug("after cathcing exception");

		Writer output;
		try {
			output = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outFile)));

			float precision = (float) matched_correctly / mapped_overall;
			float sensitivity = (float) matched_correctly / allreads;
			float mapped_perc = (float) mapped_overall / allreads;
			float f1_score = (float) 2 * (precision * sensitivity)
					/ (precision + sensitivity);

			MappingLogger.getLogger().info("Writing statistics to file.");
			output.write("matched correctly:\t" + matched_correctly
					+ "\nmapped overall:\t" + mapped_overall + "\nall reads:\t"
					+ allreads + "\nmapped overall percentage:\t" + mapped_perc
					+ "\nprecision:\t" + precision + "\nsensitivity:\t"
					+ sensitivity + "\nf1 score:\t" + f1_score);
			MappingLogger.getLogger().debug(
					"matched correctly:\t\t" + matched_correctly
							+ "\nmapped overall:\t\t\t" + mapped_overall
							+ "\nall reads:\t\t\t" + allreads
							+ "\nmapped overall percentage:\t" + mapped_perc
							+ "\nprecision:\t\t\t" + precision
							+ "\nsensitivity:\t\t\t" + sensitivity
							+ "\nf1 score:\t\t\t" + f1_score);

			output.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
