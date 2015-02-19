package main;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Not used???
 * 
 * @author akloetgen
 * 
 */
public class Help {

	public void printUsageHelp() {
		try {
			InputStream fileStream = this.getClass().getResourceAsStream(
					"help-file");
			if (fileStream == null) {
				MappingLogger.getLogger().error(
						"File stream empty for Help print. Abort.");
			}
			BufferedReader reader = new BufferedReader(new InputStreamReader(
					fileStream));
			String line;

			while ((line = reader.readLine()) != null) {
				System.out.println(line);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void printProgramModeHelp() {
		MappingLogger.getLogger().info(
				"Usage: java -jar parma.jar <PROGRAMMODE> [options]"
						+ System.getProperty("line.separator"));
		MappingLogger
				.getLogger()
				.info("Following programmodes are available through the PARMA toolkit:");
		MappingLogger
				.getLogger()
				.info("\tmap\t\tstarts the read mapping of a read set against a reference sequence");
		MappingLogger
				.getLogger()
				.info("\tcomb\t\tcombines genomic and transcriptomic alignment files");
		MappingLogger.getLogger().info(
				"\terror\t\tcalculates the errorprofile of a given alignment");
		MappingLogger
				.getLogger()
				.info("\tclust\t\tclusters aligned reads together if they overlap by at least 5 bases");
		MappingLogger.getLogger().info(
				"\tsimluate\tsimulates a realistic PAR-CLIP read dataset");
		MappingLogger
				.getLogger()
				.info("\tbenchmark\tcalculates mapping statsistics of a mapping against a simulated PAR-CLIP dataset");
		MappingLogger.getLogger().info(
				"\tsetup\t\tenables to set some PARMA-tk specific parameters");
	}

	public void printMappingToolHelp() {
		MappingLogger.getLogger().info(
				"Usage: java -jar parma.jar mapping [options]"
						+ System.getProperty("line.separator"));
		MappingLogger.getLogger().info("Following options are available");
		MappingLogger.getLogger().info(
				"\t-q FILE*\t\tread file in FASTQ-format");
		MappingLogger
				.getLogger()
				.info("\t-r FILE*\t\treference genomic sequence file, must be indexed");
		MappingLogger.getLogger().info(
				"\t-t FILE\t\treference transcriptomic sequence file");
		MappingLogger.getLogger().info("\t-o STRING*\toutput file prefix");
		MappingLogger
				.getLogger()
				.info("\t-l INT\t\tmaximal read length given [can probably be calculated??] [default: 101]");
		MappingLogger.getLogger().info("\t-p INT\t\tnum threads [default: 1]");
		MappingLogger
				.getLogger()
				.info("\t--gm INT\tmapping quality filter for genomic mapping [default: 10]");
		MappingLogger
				.getLogger()
				.info("\t--tm INT\tmapping quality filter for transcriptomic mapping [default: 1]");
		MappingLogger
				.getLogger()
				.info("\t--refine\tuses error-profile derived from selected mapping vs. genome to apply PARMA with the calculated error-profile");
		MappingLogger
				.getLogger()
				.info("\t--ref-refine FILE\tReference file for PARMA algorithm, if --refine option was activated [default: same as -r file]");
		MappingLogger
				.getLogger()
				.info("\t--mode\t\tdefines mapping algorithm used: BT2, BWA, PARMA, USER [default: BWA]");
		MappingLogger
				.getLogger()
				.info("\t--parma-mm	ONLY IF \"--mode parma\": number of average mismatches while applying PARMA algorithm [default: 2]");
		MappingLogger
				.getLogger()
				.info("\t--parma-ep	ONLY IF \"--mode parma\": filename of the error profile file");
		MappingLogger
				.getLogger()
				.info("\t--parma-indel	ONLY IF \"--mode parma\": filename of the error profile file");
		MappingLogger.getLogger().debug(
				"SOME VALUES ARE STILL MISSING IN THIS HELP!!!!");
		System.out.println();
		MappingLogger.getLogger()
				.info("Options marked with a * are requiered.");

	}

	public void printCombineToolHelp() {
		MappingLogger
				.getLogger()
				.info("Usage: java -jar parma.jar combine GENOMIC_MAPPING_FILE TRANSCRIPTOMIC_MAPPING_FILE OUTPUT_FILE");
	}

	public void printBenchmarkToolHelp() {
		MappingLogger
				.getLogger()
				.info("Usage: java -jar parma.jar benchmark MAPPING_FILE OUT_STATISTICS_FILE READS_FILE");
	}

	public void printClusterToolHelp() {
		MappingLogger
				.getLogger()
				.info("Usage: java -jar parma.jar clustering MAPPING_FILE REFERENCE_FILE OUTPUT_FILE SNP_FILE MIN_COVERAGE");
	}

	public void printErrorprofileToolHelp(boolean isQualityCalc,
			boolean showErrorPlot) {
		MappingLogger
				.getLogger()
				.info("Usage: java -jar parma.jar error MAPPING_FILE REFERENCE_FILE MAX_READ_LENGTH [options]"
						+ System.getProperty("line.separator"));
		MappingLogger.getLogger().info("Additional options are:");
		MappingLogger.getLogger().info(
				"\t-q\ttrue/false: calculates base calling quality distribution -> SLOW ["
						+ isQualityCalc + "]");
		MappingLogger.getLogger().info(
				"\t-p\ttrue/false: calculates plot for error profile. Requires X11 terminal! ["
						+ showErrorPlot + "]");
	}

	public void printWrongInputMessage(String message) {
		MappingLogger.getLogger().error(message);
		printUsageHelp();
	}

	public void printSetupHelp() {
		MappingLogger
				.getLogger()
				.info("NOT AVAILABLE IN THE CURRENT VERSION! PLEASE MANUALLY SET THE PATHS TO THE ALIGNERS!!");
		MappingLogger.getLogger().info(
				"Usage: java -jar parma.jar setup --ALIGNER PATH_TO_ALIGNER"
						+ System.getProperty("line.separator"));
		MappingLogger
				.getLogger()
				.info("Currently supported options for aligners are: parma, bwa, bowtie2");
	}

	public void printSimulateHelp() {
		MappingLogger
				.getLogger()
				.info("Usage: java -jar parma.jar simluate TRANSCRIPT_FILE OUTPUT_PREFIX ERROR_PROFILE T2C_PROFILE T2C_POSITIONS_PROFILE QUALITY_DIST INDEL_PROFILE RBP_BOUND [options]"
						+ System.getProperty("line.separator"));
		MappingLogger
				.getLogger()
				.info("\tTRANSCRIPT_FILE\t\tfasta file containing transcript sequences on which PAR-CLIP reads are simulated");
		MappingLogger
				.getLogger()
				.info("\tOUTPUT_PREFIX\t\tprefix for the output files. A fastq file as well as a log file are created");
		MappingLogger
				.getLogger()
				.info("\tERROR_PROFILE\t\tfile name to the error profile with mismatch probabilities");
		MappingLogger
				.getLogger()
				.info("\tT2C_PROFILE\t\tfile name to T-C conversion rates, ordered by the majority of occurence within a cluster");
		MappingLogger
				.getLogger()
				.info("\tT2C_POSITIONS_PROFILE\tfile name to T-C conversion site specific probabilites. Signifies each a probability how likely a T-C conversion is at the respective cluster position");
		MappingLogger
				.getLogger()
				.info("\tQUALITY_DIST\t\tfile name to base calling qualities; PHRED64 scaling");
		MappingLogger
				.getLogger()
				.info("\tINDEL_PROFILE\t\tfile name to the indel profile, giving positions-specific probabilites for insertions and deletions");
		MappingLogger
				.getLogger()
				.info("\tRBP_BOUND\t\tdouble value specifying the fraction of RBP-bound clusters that are simulated"
						+ System.getProperty("line.separator"));
		MappingLogger.getLogger().info("Additional options are:");
		MappingLogger
				.getLogger()
				.info("\t-I\tpath to cpan Math::Random library. Necessary, if perl does not find the library.");
	}
}
