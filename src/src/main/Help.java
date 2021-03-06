package main;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import enums.ProgramMode;

/**
 * Class for printing help information.
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

	public void printProgramInfo(ProgramMode mode, String version, String author) {
		System.out.println("Program:\tPARA-suite - " + mode.toString());
		System.out.println("Version:\t" + version);
		System.out.println("Author:\t\t" + author
				+ System.getProperty("line.separator"));
	}

	public void printProgramModeHelp() {
		System.out
				.println("Usage: java -jar parasuite.jar <PROGRAMMODE> [options]"
						+ System.getProperty("line.separator"));
		System.out
				.println("Following programmodes are available through the PARA-suite:");
		System.out
				.println("\tmap\t\tstarts the read mapping of a read set against a reference sequence");
		System.out
				.println("\tcomb\t\tcombines genomic and transcriptomic alignment files");
		System.out
				.println("\terror\t\tcalculates the errorprofile of a given alignment");
		System.out
				.println("\tclust\t\tclusters aligned reads together if they overlap by at least 5 bases");
		System.out
				.println("\tsimulate\tsimulates a realistic PAR-CLIP read dataset");
		System.out
				.println("\tbenchmark\tcalculates mapping statsistics of a mapping against a simulated PAR-CLIP dataset");
		System.out
				.println("\tbenchmarkClusters\tcalculates cluster statsistics of detected binding sites from a simulated PAR-CLIP dataset");
		System.out
				.println("\tsetup\t\tenables to set some PARA-suite specific parameters");
		System.out
				.println("\textract\t\textracts weakly aligned/unaligned reads from a BAM file");
	}

	public void printMappingToolHelp() {
		System.out.println("Usage: java -jar parasuite.jar map [options]"
				+ System.getProperty("line.separator"));
		System.out.println("Following options are available");
		System.out.println("\t-q FILE*\tread file in FASTQ-format");
		System.out
				.println("\t-r FILE*\treference genomic sequence file, must be indexed");
		System.out
				.println("\t-t FILE\t\treference transcriptomic sequence file");
		System.out.println("\t-o STRING*\toutput file prefix");
		System.out
				.println("\t-l INT\t\tmaximal read length given [can probably be calculated??] [default: 101]");
		System.out.println("\t-p INT\t\tnum threads [default: 1]");
		System.out
				.println("\t--gm INT\tmapping quality filter for genomic mapping [default: 10]");
		System.out
				.println("\t--tm INT\tmapping quality filter for transcriptomic mapping [default: 1]");
		System.out
				.println("\t--refine\tuses error-profile derived from selected mapping vs. genome to apply PARA-suite with the calculated error-profile");
		System.out
				.println("\t--ref-refine FILE\tReference file for PARA-suite algorithm, if --refine option was activated [default: same as -r file]");
		System.out
				.println("\t--unaligned\textracts unaligned reads after MAPQ filtering into a separate file");
		System.out
				.println("\t--mode\t\tdefines mapping algorithm used: BT2, BWA, PARAsuite, USER [default: BWA]");
		System.out
				.println("\t--parauiste-mm	ONLY IF \"--mode parasuite\": number of average mismatches while applying PARA-suite algorithm [default: 2]");
		System.out
				.println("\t--parauiste-ep	ONLY IF \"--mode parasuite\": filename of the error profile file");
		System.out
				.println("\t--parauiste-indel	ONLY IF \"--mode parasuite\": filename of the error profile file");
		System.out
				.println("\t-c COMMAND\tsets the user aligner command in \"\". Use INPUT, REFERENCE, OUTPUT and THREADS as place holders."
						+ " The TK will exchange those place holders with the respective options.");
		System.out.println("\t-o FILE\t\tfilename for the logging file.");
		System.out
				.println("\t--unaligned\t\tkeeps unaligned/weakly aligned reads in the BAM file");
		System.out.println();
		System.out.println("Options marked with a * are requiered.");

	}

	public void printCombineToolHelp() {
		System.out.println("Usage: java -jar parasuite.jar comb");
		System.out.println("\t-g FILE\talignment vs. genomic reference file");
		System.out
				.println("\t-t FILE\talignment vs. transcriptomic reference file");
		System.out.println("\t-o FILE\toutput file");
	}

	public void printBenchmarkToolHelp() {
		System.out
				.println("Usage: java -jar parasuite.jar benchmark MAPPING_FILE OUT_STATISTICS_FILE READS_FILE");
	}

	public void printBenchmarkClustersToolHelp() {
		System.out.println("Usage: java -jar parasuite.jar benchmarkClusters "
				+ "CLUSTERS_FILE OUT_STATISTICS_FILE CLUSTER_REFERENCE_FILE");

	}

	public void printClusterToolHelp() {
		System.out
				.println("Usage: java -jar parasuite.jar clust MAPPING_FILE REFERENCE_FILE OUTPUT_FILE SNP_FILE MIN_COVERAGE");
	}

	public void printErrorprofileToolHelp(boolean isQualityCalc,
			boolean showErrorPlot) {
		System.out
				.println("Usage: java -jar parasuite.jar error MAPPING_FILE REFERENCE_FILE MAX_READ_LENGTH [options]"
						+ System.getProperty("line.separator"));
		System.out.println("Additional options are:");
		System.out
				.println("\t-q\ttrue/false: calculates base calling quality distribution -> SLOW ["
						+ isQualityCalc + "]");
		System.out
				.println("\t-p\ttrue/false: calculates plot for error profile. Requires X11 terminal! ["
						+ showErrorPlot + "]");
	}

	public void printWrongInputMessage(String message) {
		MappingLogger.getLogger().error(message);
		printUsageHelp();
	}

	public void printSetupHelp() {
		System.out
				.println("Usage: java -jar parasuite.jar setup --ALIGNER PATH_TO_ALIGNER"
						+ System.getProperty("line.separator"));
		System.out
				.println("Currently supported options for aligners are: --parasuite, --bwa, --bowtie2");
	}

	public void printSimulateHelp() {
		System.out
				.println("Usage: java -jar parasuite.jar simulate TRANSCRIPT_FILE OUTPUT_PREFIX ERROR_PROFILE T2C_PROFILE T2C_POSITIONS_PROFILE QUALITY_DIST INDEL_PROFILE RBP_BOUND [options]"
						+ System.getProperty("line.separator"));
		System.out
				.println("\tTRANSCRIPT_FILE\t\tfasta file containing transcript sequences on which PAR-CLIP reads are simulated");
		System.out
				.println("\tOUTPUT_PREFIX\t\tprefix for the output files. A fastq file as well as a log file are created");
		System.out
				.println("\tERROR_PROFILE\t\tfile name to the error profile with mismatch probabilities");
		System.out
				.println("\tT2C_PROFILE\t\tfile name to T-C conversion rates, ordered by the majority of occurence within a cluster");
		System.out
				.println("\tT2C_POSITIONS_PROFILE\tfile name to T-C conversion site specific probabilites. Signifies each a probability how likely a T-C conversion is at the respective cluster position");
		System.out
				.println("\tQUALITY_DIST\t\tfile name to base calling qualities; PHRED64 scaling");
		System.out
				.println("\tINDEL_PROFILE\t\tfile name to the indel profile, giving positions-specific probabilites for insertions and deletions");
		System.out
				.println("\tRBP_BOUND\t\tdouble value specifying the fraction of RBP-bound clusters that are simulated"
						+ System.getProperty("line.separator"));
		System.out.println("Additional options are:");
		System.out
				.println("\t-I\tpath to cpan Math::Random library. Necessary, if perl does not find the library.");
	}

	public void printExtractHelp(int mapqThreshold) {
		System.out.println("Usage: java -jar parasuite.jar extract [options]"
				+ System.getProperty("line.separator"));
		System.out.println("\t-i FILE\t\tBAM file containing read alignments");
		System.out
				.println("\t-o FILE\t\tFASTQ filename where weak mappings are saved to");
		System.out
				.println("\t-t INT\t\tMAPQ used as threshold to extract weak read alignments ["
						+ mapqThreshold + "]");
	}
}
