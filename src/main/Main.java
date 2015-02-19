package main;

import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.net.URL;
import java.net.URLConnection;
import java.util.LinkedList;
import java.util.List;

import mapping.BWAMapping;
import mapping.BowtieMapping;
import mapping.Mapping;
import mapping.MappingErrorException;
import mapping.PARMAMapping;
import mapping.UserMapping;
import utils.benchmarking.ValidateBenchmarkStatisticsPARCLIP;
import utils.errorprofile.ErrorProfiling;
import utils.pileupclusters.PileupClusters;
import utils.postprocessing.CombineGenomeTranscript;
import utils.postprocessing.ExtractWeakMappingReads;
import enums.MappingMode;
import enums.PARMAPropertiesEnum;
import enums.ProgramMode;

/**
 * 
 * Main class as entry point for PARMA pipeline. Can execute several sub-utils
 * of the PARMA project.
 * 
 * @author akloetgen
 * 
 */
public class Main {

	private static final String VERSION = "0.5 alpha";
	private static final String AUTHOR = "Andreas Kloetgen";

	public static void main(String[] args) {
		MappingLogger.getLogger().info("Program:\t\tPARMA toolkit");
		MappingLogger.getLogger().info("Version:\t\t" + VERSION);
		MappingLogger.getLogger().info(
				"Author:\t\t" + AUTHOR + System.getProperty("line.separator"));
		// MappingLogger.getLogger().info(
		// "Further reading:\t" + "URL" + System.getProperty("line.separator"));

		Help help = new Help();
		ProgramMode mode = null;

		if (args.length == 0 || args[0].equals("-h")
				|| args[0].equals("--help")) {
			help.printProgramModeHelp();
			System.exit(0);
		}

		switch (args[0]) {
		case "map":
			mode = ProgramMode.MAPPING;
			break;
		case "comb":
			mode = ProgramMode.COMBINE;
			break;
		case "benchmark":
			mode = ProgramMode.BENCHMARK;
			break;
		case "error":
			mode = ProgramMode.ERRORPROFILE;
			break;
		case "clust":
			mode = ProgramMode.CLUSTERING;
			break;
		case "simulate":
			mode = ProgramMode.SIMULATE;
			break;
		case "setup":
			mode = ProgramMode.SETUP;
			break;
		default:
			MappingLogger
					.getLogger()
					.error("Programmode \""
							+ args[0]
							+ "\" not found. Please consider the following help-message for a correct usage.");
			help.printProgramModeHelp();
			System.exit(0);
			break;
		}

		if (mode.equals(ProgramMode.MAPPING)) {
			// ############## MAPPING TOOL ##############
			try {
				String mappingFileName = null;
				String referenceFileName = null;
				String referenceRefineFileName = null;
				String transcriptFileName = null;
				String outputPrefix = null;
				String readFileName = null;
				String profileFileName = null;
				String indelFileName = null;
				boolean isMapping = false;
				boolean isTranscriptMapping = false;
				int threads = 1;
				int maxReadLength = 101;
				int mappingQualityFilterGenomic = 10;
				int mappingQualityFilterTranscript = 1;
				String alignmentMode = MappingMode.BWA.toString();
				String userAlignmentCommand = "";
				boolean isRefine = true;
				boolean isFirstMapping = true;
				String additionalOptions = "";
				String parmaMismatches = "-1";
				String bwaMismatches = "2";

				if (args.length == 1) {
					help.printMappingToolHelp();
					System.exit(0);
				}

				for (int i = 1; i < args.length; i++) {
					switch (args[i]) {
					case "-q":
						i++;
						readFileName = args[i];
						isMapping = true;
						break;
					case "-r":
						i++;
						referenceFileName = args[i];
						referenceRefineFileName = referenceFileName;
						break;
					case "-t":
						i++;
						transcriptFileName = args[i];
						isTranscriptMapping = true;
						break;
					case "-p":
						i++;
						threads = Integer.parseInt(args[i]);
						break;
					case "-o":
						i++;
						outputPrefix = args[i];
						MappingLogger.setLoggingFile(outputPrefix + ".log");
						break;
					case "-l":
						// should be calculated automatically..
						i++;
						maxReadLength = Integer.parseInt(args[i]);
						break;
					case "--gm":
						i++;
						mappingQualityFilterGenomic = Integer.parseInt(args[i]);
						break;
					case "--tm":
						i++;
						mappingQualityFilterTranscript = Integer
								.parseInt(args[i]);
						break;
					case "--mode":
						i++;
						boolean found = false;
						for (MappingMode mappingMode : MappingMode.values()) {
							if (args[i].toUpperCase().equals(
									mappingMode.toString())) {
								alignmentMode = args[i];
								found = true;
							}
						}
						if (!found) {
							MappingLogger
									.getLogger()
									.error("Wrong mapping mode set. Please choose "
											+ "one of the following: bt2, bwa, parma, user");
						}
						break;
					case "--parma-mm":
						i++;
						parmaMismatches = args[i];
						break;
					case "--parma-ep":
						// COULD BE DELETED AS THIS IS A TK AND NOT A PIPELINE
						i++;
						isFirstMapping = false;
						profileFileName = args[i];
						break;
					case "--parma-indel":
						// COULD BE DELETED AS THIS IS A TK AND NOT A PIPELINE
						i++;
						indelFileName = args[i];
						break;
					case "--bwa-mm":
						i++;
						bwaMismatches = args[i];
						break;
					case "-c":
						i++;
						userAlignmentCommand = args[i];
						break;
					case "--refine":
						// COULD BE DELETED AS THIS IS A TK AND NOT A PIPELINE
						isRefine = true;
						break;
					case "--ref-refine":
						// if the first mapping is not made by BWA, a second
						// reference file for the bwa parma mod is necessary.
						i++;
						referenceRefineFileName = args[i];
						break;
					default:
						MappingLogger
								.getLogger()
								.error("Parameter not found: "
										+ args[i]
										+ ". Please consider the following help:");
						help.printMappingToolHelp();
					}
				}
				if (readFileName == null || outputPrefix == null
						|| referenceFileName == null) {
					MappingLogger.getLogger().error(
							"No valid argument for -q, -r or -o given in input. Please consider"
									+ " help for correct usage:");
					help.printMappingToolHelp();
					System.exit(0);
				}

				if (isMapping) {
					String unalignedGenomicFileName = outputPrefix
							+ ".unaligned.fastq";
					String outputPrefixBowtieGenomic = outputPrefix + "."
							+ alignmentMode + "-genomic";
					String genomicMappingFileName = outputPrefixBowtieGenomic
							+ ".bam";
					String genomicMappingFileNameNew = genomicMappingFileName
							+ ".new";
					Mapping mapping = null;
					ExtractWeakMappingReads extract;

					if (alignmentMode.equals(MappingMode.BT2.toString())) {
						mapping = new BowtieMapping();
					} else if (alignmentMode.equals(MappingMode.BWA.toString())) {
						mapping = new BWAMapping();
						additionalOptions = bwaMismatches;
					} else if (alignmentMode.equals(MappingMode.PARMA
							.toString())) {
						mapping = new PARMAMapping();
						additionalOptions = parmaMismatches;
					} else if (alignmentMode
							.equals(MappingMode.USER.toString())) {
						mapping = new UserMapping(userAlignmentCommand);
					} else {
						MappingLogger.getLogger().error(
								"Mapping mode not set. Exit.");
						System.exit(0);
					}
					if (isFirstMapping) {
						int tempMappingQualityFilterGenomic = mappingQualityFilterGenomic;
						if (!isRefine && isTranscriptMapping) {
							tempMappingQualityFilterGenomic = 0;
						}
						mapping.executeMapping(threads, referenceFileName,
								readFileName, outputPrefixBowtieGenomic,
								tempMappingQualityFilterGenomic,
								additionalOptions);
						if (!isRefine && isTranscriptMapping) {
							MappingLogger
									.getLogger()
									.info("Extract unaligned/weakly aligned reads (with MAPQ < "
											+ mappingQualityFilterGenomic
											+ ") from "
											+ "BAM file and create new FASTQ file");
							extract = new ExtractWeakMappingReads();
							extract.extractReads(genomicMappingFileName,
									genomicMappingFileNameNew,
									unalignedGenomicFileName,
									mappingQualityFilterGenomic);

							// rename and remove temporary mapping file; is now
							// filtered!!!
							mapping.renameFile(genomicMappingFileNameNew,
									genomicMappingFileName);
							mapping.removeFile(genomicMappingFileNameNew);
						}
						mapping.sortByCoordinateAndIndex(genomicMappingFileName);
					}

					if (isRefine) {
						MappingLogger
								.getLogger()
								.debug("Reinitialize mapper with PARMA mapping algorithm");
						mapping = new PARMAMapping();
						alignmentMode = MappingMode.PARMA.toString();
						referenceFileName = referenceRefineFileName;
						additionalOptions = parmaMismatches;

						if (isFirstMapping) {
							// only if error profile was not given by the user
							// and a
							// first alignment was performed
							MappingLogger
									.getLogger()
									.debug("Calculate error profile from first mapping for PARMA");
							ErrorProfiling profile = new ErrorProfiling(
									genomicMappingFileName, referenceFileName,
									maxReadLength);
							profile.inferErrorProfile(false, false);
							profileFileName = genomicMappingFileName
									+ ".errorprofile";
							indelFileName = genomicMappingFileName
									+ ".indelprofile";
						}

						((PARMAMapping) mapping)
								.setErrorProfileFilename(profileFileName);
						((PARMAMapping) mapping)
								.setIndelProfileFilename(indelFileName);

						outputPrefixBowtieGenomic = outputPrefix + "."
								+ "PARMA" + "-genomic";
						genomicMappingFileName = outputPrefixBowtieGenomic
								+ ".bam";

						MappingLogger
								.getLogger()
								.debug("Start PARMA algorithm on genomic mapping file "
										+ "using error profile calculated in the step before");

						int tempMappingQualityFilterGenomic = 0;
						if (!isTranscriptMapping) {
							tempMappingQualityFilterGenomic = mappingQualityFilterGenomic;
						}
						mapping.executeMapping(threads, referenceFileName,
								readFileName, outputPrefixBowtieGenomic,
								tempMappingQualityFilterGenomic,
								additionalOptions);
						if (isTranscriptMapping) {
							MappingLogger
									.getLogger()
									.info("Extract unaligned/weakly aligned reads (with MAPQ < "
											+ mappingQualityFilterGenomic
											+ ") from "
											+ "BAM file and create new FASTQ file");
							extract = new ExtractWeakMappingReads();
							extract.extractReads(genomicMappingFileName,
									genomicMappingFileNameNew,
									unalignedGenomicFileName,
									mappingQualityFilterGenomic);

							mapping.renameFile(genomicMappingFileNameNew,
									genomicMappingFileName);
							mapping.removeFile(genomicMappingFileNameNew);
						}
						mapping.sortByCoordinateAndIndex(genomicMappingFileName);

					}

					if (isTranscriptMapping) {
						String outputPrefixTranscript = outputPrefix + "."
								+ alignmentMode + "-transcript";
						String transcriptMappingFileName = outputPrefixTranscript
								+ ".bam";
						mappingFileName = outputPrefix + ".combined.bam";
						mapping.executeMapping(threads, transcriptFileName,
								unalignedGenomicFileName,
								outputPrefixTranscript,
								mappingQualityFilterTranscript,
								additionalOptions);
						mapping.sortByNameAndIndex(transcriptMappingFileName);

						mapping.removeFile(unalignedGenomicFileName);

						// Start postprocessing after transcriptional alignment
						MappingLogger
								.getLogger()
								.info("Combine genomic mapping results with transcriptomic mapping results.");
						CombineGenomeTranscript combiner = new CombineGenomeTranscript();
						combiner.combine(genomicMappingFileName,
								transcriptMappingFileName, mappingFileName);
						mapping.sortByCoordinateAndIndex(mappingFileName);
						MappingLogger.getLogger().info("Combine finished.");
					} else {
						mappingFileName = genomicMappingFileName;
					}
				}
			} catch (NumberFormatException e) {
				MappingLogger.getLogger().error(
						"Abort due to error: parameters read-length "
								+ "or #threads are no valid numbers.");
			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger.getLogger().error(
						"Abort due to error: wrong number of input parameters given. "
								+ "Please consider help for correct usage:");
				e.printStackTrace();
				help.printMappingToolHelp();
			} catch (MappingErrorException e) {
				MappingLogger.getLogger().error(
						"The following mapping command was executed "
								+ "and resulted in an unknown error!");
				MappingLogger.getLogger().error(e.getMappingCommand());
				e.printStackTrace();
			}

			MappingLogger.getLogger()
					.info("Algorithm finished without errors.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.COMBINE)) {
			// ############## COMBINE TOOL ##############
			String genomicMappingFileName = null;
			String transcriptMappingFileName = null;
			String mappingFileName = null;
			if (args.length == 1) {
				help.printCombineToolHelp();
				System.exit(0);
			}
			try {
				genomicMappingFileName = args[1];
				transcriptMappingFileName = args[2];
				mappingFileName = args[3];
			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger
						.getLogger()
						.error("Wrong number of input files. Please consider the following help:");
				help.printCombineToolHelp();
				System.exit(0);
			}

			MappingLogger.getLogger().debug(
					"Start combining genomic and transcriptomic mapping");
			CombineGenomeTranscript combiner = new CombineGenomeTranscript();
			combiner.combine(genomicMappingFileName, transcriptMappingFileName,
					mappingFileName);
			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.BENCHMARK)) {
			// ############## BENCHMARK TOOL ##############

			String mappingFileName = null;
			String outStatistics = null;
			String readsFile = null;

			if (args.length == 1) {
				help.printBenchmarkToolHelp();
				System.exit(0);
			}
			try {
				mappingFileName = args[1];
				outStatistics = args[2];
				readsFile = args[3];
			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger
						.getLogger()
						.error("Wrong number of input files. Please consider the following help:");
				help.printBenchmarkToolHelp();
				System.exit(0);
			}
			MappingLogger.getLogger().debug("Start benchmarking mapping file");
			ValidateBenchmarkStatisticsPARCLIP validate = new ValidateBenchmarkStatisticsPARCLIP();
			validate.calculateBenchmarkStatistics(mappingFileName,
					outStatistics, readsFile);
			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.ERRORPROFILE)) {
			// ############## ERRORPROFILE TOOL ##############

			String mappingFileName = null;
			String referenceFileName = null;
			int maxReadLength = 0;
			boolean isQualityCalc = false;
			boolean showErrorPlot = false;

			if (args.length == 1) {
				help.printErrorprofileToolHelp(isQualityCalc, showErrorPlot);
				System.exit(0);
			}
			try {
				mappingFileName = args[1];
				referenceFileName = args[2];
				maxReadLength = Integer.parseInt(args[3]);
				for (int i = 4; i < args.length; i++) {
					switch (args[i]) {
					case "-q":
						i++;
						isQualityCalc = Boolean.parseBoolean(args[i]);
						break;
					case "-p":
						i++;
						showErrorPlot = Boolean.parseBoolean(args[i]);
						break;
					default:
						MappingLogger.getLogger().error(
								"Parameter \"" + args[i] + "\" not found.");
					}
				}
			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger
						.getLogger()
						.error("Wrong number of input files. Please consider the following help:");
				help.printErrorprofileToolHelp(isQualityCalc, showErrorPlot);
				System.exit(0);
			} catch (NumberFormatException e) {
				MappingLogger.getLogger().error(
						"Wrong input, MAX_READ_LENGTH is not a number.");
				System.exit(0);
			}

			MappingLogger.getLogger().debug("Start errorprofile calculation");
			ErrorProfiling profiling = new ErrorProfiling(mappingFileName,
					referenceFileName, maxReadLength);
			profiling.inferErrorProfile(isQualityCalc, showErrorPlot);
			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.CLUSTERING)) {
			// ############## CLUSTERING TOOL ##############

			String mappingFileName = null;
			String referenceFileName = null;
			String outputFileName = null;
			String snpVCFFile = null;
			int minReadCoverage = 0;

			if (args.length == 1) {
				help.printClusterToolHelp();
				System.exit(0);
			}
			try {
				mappingFileName = args[1];
				referenceFileName = args[2];
				outputFileName = args[3];
				snpVCFFile = args[4];
				minReadCoverage = Integer.parseInt(args[5]);
			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger
						.getLogger()
						.error("Wrong number of input files. Please consider the following help:");
				help.printClusterToolHelp();
				System.exit(0);
			} catch (NumberFormatException e) {
				MappingLogger.getLogger().error(
						"Wrong input, MIN_COVERAGE is not a number.");
				System.exit(0);
			}
			MappingLogger
					.getLogger()
					.debug("Calculate read clusters of overlapping reads for mappingFile");
			PileupClusters clusters = new PileupClusters();
			clusters.calculateReadPileups(mappingFileName, referenceFileName,
					outputFileName, snpVCFFile, minReadCoverage);
			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.SETUP)) {
			// ############## SETUP ##############

			if (args.length == 1) {
				help.printSetupHelp();
				System.exit(0);
			}
			for (int i = 1; i < args.length; i++) {
				switch (args[i]) {
				case "--parma":
					PARMAProperties.setProperty(
							PARMAPropertiesEnum.PARMA_LOCATION, args[i + 1]);
					break;
				case "--bwa":
					PARMAProperties.setProperty(
							PARMAPropertiesEnum.BWA_LOCATION, args[i + 1]);
					break;
				case "--bt2":
					PARMAProperties.setProperty(
							PARMAPropertiesEnum.BT2_LOCATION, args[i + 1]);
					break;
				default:
					break;
				}
			}
		} else if (mode.equals(ProgramMode.SIMULATE)) {
			// ############## PAR-CLIP READ SIMUALTOR ##############

			if (args.length == 1) {
				help.printSimulateHelp();
				System.exit(0);
			}
			for (int i = 1; i <= 8; i++) {
				if (args[i].startsWith("-")) {
					help.printSimulateHelp();
					System.exit(0);
				}
			}
			String mathRandomPath = "";
			String transcriptDatabase;
			String outputPrefix;
			String errorProfile;
			String t2cProfile;
			String t2cPositionProfile;
			String indelProfile;
			String qualityDistribution;
			Double rbpBound;

			for (int i = 9; i < args.length; i++) {
				switch (args[i]) {
				case "-I":
					i++;
					mathRandomPath = args[i];
					break;
				default:
					break;
				}
			}

			try {
				transcriptDatabase = args[1];
				outputPrefix = args[2];
				errorProfile = args[3];
				t2cProfile = args[4];
				t2cPositionProfile = args[5];
				indelProfile = args[6];
				qualityDistribution = args[7];
				rbpBound = Double.parseDouble(args[8]);

				String scriptFileName = "createSimulatedPARCLIPDataset.pl";

				List<String> commands = new LinkedList<String>();
				commands.add("perl");
				if (!mathRandomPath.isEmpty()) {
					commands.add("-I");
					commands.add(mathRandomPath);
				}
				commands.add(scriptFileName);
				commands.add(transcriptDatabase);
				commands.add(outputPrefix);
				commands.add(errorProfile);
				commands.add(t2cProfile);
				commands.add(t2cPositionProfile);
				commands.add(indelProfile);
				commands.add(qualityDistribution);
				commands.add(rbpBound.toString());

				ProcessBuilder mappingProcessBuilder = new ProcessBuilder(
						commands);
				mappingProcessBuilder.redirectOutput(Redirect.INHERIT);
				mappingProcessBuilder.redirectError(Redirect.INHERIT);
				Process mappingProcess = mappingProcessBuilder.start();
				mappingProcess.waitFor();

			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger
						.getLogger()
						.error("Wrong number of input files. Please consider the following help:");
				help.printClusterToolHelp();
				System.exit(0);
			} catch (NumberFormatException e) {
				MappingLogger
						.getLogger()
						.error("Parameter 7 must be a valid double-value for "
								+ "the fraction of RBP-bound clusters to be simualted.");
				help.printClusterToolHelp();
				System.exit(0);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
}
