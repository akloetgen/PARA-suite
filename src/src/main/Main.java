package main;

import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.LinkedList;
import java.util.List;

import mapping.BWAMapping;
import mapping.Bowtie2Mapping;
import mapping.BowtieMapping;
import mapping.ExternalCallErrorException;
import mapping.Mapping;
import mapping.PARAsuiteMapping;
import mapping.UserMapping;
import utils.benchmarking.ValidateBenchmarkStatisticsClusters;
import utils.benchmarking.ValidateBenchmarkStatisticsPARCLIP;
import utils.errorprofile.ErrorProfiling;
import utils.pileupclusters.FetchSequencesForBEDFile;
import utils.pileupclusters.FetchSequencesForBindingSites;
import utils.pileupclusters.PileupClusters;
import utils.postprocessing.CombineGenomeTranscript;
import utils.postprocessing.ExtractNotMappedReads;
import utils.postprocessing.ExtractWeakMappingReads;
import enums.MappingMode;
import enums.PARAsuitePropertiesEnum;
import enums.ProgramMode;

/**
 * 
 * Main class as entry point for PARA-suite pipeline. Can execute several
 * sub-utils of the PARA-suite project.
 * 
 * @author akloetgen
 * 
 */
public class Main {

	private static final String VERSION = "0.65 beta";
	private static final String AUTHOR = "Andreas Kloetgen";

	public static void main(String[] args) {
		Help help = new Help();
		ProgramMode mode = null;

		if (args.length == 0 || args[0].equals("-h")
				|| args[0].equals("--help")) {
			help.printProgramInfo(ProgramMode.TOOLKIT, VERSION, AUTHOR);
			help.printProgramModeHelp();
			System.exit(1);
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
		case "extract":
			mode = ProgramMode.EXTRACT;
			break;
		case "extractAligned":
			mode = ProgramMode.EXTRACTALIGNED;
			break;
		case "benchmarkClusters":
			mode = ProgramMode.BENCHMARK_CLUSTERS;
			break;
		case "fetch":
			mode = ProgramMode.FETCHSEQUENCES;
			break;
		case "fetchBed":
			mode = ProgramMode.FETCHSEQUENCESBEDFILE;
			break;
		default:
			MappingLogger.getLogger().error(
					"Programmode \"" + args[0]
							+ "\" not found. Please consider the following "
							+ "help-message for a correct usage.");
			help.printProgramModeHelp();
			System.exit(1);
		}
		help.printProgramInfo(mode, VERSION, AUTHOR);

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
				boolean isRefine = false;
				boolean isFirstMapping = true;
				boolean keepUnaligned = false;
				String additionalOptions = "";
				String parasuiteMismatches = "-1";
				String btMismatches = "1";
				String bwaMismatches = "2";

				if (args.length == 1) {
					help.printMappingToolHelp();
					System.exit(1);
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
											+ "one of the following: bt2, bwa, parasuite, user");
							System.exit(1);
						}
						break;
					case "--parasuite-mm":
						i++;
						parasuiteMismatches = args[i];
						break;
					case "--parasuite-ep":
						// COULD BE DELETED AS THIS IS A TK AND NOT A PIPELINE
						i++;
						isFirstMapping = false;
						profileFileName = args[i];
						break;
					case "--parasuite-indel":
						// COULD BE DELETED AS THIS IS A TK AND NOT A PIPELINE
						i++;
						indelFileName = args[i];
						break;
					case "--bwa-mm":
						i++;
						bwaMismatches = args[i];
						break;
					case "--bt-mm":
						i++;
						btMismatches = args[i];
						break;
					case "-c":
						i++;
						userAlignmentCommand = args[i];
						break;
					case "--refine":
						isRefine = true;
						break;
					case "--ref-refine":
						// if the first mapping is not made by BWA, a second
						// reference file for the bwa parasuite mod is
						// necessary.
						i++;
						referenceRefineFileName = args[i];
						break;
					case "--unaligned":
						// keep unaligned reads in a separate file
						keepUnaligned = true;
						break;
					default:
						MappingLogger
								.getLogger()
								.error("Parameter not found: "
										+ args[i]
										+ ". Please consider the following help:");
						help.printMappingToolHelp();
						System.exit(1);
					}
				}
				if (readFileName == null || outputPrefix == null
						|| referenceFileName == null) {
					MappingLogger.getLogger().error(
							"No valid argument for -q, -r or -o given in input. Please consider"
									+ " help for correct usage:");
					help.printMappingToolHelp();
					System.exit(1);
				}

				if (isMapping) {
					String unalignedGenomicFileName = outputPrefix
							+ ".unaligned.fastq";
					String outputPrefixMapperGenomic = outputPrefix + "."
							+ alignmentMode + "-genomic";
					String genomicMappingFileName = outputPrefixMapperGenomic
							+ ".bam";
					String genomicMappingFileNameNew = genomicMappingFileName
							+ ".new";
					Mapping mapping = null;
					ExtractWeakMappingReads extract;
					extract = new ExtractWeakMappingReads();

					if (alignmentMode.equals(MappingMode.BT2.toString())) {
						mapping = new Bowtie2Mapping();
					} else if (alignmentMode.equals(MappingMode.BWA.toString())) {
						mapping = new BWAMapping();
						additionalOptions = bwaMismatches;
					} else if (alignmentMode.equals(MappingMode.PARAsuite
							.toString())) {
						mapping = new PARAsuiteMapping();
						additionalOptions = parasuiteMismatches;
					} else if (alignmentMode
							.equals(MappingMode.USER.toString())) {
						mapping = new UserMapping(userAlignmentCommand);
					} else if (alignmentMode.equals(MappingMode.BT.toString())) {
						mapping = new BowtieMapping();
						additionalOptions = btMismatches;
					} else {
						MappingLogger
								.getLogger()
								.error("Mapping mode not set. Program integrity not given.");
						System.exit(-1);
					}
					if (isFirstMapping) {
						int tempMappingQualityFilterGenomic = mappingQualityFilterGenomic;
						if (!isRefine && isTranscriptMapping) {
							tempMappingQualityFilterGenomic = 0;
						}
						mapping.executeMapping(threads, referenceFileName,
								readFileName, outputPrefixMapperGenomic,
								tempMappingQualityFilterGenomic,
								additionalOptions);
						if (!isRefine && isTranscriptMapping) {
							MappingLogger
									.getLogger()
									.debug("Extract unaligned/weakly aligned reads (with MAPQ < "
											+ mappingQualityFilterGenomic
											+ ") from "
											+ "BAM file and create new FASTQ file");
							extract.extractReads(genomicMappingFileName,
									genomicMappingFileNameNew,
									unalignedGenomicFileName,
									mappingQualityFilterGenomic);

							// rename temporary mapping file
							mapping.renameFile(genomicMappingFileNameNew,
									genomicMappingFileName);
							// mapping.removeFile(genomicMappingFileNameNew);
						}
						mapping.sortByCoordinateAndIndex(genomicMappingFileName);
					}

					if (isRefine) {
						MappingLogger
								.getLogger()
								.debug("Reinitialize mapper with PARA-suite mapping algorithm");
						mapping = new PARAsuiteMapping();
						alignmentMode = MappingMode.PARAsuite.toString();
						referenceFileName = referenceRefineFileName;
						additionalOptions = parasuiteMismatches;

						if (isFirstMapping) {
							// only if error profile was not given by the user
							// and a first alignment was performed
							MappingLogger
									.getLogger()
									.debug("Calculate error profile from first mapping for PARA-suite");
							ErrorProfiling profile = new ErrorProfiling(
									genomicMappingFileName, referenceFileName,
									maxReadLength);
							profile.inferErrorProfile(false, false);
							profileFileName = genomicMappingFileName
									+ ".errorprofile";
							indelFileName = genomicMappingFileName
									+ ".indelprofile";
						}

						((PARAsuiteMapping) mapping)
								.setErrorProfileFilename(profileFileName);
						((PARAsuiteMapping) mapping)
								.setIndelProfileFilename(indelFileName);

						outputPrefixMapperGenomic = outputPrefix + "."
								+ "PARAsuite" + "-genomic";
						genomicMappingFileName = outputPrefixMapperGenomic
								+ ".bam";
						genomicMappingFileNameNew = genomicMappingFileName
								+ ".new";

						MappingLogger
								.getLogger()
								.debug("Start PARA-suite algorithm on genomic mapping file "
										+ "using error profile calculated in the step before");

						int tempMappingQualityFilterGenomic = 0;
						if (!isTranscriptMapping && !keepUnaligned) {
							tempMappingQualityFilterGenomic = mappingQualityFilterGenomic;
						}

						mapping.executeMapping(threads, referenceFileName,
								readFileName, outputPrefixMapperGenomic,
								tempMappingQualityFilterGenomic,
								additionalOptions);
						if (isTranscriptMapping) {
							MappingLogger
									.getLogger()
									.debug("Extract unaligned/weakly aligned reads (with MAPQ < "
											+ mappingQualityFilterGenomic
											+ ") from "
											+ "BAM file and create new FASTQ file");

							extract.extractReads(genomicMappingFileName,
									genomicMappingFileNameNew,
									unalignedGenomicFileName,
									mappingQualityFilterGenomic);

							mapping.renameFile(genomicMappingFileNameNew,
									genomicMappingFileName);
							// mapping.removeFile(genomicMappingFileNameNew);
						}
						mapping.sortByCoordinateAndIndex(genomicMappingFileName);
					}
					if (keepUnaligned) {
						extract.extractReads(genomicMappingFileName,
								genomicMappingFileNameNew,
								unalignedGenomicFileName,
								mappingQualityFilterGenomic);
						mapping.renameFile(genomicMappingFileNameNew,
								genomicMappingFileName);
						// mapping.removeFile(genomicMappingFileNameNew);
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

						if (!keepUnaligned) {
							mapping.removeFile(unalignedGenomicFileName);
						}
						// Start postprocessing after transcriptional alignment
						MappingLogger
								.getLogger()
								.info("Combine genomic mapping results with transcriptomic mapping results.");
						CombineGenomeTranscript combiner = new CombineGenomeTranscript();
						combiner.combine(genomicMappingFileName,
								transcriptMappingFileName, mappingFileName);
						mapping.sortByCoordinateAndIndex(mappingFileName);
						MappingLogger.getLogger().debug("Combine finished.");
					} else {
						mappingFileName = genomicMappingFileName;
					}
				}
			} catch (NumberFormatException e) {
				MappingLogger.getLogger().error(
						"Abort due to error: parameters read-length "
								+ "or #threads are no valid numbers.");
				System.exit(1);
			} catch (ArrayIndexOutOfBoundsException e) {
				e.printStackTrace();
				MappingLogger.getLogger().error(
						"Abort due to error: wrong number of input parameters given. "
								+ "Please consider help for correct usage:");
				help.printMappingToolHelp();
				System.exit(1);
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
				System.exit(1);
			}
			try {

				for (int i = 1; i < args.length; i++) {
					switch (args[i]) {
					case "-g":
						i++;
						genomicMappingFileName = args[i];
						break;
					case "-t":
						i++;
						transcriptMappingFileName = args[i];
						break;
					case "-o":
						i++;
						mappingFileName = args[i];
						break;
					default:
						MappingLogger
								.getLogger()
								.error("Input " + args[i]
										+ " invalid. Please consider the help:");
						help.printCombineToolHelp();
						System.exit(1);
					}
				}

			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger
						.getLogger()
						.error("Wrong number of input files. Please consider the following help:");
				help.printCombineToolHelp();
				System.exit(1);
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
			boolean onlyBoundClusters = false;

			if (args.length == 1) {
				help.printBenchmarkToolHelp();
				System.exit(1);
			}
			try {
				mappingFileName = args[1];
				outStatistics = args[2];
				readsFile = args[3];
				if (args.length >= 5 && args[4].equals("--only-bound")) {
					onlyBoundClusters = true;
				}
			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger.getLogger().error(
						"Wrong number of input files. "
								+ "Please consider the following help:");
				help.printBenchmarkToolHelp();
				System.exit(1);
			}
			MappingLogger.getLogger().debug("Start benchmarking mapping file");
			ValidateBenchmarkStatisticsPARCLIP validate = new ValidateBenchmarkStatisticsPARCLIP();
			validate.calculateBenchmarkStatistics(mappingFileName,
					outStatistics, readsFile, onlyBoundClusters);

			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.BENCHMARK_CLUSTERS)) {
			// ############## BENCHMARK TOOL ##############

			String clustersFileName = null;
			String outStatistics = null;
			String clustersReferenceFile = null;

			if (args.length == 1) {
				help.printBenchmarkToolHelp();
				System.exit(1);
			}
			try {
				clustersFileName = args[1];
				outStatistics = args[2];
				clustersReferenceFile = args[3];
			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger.getLogger().error(
						"Wrong number of input files. "
								+ "Please consider the following help:");
				help.printBenchmarkToolHelp();
				System.exit(1);
			}
			MappingLogger.getLogger().debug("Start benchmarking mapping file");
			ValidateBenchmarkStatisticsClusters validate = new ValidateBenchmarkStatisticsClusters();
			validate.calculateBenchmarkStatistics(clustersFileName,
					outStatistics, clustersReferenceFile);

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
				System.exit(1);
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
				MappingLogger.getLogger().error(
						"Wrong number of input files. "
								+ "Please consider the following help:");
				help.printErrorprofileToolHelp(isQualityCalc, showErrorPlot);
				System.exit(1);
			} catch (NumberFormatException e) {
				MappingLogger.getLogger().error(
						"Wrong input, MAX_READ_LENGTH is not a number.");
				System.exit(1);
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
				System.exit(1);
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
				System.exit(1);
			} catch (NumberFormatException e) {
				MappingLogger.getLogger().error(
						"Wrong input, MIN_COVERAGE is not a number.");
				System.exit(1);
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
				System.exit(1);
			}
			for (int i = 1; i < args.length; i++) {
				switch (args[i]) {
				case "--parasuite":
					PARAsuiteProperties.setProperty(
							PARAsuitePropertiesEnum.PARAsuite_LOCATION,
							args[i + 1]);
					break;
				case "--bwa":
					PARAsuiteProperties.setProperty(
							PARAsuitePropertiesEnum.BWA_LOCATION, args[i + 1]);
					break;
				case "--bt2":
					PARAsuiteProperties.setProperty(
							PARAsuitePropertiesEnum.BT2_LOCATION, args[i + 1]);
					break;
				case "--bowtie":
					PARAsuiteProperties.setProperty(
							PARAsuitePropertiesEnum.BT_LOCATION, args[i + 1]);
					break;
				default:
					MappingLogger
							.getLogger()
							.error("Aligner \""
									+ args[0]
									+ "\" not found. Please consider the following "
									+ "help-message for setting up the PARA-suite properlys.");
					help.printSetupHelp();
					System.exit(1);
				}
			}

			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.SIMULATE)) {
			// ############## PAR-CLIP READ SIMUALTOR ##############

			if (args.length == 1) {
				help.printSimulateHelp();
				System.exit(1);
			}
			for (int i = 1; i <= 8; i++) {
				if (args[i].startsWith("-")) {
					help.printSimulateHelp();
					System.exit(1);
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
				int simulatorExitStatus = mappingProcess.waitFor();

				if (simulatorExitStatus != 0) {
					String commandString = "";
					for (int i = 0; i < commands.size(); i++) {
						commandString += commands.get(i) + " ";
					}
					MappingLogger.getLogger().debug(commandString);
					throw (new ExternalCallErrorException(commandString));
				}

			} catch (ArrayIndexOutOfBoundsException e) {
				MappingLogger
						.getLogger()
						.error("Wrong number of input files. Please consider the following help:");
				help.printClusterToolHelp();
				System.exit(1);
			} catch (NumberFormatException e) {
				MappingLogger
						.getLogger()
						.error("Parameter 7 must be a valid double-value for "
								+ "the fraction of RBP-bound clusters to be simualted.");
				help.printClusterToolHelp();
				System.exit(1);
			} catch (IOException e) {
				MappingLogger
						.getLogger()
						.error("File not found. Please consider that all "
								+ "input files exist and that the tool has "
								+ "permissions to write in the output directory.");
				System.exit(1);
			} catch (InterruptedException e) {
				MappingLogger
						.getLogger()
						.error("Simulator failed to return. Please consider "
								+ "the .log file of the PARA-suite's simulator. "
								+ "The simulator has to run properly before the "
								+ "pipeline can continue.");
				System.exit(1);
			} catch (ExternalCallErrorException e) {
				MappingLogger
						.getLogger()
						.error("Simulator exited with non-zero exit status. Please consider "
								+ "the .log file of the PARA-suite's simulator. "
								+ "The simulator has to run properly before the "
								+ "pipeline can continue.");
				MappingLogger.getLogger().error(e.getMappingCommand());
				System.exit(1);
			}

			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.EXTRACT)) {
			// ############## EXTRACTING WEAK MAPPINGS FROM BAM ##############
			int mapqThreshold = 10;
			if (args.length == 1) {
				help.printExtractHelp(mapqThreshold);
				System.exit(1);
			}
			String inputBAMFile = null;
			String outputFastqFile = null;
			String tempBAMFile = "tempfile.bam";

			for (int i = 1; i < args.length; i++) {
				switch (args[i]) {
				case "-i":
					inputBAMFile = args[i + 1];
					break;
				case "-o":
					outputFastqFile = args[i + 1];
					break;
				case "-t":
					mapqThreshold = Integer.parseInt(args[i + 1]);
					break;
				default:
					MappingLogger
							.getLogger()
							.error("Input \""
									+ args[0]
									+ "\" invalid. Please consider the following "
									+ "help-message.");
					help.printExtractHelp(mapqThreshold);
					System.exit(1);
				}
			}
			if (inputBAMFile == null || outputFastqFile == null) {
				help.printExtractHelp(mapqThreshold);
				System.exit(1);
			}

			ExtractWeakMappingReads extract = new ExtractWeakMappingReads();
			extract.extractReads(inputBAMFile, tempBAMFile, outputFastqFile,
					mapqThreshold);
			Mapping mapping = new BWAMapping();
			// rename and remove temporary mapping file; is now
			// filtered!!!
			mapping.renameFile(tempBAMFile, inputBAMFile);
			mapping.removeFile(tempBAMFile);

			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.EXTRACTALIGNED)) {
			// ####### EXTRACTING ALL ALIGNED READ NAMES FROM BAM #######
			String inputBAMFile = null;
			String outputReads = null;
			for (int i = 1; i < args.length; i++) {
				switch (args[i]) {
				case "-a":
					inputBAMFile = args[i + 1];
					break;
				case "-o":
					outputReads = args[i + 1];
					break;

				default:
					MappingLogger.getLogger().error(
							"Input \"" + args[0]
									+ "\" invalid. -a for valid BAM input "
									+ "and -o for writeable output.");
					System.exit(1);
				}
			}
			if (inputBAMFile == null || outputReads == null) {
				MappingLogger.getLogger().error("error message missing!!!");
				System.exit(1);
			}

			ExtractNotMappedReads extract = new ExtractNotMappedReads();
			extract.extractReads(inputBAMFile, outputReads);

			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.FETCHSEQUENCES)) {
			// ####### FETCHING SEQUENCES FOR CLUSTERS #######
			String inputBindingSites = null;
			String referenceFile = null;
			String outputFile = null;
			for (int i = 1; i < args.length; i++) {
				switch (args[i]) {
				case "-i":
					inputBindingSites = args[i + 1];
					break;
				case "-r":
					referenceFile = args[i + 1];
					break;
				case "-o":
					outputFile = args[i + 1];
					break;
				}
			}
			if (inputBindingSites == null || referenceFile == null
					|| outputFile == null) {
				// help.printExtractHelp(mapqThreshold);
				MappingLogger.getLogger().error(
						"-i inputBindingSites; -r reference; -o outputFile");
				System.exit(1);
			}

			FetchSequencesForBindingSites fetch = new FetchSequencesForBindingSites();
			fetch.fetchSequences(referenceFile, inputBindingSites, outputFile);

			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		} else if (mode.equals(ProgramMode.FETCHSEQUENCESBEDFILE)) {
			// ####### FETCHING SEQUENCES FOR CLUSTERS #######
			String inputBindingSites = null;
			String referenceFile = null;
			String outputFile = null;
			for (int i = 1; i < args.length; i++) {
				switch (args[i]) {
				case "-i":
					inputBindingSites = args[i + 1];
					break;
				case "-r":
					referenceFile = args[i + 1];
					break;
				case "-o":
					outputFile = args[i + 1];
					break;
				}
			}
			if (inputBindingSites == null || referenceFile == null
					|| outputFile == null) {
				// help.printExtractHelp(mapqThreshold);
				MappingLogger.getLogger().error(
						"-i inputBindingSites; -r reference; -o outputFile");
				System.exit(1);
			}

			FetchSequencesForBEDFile fetch = new FetchSequencesForBEDFile();
			fetch.fetchSequences(referenceFile, inputBindingSites, outputFile);

			MappingLogger.getLogger().debug("Finished.");
			System.exit(0);
		}

	}
}
