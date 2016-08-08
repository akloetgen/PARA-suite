package mapping;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import main.MappingLogger;
import main.PARAsuiteProperties;
import enums.PARAsuitePropertiesEnum;
import enums.StreamRedirect;

/**
 * Read alignment using Bowtie.
 * 
 * @author akloetgen
 * 
 */
public class BowtieMapping extends Mapping {
	public void executeMapping(int threads, String reference, String input,
			String outputPrefix, int mappingQualityFilter,
			String additionalOptions) throws MappingErrorException {
		try {
			List<String> bowtieCommandList = new LinkedList<String>();
			String btLocation = PARAsuiteProperties
					.getProperty(PARAsuitePropertiesEnum.BT_LOCATION);
			if (btLocation == null) {
				btLocation = "";
			} else if (!btLocation.endsWith("/")) {
				btLocation += "/";
			}

			// check whether BWA index exists for reference file, else calculate
			// it
			File indexFile = new File(reference + ".1.ebwt");
			if (!indexFile.exists()) {
				// create index
				bowtieCommandList.add(btLocation + "bowtie-build");
				bowtieCommandList.add(reference);
				bowtieCommandList.add(reference);
				MappingLogger
						.getLogger()
						.info("Creating Bowtie-index for reference file using Bowtie");
				if (executeCommand(bowtieCommandList, StreamRedirect.ERROR) != 0) {
					MappingErrorException e = new MappingErrorException();
					e.setMappingCommand(bowtieCommandList);
					throw e;
				}
				bowtieCommandList.clear();
			}

			setTimeStart();
			MappingLogger.getLogger().info(
					"Starting Bowtie mapping to investigate error-profile");

			bowtieCommandList.add(btLocation + "bowtie");
			// -v for additional options? is mismatch counter...
			bowtieCommandList.add("-S");
			bowtieCommandList.add("-p");
			bowtieCommandList.add(threads + "");
			bowtieCommandList.add("-v");
			bowtieCommandList.add(additionalOptions);
			bowtieCommandList.add("--best");
			bowtieCommandList.add(reference);
			bowtieCommandList.add("-q");
			bowtieCommandList.add(input);
			bowtieCommandList.add(outputPrefix + ".sam");
			if (executeCommand(bowtieCommandList, StreamRedirect.ALL) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bowtieCommandList);
				throw e;
			}

			MappingLogger.getLogger().debug(
					"Time passed for Bowtie algorithm: "
							+ calculatePassedTime()
							+ " seconds elapsed for Bowtie alignment");

			MappingLogger
					.getLogger()
					.info("Convert SAM-file of Bowtie2 mapped reads to BAM-file, filter, sort, index, remove temp files");
			MappingLogger.getLogger().debug(
					"Convert SAM-file of Bowtie mapped reads to BAM-file");
			List<String> cleanUpCommandsList = new LinkedList<String>();
			cleanUpCommandsList.add("samtools");
			cleanUpCommandsList.add("view");
			cleanUpCommandsList.add("-bS");
			cleanUpCommandsList.add("-t");
			cleanUpCommandsList.add(reference);
			cleanUpCommandsList.add(outputPrefix + ".sam");
			cleanUpCommandsList.add("-o");
			cleanUpCommandsList.add(outputPrefix + ".bam");
			executeCommand(cleanUpCommandsList, StreamRedirect.ALL);

			MappingLogger.getLogger().debug(
					"Filtering mapped reads with MAPQ lower than "
							+ mappingQualityFilter);
			cleanUpCommandsList.clear();
			cleanUpCommandsList.add("samtools");
			cleanUpCommandsList.add("view");
			cleanUpCommandsList.add("-q");
			cleanUpCommandsList.add("" + mappingQualityFilter);
			cleanUpCommandsList.add("-b");
			cleanUpCommandsList.add(outputPrefix + ".bam");
			cleanUpCommandsList.add("-o");
			cleanUpCommandsList.add(outputPrefix + ".unique.bam");
			executeCommand(cleanUpCommandsList, StreamRedirect.ALL);

			MappingLogger.getLogger().debug("Removing temporary files");
			cleanUpCommandsList.clear();
			cleanUpCommandsList.add("rm");
			cleanUpCommandsList.add(outputPrefix + ".sam");
			executeCommand(cleanUpCommandsList, StreamRedirect.ALL);
			cleanUpCommandsList.clear();
			cleanUpCommandsList.add("rm");
			cleanUpCommandsList.add(outputPrefix + ".bam");
			executeCommand(cleanUpCommandsList, StreamRedirect.ALL);
			cleanUpCommandsList.clear();
			cleanUpCommandsList.add("mv");
			cleanUpCommandsList.add(outputPrefix + ".unique.bam");
			cleanUpCommandsList.add(outputPrefix + ".bam");
			executeCommand(cleanUpCommandsList, StreamRedirect.ALL);

		} catch (IOException e) {
			MappingLogger.getLogger().error(
					"IO exception thrown in FirstMapping with "
							+ "following error message:"
							+ System.getProperty("line.separator")
							+ e.getMessage());
			e.printStackTrace();
		} catch (InterruptedException e) {
			MappingLogger.getLogger().error(
					"Interrupted exception thrown in FirstMapping "
							+ "with following error message:"
							+ System.getProperty("line.separator")
							+ e.getMessage());
			e.printStackTrace();
		}
	}
}
