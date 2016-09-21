package mapping;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import main.MappingLogger;
import main.PARAsuiteProperties;
import enums.PARAsuitePropertiesEnum;
import enums.StreamRedirect;

/**
 * Read alignment using Bowtie2.
 * 
 * @author akloetgen
 * 
 */
public class Bowtie2Mapping extends Mapping {

	public void executeMapping(int threads, String reference, String input,
			String outputPrefix, int mappingQualityFilter,
			String additionalOptions) {

		List<String> bowtieCommandList = new LinkedList<String>();
		String bt2Location = PARAsuiteProperties
				.getProperty(PARAsuitePropertiesEnum.BT2_LOCATION);
		if (bt2Location == null) {
			bt2Location = "";
		} else if (!bt2Location.endsWith("/")) {
			bt2Location += "/";
		}

		// check whether BWA index exists for reference file, else calculate
		// it with the PARMA-extension of BWA
		File indexFile = new File(reference + ".bwt");
		if (!indexFile.exists()) {
			// create index
			bowtieCommandList.add(bt2Location + "bowtie2-build");
			bowtieCommandList.add(reference);
			bowtieCommandList.add(reference);
			MappingLogger.getLogger().info(
					"Creating BWA-index for reference file using Bowtie2");
			executeCommand(bowtieCommandList, StreamRedirect.ERROR);
			bowtieCommandList.clear();
		}

		setTimeStart();
		MappingLogger.getLogger().info(
				"Starting Bowtie2 mapping to investigate error-profile");

		bowtieCommandList.add(bt2Location + "bowtie2");
		bowtieCommandList.add("-p " + threads + " -x " + reference + " -U "
				+ input + " -S " + outputPrefix + ".sam" + additionalOptions);
		executeCommand(bowtieCommandList, StreamRedirect.ALL);

		MappingLogger.getLogger().debug(
				"Time passed for Bowtie2 algorithm: " + calculatePassedTime()
						+ " seconds elapsed for Bowtie2 alignment");

		MappingLogger
				.getLogger()
				.info("Convert SAM-file of Bowtie2 mapped reads to BAM-file, filter, sort, index, remove temp files");
		MappingLogger.getLogger().debug(
				"Convert SAM-file of Bowtie2 mapped reads to BAM-file");
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

	}
}
