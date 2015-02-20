package mapping;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import enums.PARMAPropertiesEnum;
import enums.StreamRedirect;
import main.MappingLogger;
import main.PARMAProperties;

/**
 * 
 * @author akloetgen
 * 
 */
public class BowtieMapping extends Mapping {

	public void executeMapping(int threads, String reference, String input,
			String outputPrefix, int mappingQualityFilter,
			String additionalOptions) throws MappingErrorException {
		try {
			setTimeStart();
			MappingLogger.getLogger().info(
					"Starting Bowtie2 mapping to investigate error-profile");
			List<String> bowtieCommandList = new LinkedList<String>();
			String bt2Location = PARMAProperties
					.getProperty(PARMAPropertiesEnum.BT2_LOCATION);
			if (bt2Location == null) {
				bt2Location = "";
			}
			bowtieCommandList.add(bt2Location + "bowtie2");
			bowtieCommandList.add("-p " + threads + " -x " + reference + " -U "
					+ input + " -S " + outputPrefix + ".sam"
					+ additionalOptions);
			// MappingLogger.getLogger().debug(bowtie2Command);
			if (executeCommand(bowtieCommandList, StreamRedirect.ALL) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bowtieCommandList);
				throw e;
			}

			MappingLogger.getLogger().info(
					"Time passed for Bowtie2 algorithm: "
							+ calculatePassedTime()
							+ " seconds elapsed for Bowtie2 alignment");

			MappingLogger.getLogger().info(
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

			MappingLogger.getLogger().info(
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

			MappingLogger.getLogger().info("Removing temporary files");
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

			// sort and index
			// sortByCoordinateAndIndex(outputPrefixBowtie + ".bam");

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