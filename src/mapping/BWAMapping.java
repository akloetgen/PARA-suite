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
public class BWAMapping extends Mapping {

	public void executeMapping(int threads, String reference, String input,
			String outputPrefix, int mappingQualityFilter,
			String additionalOptions) throws MappingErrorException {
		try {
			// String mismatches = "2";
			setTimeStart();
			MappingLogger.getLogger().info(
					"Starting BWA mapping to investigate error-profile");
			List<String> bwaCommandList = new LinkedList<String>();
			String bwaLocation = PARMAProperties
					.getProperty(PARMAPropertiesEnum.BWA_LOCATION);
			if (bwaLocation == null) {
				bwaLocation = "";
			}
			bwaCommandList.add(bwaLocation + "bwa");
			bwaCommandList.add("aln");
			bwaCommandList.add("-t");
			bwaCommandList.add(threads + "");
			bwaCommandList.add("-n");
			bwaCommandList.add(additionalOptions);
			bwaCommandList.add(reference);
			bwaCommandList.add(input);
			bwaCommandList.add("-f");
			bwaCommandList.add(outputPrefix + ".sai");
			// MappingLogger.getLogger().debug(bwaCommand + bwaOptions);
			if (executeCommand(bwaCommandList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandList);
				throw e;
			}
			MappingLogger.getLogger().info("Convert BWA mapping to SAM file");
			bwaCommandList.clear();
			bwaCommandList.add("bwa");
			bwaCommandList.add("samse");
			bwaCommandList.add(reference);
			bwaCommandList.add(outputPrefix + ".sai");
			bwaCommandList.add(input);
			bwaCommandList.add("-f");
			bwaCommandList.add(outputPrefix + ".sam");
			// MappingLogger.getLogger().debug(bwaCommand);
			if (executeCommand(bwaCommandList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandList);
				throw e;
			}
			MappingLogger.getLogger().info(
					"Time passed for BWA algorithm: " + calculatePassedTime()
							+ " seconds elapsed for BWA alignment");

			MappingLogger.getLogger().info(
					"Convert SAM-file of BWA mapped reads to BAM-file");
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

			// MAYBE CHECK FOR MAPQFILTER EQUALS 0 AND THEN SKIP THE FOLLOWING
			// STEP!!!
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
			cleanUpCommandsList.add(outputPrefix + ".sai");
			executeCommand(cleanUpCommandsList, StreamRedirect.ALL);
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
