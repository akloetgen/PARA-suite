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
 * Read alignment using BWA.
 * 
 * @author akloetgen
 * 
 */
public class BWAMapping extends Mapping {

	public void executeMapping(int threads, String reference, String input,
			String outputPrefix, int mappingQualityFilter,
			String additionalOptions) throws MappingErrorException {
		try {
			List<String> bwaCommandList = new LinkedList<String>();
			String bwaLocation = PARAsuiteProperties
					.getProperty(PARAsuitePropertiesEnum.BWA_LOCATION);
			if (bwaLocation == null) {
				bwaLocation = "";
			} else if (!bwaLocation.endsWith("/")) {
				bwaLocation += "/";
			}

			// check whether BWA index exists for reference file, else calculate
			// it with the PARMA-extension of BWA
			File indexFile = new File(reference + ".bwt");
			if (!indexFile.exists()) {
				// create index
				bwaCommandList.add(bwaLocation + "bwa");
				bwaCommandList.add("index");
				bwaCommandList.add(reference);
				MappingLogger.getLogger().info(
						"Creating BWA-index for reference file using BWA");
				if (executeCommand(bwaCommandList, StreamRedirect.ERROR) != 0) {
					MappingErrorException e = new MappingErrorException();
					e.setMappingCommand(bwaCommandList);
					throw e;
				}
				bwaCommandList.clear();
			}

			setTimeStart();
			MappingLogger.getLogger().info(
					"Starting BWA mapping to investigate error-profile");

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
			if (executeCommand(bwaCommandList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandList);
				throw e;
			}

			MappingLogger
					.getLogger()
					.info("Convert SAM-file of Bowtie2 mapped reads to BAM-file, filter, sort, index, remove temp files");
			MappingLogger.getLogger().debug("Convert BWA mapping to SAM file");
			bwaCommandList.clear();
			bwaCommandList.add(bwaLocation + "bwa");
			bwaCommandList.add("samse");
			bwaCommandList.add(reference);
			bwaCommandList.add(outputPrefix + ".sai");
			bwaCommandList.add(input);
			bwaCommandList.add("-f");
			bwaCommandList.add(outputPrefix + ".sam");
			if (executeCommand(bwaCommandList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandList);
				throw e;
			}
			MappingLogger.getLogger().debug(
					"Time passed for BWA algorithm: " + calculatePassedTime()
							+ " seconds elapsed for BWA alignment");

			MappingLogger.getLogger().debug(
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
