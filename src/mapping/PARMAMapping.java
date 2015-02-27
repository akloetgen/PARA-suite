package mapping;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import main.MappingLogger;
import main.PARMAProperties;
import enums.PARMAPropertiesEnum;
import enums.StreamRedirect;

/**
 * 
 * @author akloetgen
 * 
 */
public class PARMAMapping extends Mapping {

	private String errorProfileFilename;
	private String indelProfileFilename;

	public void setErrorProfileFilename(String errorProfileFilename) {
		this.errorProfileFilename = errorProfileFilename;
	}

	public void setIndelProfileFilename(String indelProfileFilename) {
		this.indelProfileFilename = indelProfileFilename;
	}

	public void executeMapping(int threads, String reference, String input,
			String outputPrefix, int mappingQualityFilter,
			String additionalOptions) throws MappingErrorException {
		try {
			List<String> bwaCommandsList = new LinkedList<String>();
			String parmaLocation = PARMAProperties
					.getProperty(PARMAPropertiesEnum.PARMA_LOCATION);
			if (parmaLocation == null) {
				parmaLocation = "";
			} else if (!parmaLocation.endsWith("/")) {
				parmaLocation += "/";
			}

			// check whether BWA index exists for reference file, else calculate
			// it with the PARMA-extension of BWA
			File indexFile = new File(reference + ".bwt");
			if (!indexFile.exists()) {
				// create index
				bwaCommandsList.add(parmaLocation + "bwa");
				bwaCommandsList.add("index");
				bwaCommandsList.add(reference);
				MappingLogger
						.getLogger()
						.info("Creating BWA-index for reference file using BWA 0.7.8");
				if (executeCommand(bwaCommandsList, StreamRedirect.ERROR) != 0) {
					MappingErrorException e = new MappingErrorException();
					e.setMappingCommand(bwaCommandsList);
					throw e;
				}
				bwaCommandsList.clear();
			}

			setTimeStart();

			// String ep_filename =
			// "/home/akloetgen/read_mapper/bwa-0.7.8_ep/test/error_profile_MSI1+MSI2_bowtie2.tsv";
			MappingLogger.getLogger().info(
					"Starting PARMA mapping to investigate error-profile");

			bwaCommandsList.add(parmaLocation + "bwa");
			bwaCommandsList.add("parma");
			bwaCommandsList.add("-t");
			bwaCommandsList.add(threads + "");
			bwaCommandsList.add("-n");
			bwaCommandsList.add(additionalOptions);
			bwaCommandsList.add("-p");
			bwaCommandsList.add(errorProfileFilename);
			bwaCommandsList.add("-g");
			bwaCommandsList.add(indelProfileFilename);
			bwaCommandsList.add(reference);
			bwaCommandsList.add(input);
			bwaCommandsList.add("-f");
			bwaCommandsList.add(outputPrefix + ".sai");
			// MappingLogger.getLogger().debug(bwaCommand);
			if (executeCommand(bwaCommandsList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandsList);
				throw e;
			}

			MappingLogger
					.getLogger()
					.info("Convert SAM-file of PARMA mapped reads to BAM-file, filter, sort, index, remove temp files");
			MappingLogger.getLogger()
					.debug("Convert PARMA mapping to SAM file");
			bwaCommandsList.clear();
			bwaCommandsList.add(parmaLocation + "bwa");
			bwaCommandsList.add("samse");
			bwaCommandsList.add(reference);
			bwaCommandsList.add(outputPrefix + ".sai");
			bwaCommandsList.add(input);
			bwaCommandsList.add("-f");
			bwaCommandsList.add(outputPrefix + ".sam");
			// MappingLogger.getLogger().debug(bwaCommand);
			if (executeCommand(bwaCommandsList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandsList);
				throw e;
			}

			MappingLogger.getLogger().debug(
					"Time passed for PARMA algorithm: " + calculatePassedTime()
							+ " seconds elapsed for PARMA alignment");

			MappingLogger.getLogger().debug(
					"Convert SAM-file of PARMA mapped reads to BAM-file");
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
					"IO exception thrown in PARMA with "
							+ "following error message:"
							+ System.getProperty("line.separator")
							+ e.getMessage());
			e.printStackTrace();
		} catch (InterruptedException e) {
			MappingLogger.getLogger().error(
					"Interrupted exception thrown in PARMA "
							+ "with following error message:"
							+ System.getProperty("line.separator")
							+ e.getMessage());
			e.printStackTrace();
		}
	}
}
