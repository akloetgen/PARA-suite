package mapping;

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

	// private Properties properties;

	public PARMAMapping() {
		// try {
		// File newPropertiesFile = new File("parma.properties");
		// if (!newPropertiesFile.exists()) {
		// newPropertiesFile.createNewFile();
		// }
		// properties = new Properties();
		// BufferedInputStream stream;
		// stream = new BufferedInputStream(new FileInputStream(
		// newPropertiesFile));
		// properties.load(stream);
		// stream.close();
		// } catch (IOException e) {
		// // TODO Auto-generated catch block
		// e.printStackTrace();
		// }
	}

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
			setTimeStart();

			// String ep_filename =
			// "/home/akloetgen/read_mapper/bwa-0.7.8_ep/test/error_profile_MSI1+MSI2_bowtie2.tsv";
			MappingLogger.getLogger().info(
					"Starting PARMA mapping to investigate error-profile");
			List<String> bwaCommandsList = new LinkedList<String>();
			String parmaLocation = PARMAProperties
					.getProperty(PARMAPropertiesEnum.PARMA_LOCATION);
			if (parmaLocation == null) {
				parmaLocation = "";
			}
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
			MappingLogger.getLogger().info("Convert PARMA mapping to SAM file");
			bwaCommandsList.clear();
			bwaCommandsList.add("/home/akloetgen/read_mapper/bwa-0.7.8_ep/bwa");
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

			MappingLogger.getLogger().info(
					"Time passed for PARMA algorithm: " + calculatePassedTime()
							+ " seconds elapsed for PARMA alignment");

			MappingLogger.getLogger().info(
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
