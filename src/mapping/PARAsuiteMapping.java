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
 * 
 * @author akloetgen
 * 
 */
public class PARAsuiteMapping extends Mapping {

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
			String parasuiteLocation = PARAsuiteProperties
					.getProperty(PARAsuitePropertiesEnum.PARAsuite_LOCATION);
			if (parasuiteLocation == null) {
				parasuiteLocation = "";
			} else if (!parasuiteLocation.endsWith("/")) {
				parasuiteLocation += "/";
			}

			// check whether BWA index exists for reference file, else calculate
			// it with the PARA-suite-extension of BWA
			File indexFile = new File(reference + ".bwt");
			if (!indexFile.exists()) {
				// create index
				bwaCommandsList.add(parasuiteLocation + "bwa");
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

			MappingLogger.getLogger().info(
					"Starting PARA-suite mapping to investigate error-profile");

			//TODO rename the BWA extension!!!
			bwaCommandsList.add(parasuiteLocation + "bwa");
			bwaCommandsList.add("parasuite");
			bwaCommandsList.add("-t");
			bwaCommandsList.add(threads + "");
			bwaCommandsList.add("-X");
			bwaCommandsList.add(additionalOptions);
			bwaCommandsList.add("-p");
			bwaCommandsList.add(errorProfileFilename);
			bwaCommandsList.add("-g");
			bwaCommandsList.add(indelProfileFilename);
			bwaCommandsList.add(reference);
			bwaCommandsList.add(input);
			bwaCommandsList.add("-f");
			bwaCommandsList.add(outputPrefix + ".sai");
			if (executeCommand(bwaCommandsList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandsList);
				throw e;
			}

			MappingLogger
					.getLogger()
					.info("Convert SAM-file of PARA-suite mapped reads to BAM-file, filter, sort, index, remove temp files");
			MappingLogger.getLogger()
					.debug("Convert PARA-suite mapping to SAM file");
			bwaCommandsList.clear();
			bwaCommandsList.add(parasuiteLocation + "bwa");
			bwaCommandsList.add("samse");
			bwaCommandsList.add(reference);
			bwaCommandsList.add(outputPrefix + ".sai");
			bwaCommandsList.add(input);
			bwaCommandsList.add("-f");
			bwaCommandsList.add(outputPrefix + ".sam");
			if (executeCommand(bwaCommandsList, StreamRedirect.ERROR) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(bwaCommandsList);
				throw e;
			}

			MappingLogger.getLogger().debug(
					"Time passed for PARA-suite algorithm: " + calculatePassedTime()
							+ " seconds elapsed for PARMA alignment");

			// folgendes verschieben!!!!!
			MappingLogger.getLogger().debug(
					"Convert SAM-file of PARA-suite mapped reads to BAM-file");
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

			// vorübergehende lösung. die samtools commands müssen in Mapping
			// aufgenommen und aus Main ausgeführt werden. das extract sollte
			// vorher immer erfolgen, damit man die FASTQ weiter nutzen kann
//			ExtractWeakMappingReads extract = new ExtractWeakMappingReads();
//			extract.extractReads(outputPrefix + ".bam", outputPrefix
//					+ ".new.bam", outputPrefix + ".unaligned.fastq", 10);
//			renameFile(outputPrefix + ".new.bam", outputPrefix + ".bam");

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
					"IO exception thrown in PARA-suite with "
							+ "following error message:"
							+ System.getProperty("line.separator")
							+ e.getMessage());
			e.printStackTrace();
		} catch (InterruptedException e) {
			MappingLogger.getLogger().error(
					"Interrupted exception thrown in PARA-suite "
							+ "with following error message:"
							+ System.getProperty("line.separator")
							+ e.getMessage());
			e.printStackTrace();
		}
	}
}
