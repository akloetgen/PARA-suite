package mapping;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import main.MappingLogger;
import enums.StreamRedirect;

/**
 * 
 * @author akloetgen
 * 
 */
public class UserMapping extends Mapping {

	private String userAlignerEntireCommand;

	public UserMapping(String userAlignerCommand) {
		this.userAlignerEntireCommand = userAlignerCommand;
	}

	public void executeMapping(int threads, String reference, String input,
			String outputPrefix, int mappingQualityFilter,
			String additionalOptions) throws MappingErrorException {
		try {
			MappingLogger.getLogger().info(
					"Starting user mapping to investigate error-profile");

			String[] userAlignerCommands = userAlignerEntireCommand
					.split("\\s");
			List<String> userAlignerCommandsList = new LinkedList<String>();
			for (int i = 0; i < userAlignerCommands.length; i++) {
				userAlignerCommandsList.add(userAlignerCommands[i]);
			}

			// MappingLogger.getLogger().debug(bowtie2Command);
			if (executeCommand(userAlignerCommandsList, StreamRedirect.NOTHING) != 0) {
				MappingErrorException e = new MappingErrorException();
				e.setMappingCommand(userAlignerCommandsList);
				throw e;
			}

			MappingLogger.getLogger().info(
					"Convert SAM-file of usermapping mapped reads to BAM-file");
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
