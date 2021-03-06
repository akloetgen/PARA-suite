package mapping;

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
			String additionalOptions) {

		MappingLogger.getLogger().info(
				"Starting user mapping to investigate error-profile");
		userAlignerEntireCommand = userAlignerEntireCommand.replace(
				"REFERENCE", reference);
		userAlignerEntireCommand = userAlignerEntireCommand.replace("INPUT",
				input);
		userAlignerEntireCommand = userAlignerEntireCommand.replace("OUTPUT",
				outputPrefix + ".sam");
		userAlignerEntireCommand = userAlignerEntireCommand.replace("THREADS",
				"" + threads);

		MappingLogger.getLogger().debug(
				"USER COMMAND:" + userAlignerEntireCommand);

		String[] userAlignerCommands = userAlignerEntireCommand.split("\\s");

		List<String> userAlignerCommandsList = new LinkedList<String>();
		for (int i = 0; i < userAlignerCommands.length; i++) {
			userAlignerCommandsList.add(userAlignerCommands[i]);
		}

		executeCommand(userAlignerCommandsList, StreamRedirect.ALL);

		MappingLogger
				.getLogger()
				.info("Convert SAM-file of user mapped reads to BAM-file, filter, sort, index, remove temp files");
		MappingLogger.getLogger().debug(
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
