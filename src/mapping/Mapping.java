package mapping;

import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.LinkedList;
import java.util.List;

import main.MappingLogger;
import enums.StreamRedirect;

/**
 * Some general commands useful for read alignment post-procesing, such as
 * sorting and indexing.
 * 
 * @author akloetgen
 * 
 */
public abstract class Mapping {

	protected Process mappingProcess;
	protected ProcessBuilder mappingProcessBuilder;
	protected long timeStart;

	/**
	 * executes mapping program and all following commands to sort, filter and
	 * index the resulting BAM-file so that it can be used for further analysis
	 * such as error-profiling.
	 * 
	 * @param threads
	 *            number of threads that can be used for bowtie2
	 * @param reference
	 *            filename of reference sequence, must be indexed
	 * @param input
	 *            filename of read file for input
	 * @param outputPrefix
	 *            prefix for filename of the mapping result
	 * @throws MappingErrorException
	 *             if mapper stops with an error, an exception is thrown.
	 */
	public abstract void executeMapping(int threads, String reference,
			String input, String outputPrefix, int mappingQualityFilter,
			String additionalOptions) throws MappingErrorException;

	/**
	 * Renames a file in the file system.
	 * 
	 * @param oldFileName
	 *            old file name
	 * @param newFileName
	 *            new file name
	 */
	public void renameFile(String oldFileName, String newFileName) {
		try {
			List<String> renameCommandList = new LinkedList<String>();
			renameCommandList.add("mv");
			renameCommandList.add(oldFileName);
			renameCommandList.add(newFileName);
			executeCommand(renameCommandList, StreamRedirect.STDOUT);
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Removes a file from the file system.
	 * 
	 * @param removeFileName
	 *            filename of the file which should be removed
	 */
	public void removeFile(String removeFileName) {
		try {
			List<String> removeCommandList = new LinkedList<String>();
			removeCommandList.add("rm");
			removeCommandList.add(removeFileName);
			executeCommand(removeCommandList, StreamRedirect.STDOUT);
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Sorts a BAM file by alignment coordinate and creates an index. Requires
	 * samtools to be installed.
	 * 
	 * @param bamFileName
	 *            BAM file which should be sorted and indexed.
	 * 
	 */
	public void sortByCoordinateAndIndex(String bamFileName) {
		try {
			MappingLogger.getLogger().debug(
					"Sorting BAM-file for faster access");
			List<String> sortCommandsList = new LinkedList<String>();
			sortCommandsList.add("samtools");
			sortCommandsList.add("sort");
			sortCommandsList.add(bamFileName);
			sortCommandsList.add(bamFileName + "sort");
			executeCommand(sortCommandsList, StreamRedirect.STDOUT);
			sortCommandsList.clear();
			sortCommandsList.add("mv");
			sortCommandsList.add(bamFileName + "sort.bam");
			sortCommandsList.add(bamFileName);
			executeCommand(sortCommandsList, StreamRedirect.STDOUT);
			MappingLogger.getLogger().debug(
					"Indexing BAM-file for faster access");
			sortCommandsList.clear();
			sortCommandsList.add("samtools");
			sortCommandsList.add("index");
			sortCommandsList.add(bamFileName);
			executeCommand(sortCommandsList, StreamRedirect.STDOUT);
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Sorts a BAM file by read name and creates an index. Requires samtools to
	 * be installed.
	 * 
	 * @param bamFileName
	 *            BAM file which should be sorted and indexed.
	 * 
	 */
	public void sortByNameAndIndex(String bamFileName) {
		try {
			MappingLogger.getLogger().debug(
					"Sorting BAM-file for faster access");
			List<String> sortCommandsList = new LinkedList<String>();
			sortCommandsList.add("samtools");
			sortCommandsList.add("sort");
			sortCommandsList.add("-n");
			sortCommandsList.add(bamFileName);
			sortCommandsList.add(bamFileName + "sort");
			executeCommand(sortCommandsList, StreamRedirect.STDOUT);
			sortCommandsList.clear();
			sortCommandsList.add("mv");
			sortCommandsList.add(bamFileName + "sort.bam");
			sortCommandsList.add(bamFileName);
			executeCommand(sortCommandsList, StreamRedirect.STDOUT);
			MappingLogger.getLogger().debug(
					"Indexing BAM-file for faster access");
			sortCommandsList.clear();
			sortCommandsList.add("samtools");
			sortCommandsList.add("index");
			sortCommandsList.add(bamFileName);
			executeCommand(sortCommandsList, StreamRedirect.STDOUT);
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Executes a command on the shell.
	 * 
	 * @param commands
	 *            List of commands and arguments
	 * @param streamRedirect
	 *            if a stream should be redirected from the program output to
	 *            the console
	 * @return returns the return-value of the respective program
	 * @throws InterruptedException
	 * @throws IOException
	 */
	protected int executeCommand(List<String> commands,
			StreamRedirect streamRedirect) throws InterruptedException,
			IOException {
		String commandString = "";
		for (int i = 0; i < commands.size(); i++) {
			commandString += commands.get(i) + " ";
		}
		MappingLogger.getLogger().debug(commandString);
		mappingProcessBuilder = new ProcessBuilder(commands);
		if (streamRedirect.equals(StreamRedirect.STDOUT)
				|| streamRedirect.equals(StreamRedirect.ALL)) {
			mappingProcessBuilder.redirectOutput(Redirect.INHERIT);
		}
		if (streamRedirect.equals(StreamRedirect.ERROR)
				|| streamRedirect.equals(StreamRedirect.ALL)) {
			mappingProcessBuilder.redirectError(Redirect.INHERIT);
		}
		mappingProcess = mappingProcessBuilder.start();
		return mappingProcess.waitFor();
	}

	protected void setTimeStart() {
		timeStart = System.currentTimeMillis();
	}

	protected long calculatePassedTime() {
		return (long) (System.currentTimeMillis() - timeStart) / 1000;
	}
}
