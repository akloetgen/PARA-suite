package utils.postprocessing;

import main.MappingLogger;
import mapping.Bowtie2Mapping;

public class MainPostprocessingTest {

	public static void main(String args[]) {

		MappingLogger.getLogger().debug("bla bla test start");
		CombineGenomeTranscript test = new CombineGenomeTranscript();
		
		// genomeMappingFile, transcriptMappingFile, combinedMappingFile
		test.combine(args[0], args[1], args[2]);
		Bowtie2Mapping bowtieMapping = new Bowtie2Mapping();
		bowtieMapping.sortByCoordinateAndIndex(args[2]);
	}

}
