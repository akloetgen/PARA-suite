package utils.postprocessing;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;

/**
 * Creates a new unaligned fastq file. uses a BAM file and a fastq file as input
 * and extracts all reads from the fastq file that are not aligned within the
 * BAM file.
 * 
 * @author akloetgen
 * 
 */
public class ExtractNotMappedReads {

	/**
	 * 
	 * @param mappingFileName
	 *            mapping file containing no unaligned sequencing reads
	 * @param alignedReadFileName
	 *            file name where all aligned reads are saved to
	 */
	public void extractReads(String mappingFileName, String alignedReadFileName) {
		final SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
				.validationStringency(ValidationStringency.LENIENT);
		SamReader mappingFileReader = factory.open(new File(mappingFileName));
		// File readFile = new File(readFileName);

		File alignedReadFile = new File(alignedReadFileName);
		// List<String> alignedReadNames = new LinkedList<String>();

		// SAMFileWriterFactory fac = new SAMFileWriterFactory();
		// SAMFileWriter mappingFileNew = fac.makeBAMWriter(mappingFileReader
		// .getFileHeader(), false, new File(mappingFileNameNew));

		try {
			Writer unalignedReadOutput = new BufferedWriter(
					new OutputStreamWriter(
							new FileOutputStream(alignedReadFile)));

			for (SAMRecord readHit : mappingFileReader) {
				// if (readHit.getMappingQuality() < mapqThreshold) {
				// print to FASTQ-File

				// alignedReadNames.add("@" + readHit.getReadName());

				String readHeader = "@" + readHit.getReadName();
				// String readSequence = readHit.getReadString();
				// String plus = "+";
				// String readQuality = readHit.getBaseQualityString();

				// if (readHit.getReadNegativeStrandFlag()) {
				// String tempReadSequence = htsjdk.samtools.util.SequenceUtil
				// .reverseComplement(readSequence);
				// readSequence = tempReadSequence;
				// String tempQuality = new StringBuffer(readQuality)
				// .reverse().toString();
				// readQuality = tempQuality;
				// }
				unalignedReadOutput.write(readHeader
						+ System.getProperty("line.separator"));

				// } else {
				// // good mapping, add to new mapping file with high MAPQ
				// mappingFileNew.addAlignment(readHit);
				// }
			}

			// mappingFileNew.close();
			unalignedReadOutput.close();
			mappingFileReader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
