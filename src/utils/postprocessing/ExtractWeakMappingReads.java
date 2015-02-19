package utils.postprocessing;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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
 * 
 * 
 * 
 * @author akloetgen
 * 
 */
public class ExtractWeakMappingReads {

	public ExtractWeakMappingReads() {

	}

	public void extractReads(String mappingFileName, String mappingFileNameNew,
			String unalignedReadFileName, int mapqThreshold) {
		// File mappingFile = new File(mappingFileName);
		final SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
				.validationStringency(ValidationStringency.LENIENT);
		SamReader mappingFileReader = factory.open(new File(mappingFileName));

		// SAMFileReader mappingFileReader = new SAMFileReader(mappingFile);
		File unalignedReadFile = new File(unalignedReadFileName);
		SAMFileWriterFactory fac = new SAMFileWriterFactory();
		SAMFileWriter mappingFileNew = fac.makeBAMWriter(mappingFileReader
				.getFileHeader(), false, new File(mappingFileNameNew));

		try {
			Writer unalignedReadOutput = new BufferedWriter(
					new OutputStreamWriter(new FileOutputStream(
							unalignedReadFile)));

			for (SAMRecord readHit : mappingFileReader) {
				if (readHit.getMappingQuality() < mapqThreshold) {
					// print to FASTQ-File
					String readHeader = "@" + readHit.getReadName();
					String readSequence = readHit.getReadString();
					String plus = "+";
					String readQuality = readHit.getBaseQualityString();

					if (readHit.getReadNegativeStrandFlag()) {
						String tempReadSequence = htsjdk.samtools.util.SequenceUtil
								.reverseComplement(readSequence);
						readSequence = tempReadSequence;
						String tempQuality = new StringBuffer(readQuality)
								.reverse().toString();
						readQuality = tempQuality;
					}

					unalignedReadOutput.write(readHeader
							+ System.getProperty("line.separator")
							+ readSequence
							+ System.getProperty("line.separator") + plus
							+ System.getProperty("line.separator")
							+ readQuality
							+ System.getProperty("line.separator"));
				} else {
					// good mapping, add to new mapping file with high MAPQ
					mappingFileNew.addAlignment(readHit);
				}
			}
			mappingFileNew.close();
			unalignedReadOutput.close();
			mappingFileReader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
