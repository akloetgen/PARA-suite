package utils.pileupclusters;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FetchSequencesForBEDFile {

	private IndexedFastaSequenceFile reference;

	public void fetchSequences(String referenceFile, String bindingSitesFile,
			String outputFile) {

		try {
			reference = new IndexedFastaSequenceFile(new File(referenceFile));
			FileWriter fileWriter = new FileWriter(outputFile);
			// File bindingSites = new File(bindingSitesFile);
			BufferedReader br = new BufferedReader(new FileReader(
					bindingSitesFile));
			String header = br.readLine();
			fileWriter.write(header + System.getProperty("line.separator"));

			String bindingSite;
			int passed = 0;

			while ((bindingSite = br.readLine()) != null) {

				String[] splittedLine = bindingSite.split("\t");

				// System.out.println("1=" + splittedLine[1] + "; 2="
				// + splittedLine[2] + "; 3=" + splittedLine[3] + "; 4="
				// + splittedLine[4]);

				String chr = splittedLine[0];
				if (!chr.startsWith("chr")) {
					chr = "chr" + chr;
				}
				int start = Integer.parseInt(splittedLine[1]);
				int end = Integer.parseInt(splittedLine[2]);
				// boolean isReverse = false;
				StrandOrientation isReverse = new StrandOrientation();
				isReverse.setReverse(false);
				if (splittedLine[4].equals("-")) {
					isReverse.setReverse(true);
				}

				byte[] bindingSiteSequenceBytes;
				try {
					bindingSiteSequenceBytes = reference.getSubsequenceAt(chr,
							start, end).getBases();
					if (isReverse.getStrandOrientation().equals("-")) {
						htsjdk.samtools.util.SequenceUtil
								.reverseComplement(bindingSiteSequenceBytes);
					}
				} catch (SAMException e) {
					bindingSiteSequenceBytes = new byte[0];
				}

				StringBuffer bindingSiteSequence = new StringBuffer();
				for (int basePosition = 0; basePosition < bindingSiteSequenceBytes.length; basePosition++) {
					bindingSiteSequence
							.append(Character
									.toUpperCase((char) bindingSiteSequenceBytes[basePosition]));
				}

				// splittedLine[11] = bindingSiteSequence.toString();

				fileWriter.write(">" + splittedLine[3]
						+ System.getProperty("line.separator")
						+ bindingSiteSequence.toString()
						+ System.getProperty("line.separator"));
				passed++;
				// if (passed % 1000 == 0) {
				// System.out.println("passed lines = " + passed);
				//
				// }

			}
			br.close();
			fileWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
