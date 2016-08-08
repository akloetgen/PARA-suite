package utils.pileupclusters;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;

/**
 * 
 * Wrapper for fast lookup of known SNPs in a VCF file.
 * 
 * @author akloetgen
 * 
 */
public class SNPCalling {

	VCFFileReader vcfFileReader;

	/**
	 * Creates object by already loading the VCF file and its index for a fast
	 * lookup of known SNVs.
	 * 
	 * @param vcfFile
	 *            VCF file in VCF format. Tabbed-index must exist with file name
	 *            prefix.tbi.
	 */
	public SNPCalling(String vcfFile) {
		String vcfFileIndex = vcfFile + ".tbi";

		vcfFileReader = new VCFFileReader(new File(vcfFile), new File(
				vcfFileIndex), true);
	}

	/**
	 * Checks whether a given mutation is known as a SNP.
	 * 
	 * @param chr
	 *            Chr of given mutation
	 * @param position
	 *            chromosomal position of given mutation
	 * @param refBase
	 *            base occurring in the reference sequence
	 * @param alternativeBase
	 *            base occurring in the read sequence
	 * @return returns true if the given mutation is also annotated as a known
	 *         SNP
	 */
	public boolean querySNP(String chr, int position, String refBase,
			String alternativeBase) {
		if (chr.startsWith("chr")) {
			String tempChr = chr.substring(3);
			chr = tempChr;
		}

		CloseableIterator<VariantContext> queryResults = vcfFileReader.query(
				chr, position, position + 1);
		while (queryResults.hasNext()) {
			VariantContext variant = queryResults.next();
			if (variant.getChr().equals(chr)
					&& variant.getStart() == position
					&& variant.getReference().getBaseString().contains(refBase)
					&& variant.getAlternateAllele(0).getBaseString()
							.contains(alternativeBase)) {
				return true;
			}
		}
		return false;
	}
}
