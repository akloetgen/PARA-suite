#!/bin/bash

# if you have installed the CPAN Math::Random library somewhere else, please specify the location here
myPathToMathRandom=/etc/perl/

if [ -s ../examples/simulation/testout_simulation.fastq ]
then
	rm ../examples/simulation/testout_simulation.fastq
fi

touch ../examples/simulation/testout_simulation.fastq

until [[ -s ../examples/simulation/testout_simulation.fastq ]]
do
	# executing the PAR-CLIP read simulator script using the example files in ../examples/simulation/
	java -jar parasuite.jar simulate ../examples/references/reference_chr1_transcripts.fa ../examples/simulation/testout_simulation ../examples/simulation/example.errorprofile ../examples/simulation/example.sitefrequency ../examples/simulation/example.sitepositions ../examples/simulation/example.qualities ../examples/simulation/example.indels 0.6 -I $myPathToMathRandom
	
	if [ $? != 0 ]; then
		echo "Simulation failed. Please refer to error message of the simulator and re-run this pipeline when the error is fixed."
		exit 1
	fi
	if [ ! -s ../examples/simulation/testout_simulation.fastq ]; then  
		rm ../examples/simulation/testout_simulation.fastq
		touch ../examples/simulation/testout_simulation.fastq
	fi
done

# run PARAsuite pipeline on simulated PAR-CLIP reads
if [ ! -s ../examples/references/reference_chr1.fa.fai ]; then
	samtools faidx ../examples/references/reference_chr1.fa
	if [ $? != 0 ]; then
		echo "FASTA index creation failed. Please make sure samtools is properly installed and re-run this pipeline."
	exit 1
fi
fi

java -jar parasuite.jar map -q ../examples/simulation/testout_simulation.fastq -r ../examples/references/reference_chr1.fa -t ../examples/references/reference_chr1_transcripts.fa -o ../examples/mapping/testout_simulation_mapped --refine
if [ $? != 0 ]; then
	echo "Mapping failed. Please refer to the mappings log-file and re-run this pipeline when the error is fixed."
	exit 1
fi

# benchmark PARAsuite alignment on simulated PAR-CLIP reads
java -jar parasuite.jar benchmark ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/mapping/testout_simulation_mapped.combined.stats ../examples/simulation/testout_simulation.fastq
if [ $? != 0 ]; then
	echo "Benchmarking failed. Please refer to the last error message and re-run this pipeline when the error is fixed."
	exit 1
fi

# calculate error profile
java -jar parasuite.jar error ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/references/reference_chr1.fa 51
if [ $? != 0 ]; then
	echo "Error profile calculation failed. Please refer to the last error message and re-run this pipeline when the error is fixed."
	exit 1
fi

# combine genomic and transcriptomic mapping
java -jar parasuite.jar comb ../examples/mapping/testout_simulation_mapped.PARAsuite-genomic.bam ../examples/mapping/testout_simulation_mapped.PARAsuite-transcript.bam ../examples/mapping/testout_simulation_mapped.PARAsuite.combined.bam
if [ $? != 0 ]; then
	echo "Combining results failed. Please refer to the last error message and re-run this pipeline when the error is fixed."
	exit 1
fi

# clustering aligned PAR-CLIP reads to obtain RBP-bound regions
java -jar parasuite.jar clust ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/references/reference_chr1.fa ../examples/mapping/testout_simulation_mapped.combined.clusters ../examples/references/snp_db.vcf.gz 1
if [ $? != 0 ]; then
	echo "Clustering failed. Please refer to the last error message and re-run this pipeline when the error is fixed."
	exit 1
fi
