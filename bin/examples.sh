#!/bin/bash

# if you have installed the CPAN Math::Random library somewhere else, please specify the location here
myPathToMathRandom=/etc/perl/

# executing the PAR-CLIP read simulator script using the example files in /../examples/simulation/
java -jar parasuite.jar simulate ../examples/references/reference_chr1_transcripts.fa ../examples/simulation/testout_simulation ../examples/simulation/example.errorprofile ../examples/simulation/example.sitefrequency ../examples/simulation/example.sitepositions ../examples/simulation/example.qualities ../examples/simulation/example.indels 0.6 -I $myPathToMathRandom

if ![[ -s ../examples/simulation/testout_simulation.fastq ]] ; then
	echo "../examples/simulation/testout_simulation.fastq is empty."
	echo "please run examples_remove_temp.sh and re-run examples.sh"
	echo "this error can be caused when no read was created due to the random nature of the simulation"
	exit
fi ;

# run PARAsuite pipeline on simulated PAR-CLIP reads
samtools faidx ../examples/references/reference_chr1.fa
#ln -s reference_chr1.fa.fai ../examples/references/reference_chr1_bwa.fa.fai
#ln -s reference_chr1.fa.fai ../examples/references/reference_chr1_PARAsuite.fa.fai

#ln -s reference_chr1.fa ../examples/references/reference_chr1_bwa.fa
#ln -s reference_chr1.fa ../examples/references/reference_chr1_PARAsuite.fa

java -jar parasuite.jar map -q ../examples/simulation/testout_simulation.fastq -r ../examples/references/reference_chr1.fa -t ../examples/references/reference_chr1_transcripts.fa -o ../examples/mapping/testout_simulation_mapped --refine

# benchmark PARAsuite alignment on simulated PAR-CLIP reads
java -jar parasuite.jar benchmark ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/mapping/testout_simulation_mapped.combined.stats ../examples/simulation/testout_simulation.fastq

# calculate error profile
java -jar parasuite.jar error ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/references/reference_chr1.fa 51 -q

# combine genomic and transcriptomic mapping
java -jar parasuite.jar comb ../examples/mapping/testout_simulation_mapped.PARAsuite-genomic.bam ../examples/mapping/testout_simulation_mapped.PARAsuite-transcript.bam ../examples/mapping/testout_simulation_mapped.PARAsuite.combined.bam

# clustering aligned PAR-CLIP reads to obtain RBP-bound regions
java -jar parasuite.jar clust ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/references/reference_chr1.fa ../examples/mapping/testout_simulation_mapped.combined.clusters ../examples/references/snp_db.vcf.gz 1
