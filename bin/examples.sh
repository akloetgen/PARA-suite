#!/bin/bash

# if you have installed the CPAN Math::Random library somewhere else, please specify the location here
#myPathToMathRandom=/etc/perl/
myPathToMathRandom=/home/akloetgen/.cpan/installed/lib/perl/5.14.2/

# executing the PAR-CLIP read simulator script using the example files in /../examples/simulation/
java -jar parma.jar simulate ../examples/simulation/reference_chr1_transcripts.fa ../examples/simulation/testout_simulation ../examples/simulation/example.errorprofile ../examples/simulation/example.sitefrequency ../examples/simulation/example.sitepositions ../examples/simulation/example.qualities ../examples/simulation/example.indels 0.6 -I $myPathToMathRandom

# run PARMA pipeline on simulated PAR-CLIP reads
bwa index ../examples/references/reference_chr1.fa
samtools faidx ../examples/references/reference_chr1.fa
bwa index ../examples/references/reference_chr1_transcripts.fa
java -jar parma.jar map -q ../examples/simulation/testout_simulation.fastq -r ../examples/references/reference_chr1.fa -t ../examples/references/reference_chr1_transcripts.fa -o ../examples/mapping/testout_simulation_mapped --refine

# benchmark PARMA alignment on simulated PAR-CLIP reads
java -jar parma.jar benchmark ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/mapping/testout_simulation_mapped.combined.stats ../examples/simulation/testout_simulation.fastq

# calculate error profile
java -jar parma.jar error ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/references/reference_chr1.fa 51 -q

# combine genomic and transcriptomic mapping
java -jar parma.jar comb ../examples/mapping/testout_simulation_mapped.PARMA-genomic.bam ../examples/mapping/testout_simulation_mapped.PARMA-transcript.bam ../examples/mapping/testout_simulation_mapped.PARMA.combined.bam

# clustering aligned PAR-CLIP reads to obtain RBP-bound regions
java -jar parma.jar clust ../examples/mapping/testout_simulation_mapped.combined.bam ../examples/references/reference_chr1.fa ../examples/mapping/testout_simulation_mapped.combined.clusters ../examples/references/snp_db.vcf.gz 1
