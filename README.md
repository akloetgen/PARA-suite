# PARMA toolkit
PAR-CLIP Read MApper toolkit. Useful tools for short and error prone sequence read handling. Note, that the PARMA addon of the Burrows-Wheeler Aligner (BWA) is necessary for the mapping tool of the PARMA-tk.

For more information on the usage please see "Manual.pdf".

## Installation
	git clone https://github.com/akloetgen/PARMA_tk.git
	cd PARMA_tk/bin/
	java -jar parma.jar setup --parma $myPATH_TO_PARMA

The 3rd installation step sets the path to the binary-folder of the PARMA algorithm to make it accessible for the PARMA toolkit. Therefore, make sure the PARMA algorithm (https://github.com/akloetgen/PARMA) is installed correctly.

### Requirements
- Java (7 or 8 should work)
- Perl 5
- samtools (Version 1.0 or newer: https://github.com/samtools/samtools) 
- Perl CPAN Math::Random package (http://search.cpan.org/~grommel/Math-Random-0.71/)

### Getting help
	java -jar parma.jar --help

## Testdata
	chmod u+x examples.sh examples_remove_temp.sh
	./examples.sh
	./examples_remove_temp.sh

If you have installed the Perl CPAN Math::Random package locally, you have to specify the path to to installed package in the examples.sh script at the top.

