# PARA-suite toolkit
PAR-CLIP Analyzing suite. Useful tools for short and error prone sequence read handling. Note, that the PARA-suite addon of the Burrows-Wheeler Aligner (BWA) is necessary for the mapping tool of the PARA-suite.

For more information on the usage please see "Manual.pdf".

## Installation
	git clone https://github.com/akloetgen/PARA-suite.git
	cd PARA-suite/bin/
	java -jar parasuite.jar setup --parasuite $myPATH_TO_PARASUITE

The 3rd installation step sets the path to the binary-folder of the PARA-suite alignment algorithm to make it accessible for the PARA-suite toolkit. Therefore, make sure the PARA-suite alignment algorithm (https://github.com/akloetgen/PARA-suite_aligner) is installed correctly.

### Requirements
- Java (7 or 8 should work)
- Perl 5
- samtools (Version 1.0 or newer: https://github.com/samtools/samtools) 
- Perl CPAN Math::Random package (http://search.cpan.org/~grommel/Math-Random-0.71/)

### Getting help
	java -jar parasuite.jar --help

## Testdata
	chmod u+x examples.sh examples_remove_temp.sh
	./examples.sh
	./examples_remove_temp.sh

If you have installed the Perl CPAN Math::Random package locally, you have to specify the path to the installed package in the examples.sh script at the top like in the following example:
	myPathToMathRandom=~/.cpan/installed/lib/perl
	
