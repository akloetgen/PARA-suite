# PARMA-tk
PAR-CLIP Read MApper toolkit. Useful tools for short and error prone sequence read handling. Note, that the PARMA addon of the Burrows-Wheeler Aligner (BWA) is necessary for the mapping tool of the PARMA-tk.

For more information on the usage please see "Manual.pdf".

###Installation
	git clone https://github.com/akloetgen/PARMA_tk.git
	cd PARMA_tk/bin/
	java -jar parma.jar setup --parma $myPATH_TO_PARMA
	
Make sure the PARMA algorithm (https://github.com/akloetgen/PARMA) is installed correctly
	