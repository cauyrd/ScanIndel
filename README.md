Introduction
------------
BBF pipeline find small polymorphisms, specially indels, SNPs or MNPs for cancer requencing samples. BBF is named by the first letter of three major software: **B**wa, **B**lat and **F**reebayes used.

Pre-installtalation
-------------------
Python Packages:
* samtools (http://samtools.sourceforge.net/)
* bamtools (https://github.com/pezmaster31/bamtools)
* IGV (https://www.broadinstitute.org/igv/home)
* BWA (http://bio-bwa.sourceforge.net/) 
* BLAT (http://genome.ucsc.edu/FAQ/FAQblat.html)
* freebayes (https://github.com/ekg/freebayes)

Running SVfinder
----------------
### command-line usage
	python bbf.py sample.txt(optional) config.txt(optional)
#### Options:
	sample.txt    :this file contains the name of sample and the name of input raw read files (default name is sample.txt)
	config.txt    :this file contains the path of reference file for each BWA, BLAT and Freebayes (default name is config.txt)

Output
-------------
The output files include the VCF file for detected variant and BAM files for BWA and BLAT mapping (\*.blat.bam).

Example
-------------
The folder example contains the examples of sample.txt and config.txt and VCF output.

