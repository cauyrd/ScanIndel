Introduction
------------
ScanIndel pipeline find indels by re-align soft clipped reads. ScanIndel workflow: BWA-mem -> BLAT remap softclipped read -> FreeBayes for indel calling.

Pre-installtalation
-------------------
Softwares:
* samtools (http://samtools.sourceforge.net/)
* bedtools (https://code.google.com/p/bedtools/)
* BWA (http://bio-bwa.sourceforge.net/) 
* BLAT (http://genome.ucsc.edu/FAQ/FAQblat.html)
* freebayes (https://github.com/ekg/freebayes)
* vcffilter from vcflib (https://github.com/ekg/vcflib) 
* pysam python package (https://code.google.com/p/pysam/)

All softwares above are assumed to be installed in your searching path. Ask your admistrator for assistance if necessary.

Running SVfinder
----------------
### command-line usage
	python ScanIndel.py sample.txt(optional) config.txt(optional)
#### Options:
	sample.txt    :this file contains the name of sample and the name of input raw read files (default name is sample.txt)
	config.txt    :this file contains the path of reference file for each BWA, BLAT and Freebayes (default name is config.txt)

Output
-------------
The output files include the VCF file for detected variant and BAM files for BWA and BLAT mapping (\*.blat.bam).

Example
-------------
The folder example contains the examples of sample.txt and config.txt and VCF output.

