Introduction
------------
ScanIndel pipeline find indels by re-align soft clipped reads. ScanIndel workflow: BWA-mem -> BLAT remap softclipped read -> FreeBayes for indel calling.

Pre-installation
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

Getting Soure Code
-------------------
	git clone git://github.com/cauyrd/ScanIndel.git
	cd ScanIndel
Running SVfinder
----------------
### command-line usage
	python ScanIndel.py -i sample.txt -c config.txt [options]
#### Options:
	 -f				:min-alternate-fraction for FreeBayes (default 0.2)
	 -s  			:softclipping percentage triggering BLAT re-alignment (default 0.2)
	 -l  			:minmal length of indels to be included in final output (default 4)
	 -v  			:minmal alternate fraction for other type of variants [snp, mnp, complex] except indels (default 0.1)
	 -h  			:produce this menu
#### Input:
	sample.txt    	:this file contains the name of sample and the name of input raw read files (default name is sample.txt)
	config.txt    	:this file contains the path of reference file for each BWA, BLAT and Freebayes (default name is config.txt)
#### Output:
The output files include the VCF file for detected variant and BAM files for BWA-MEM and BLAT mapping.

	*.indels.exon.vcf	:VCF file includes the indels in targeted exonic regions
	*.others.exon.vcf	:VCF file includes snp, mnp and complex events in targeted exonic regions
	*.sorted.bam		:BAM file from BWA-MEM
	*.sorted.bam.blat.bam : BAM file after BLAT re-aligning soft-clipped reads in BWA-MEM bam file
Example
-------------
The folder example contains the examples of sample.txt and config.txt and VCF output.

