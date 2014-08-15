Introduction
------------
ScanIndel finds indels (insertions and deletions) smaller than the length of a short read by re-align soft clipped reads. The workflow of ScanIndel is BWA-mem -> BLAT remap softclipped read -> FreeBayes for indel calling. In addition, ScanIndel can also predict SNPs (single-nucleotide polymorphisms) and MNPs (multi-nucleotide polymorphisms) in a seperate VCF output.

Pre-installation
----------------
Softwares:
* SAMtools (http://samtools.sourceforge.net/)
* BEDtools (https://code.google.com/p/bedtools/)
* BWA (http://bio-bwa.sourceforge.net/) 
* BLAT (http://genome.ucsc.edu/FAQ/FAQblat.html)
* freebayes (https://github.com/ekg/freebayes)
* vcflib (https://github.com/ekg/vcflib) 
* vt (https://github.com/atks/vt)
* pysam python package (https://code.google.com/p/pysam/)

All softwares above are assumed to be installed in your searching path. Ask your admistrator for assistance if necessary.

__Note__: if using BAM file as input, BWA is not required for running ScanIndel

Getting Soure Code
------------------
	git clone git://github.com/cauyrd/ScanIndel.git
	cd ScanIndel
Running ScanIndel
-----------------
### command-line usage
	python ScanIndel.py -i sample.txt -c config.txt [options]
#### Options:
	 -f				:min-alternate-fraction for FreeBayes (default 0.2)
	 -s  			:softclipping percentage triggering BLAT re-alignment (default 0.2)
	 -l  			:minimal length of indels to be included in final output (default 4)
	 -d  			:minimal sequencing depth to indetify variants (default 20)
	 -t  			:limit analysis to targets listed in a provided BED-format file
	 --bam 			:the input file is BAM format
	 -h --help 		:produce this menu
#### Input:
	sample.txt    	:this file contains the listed samples to be analyzed (one per line), the input can be raw read FastQ file or aligned BAM file and use --bam when running (default name is sample.txt)
	config.txt    	:this file contains the path of reference file for each BWA, BLAT and Freebayes (default name is config.txt)
#### Output:
The output files include the VCF file for detected variant and BAM files for BWA-MEM and BLAT mapping.

	*.indel.vcf	:VCF file includes putative INDELs or complex event (composite indel and substitution events).
	*.snp.vcf	:VCF file includes putative SNPs
	*.sorted.bam	:BAM file from BWA-MEM if using FastQ as input
	*.softclip_realigned.bam	:BAM file after BLAT re-aligning soft-clipped reads from BWA-MEM or input BAM file
Example Data
------------
The folder example contains the examples of sample.txt and config.txt and VCF output by ScanIndel.

