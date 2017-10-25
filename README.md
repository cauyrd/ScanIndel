Introduction
------------
ScanIndel is a python program to detect indels (insertions and deletions) from NGS data by re-align and de novo assemble soft clipped reads. 

Prerequisites
----------------
Softwares and Python packages:
* BedTools/2.17.0 (https://github.com/arq5x/bedtools2/releases)
* SAMtools/1.0 or less than 1.0 (http://samtools.sourceforge.net/)
* BWA/0.7.10 (http://bio-bwa.sourceforge.net/) 
* BLAT (gfServer and gfClient)/34+ (http://genome.ucsc.edu/FAQ/FAQblat.html)
* freebayes/0.9.18 (https://github.com/ekg/freebayes)
* Inchworm assembler (http://inchworm.sourceforge.net)
* Python/2.7 (https://www.python.org/)
* Pysam/0.7.7 (https://code.google.com/p/pysam/)
* PyVCF/0.6.7 (https://github.com/jamescasbon/PyVCF)
* Biopython/1.64 (http://biopython.org/wiki/Main_Page)
* SciPy/0.14.0 and NumPy/1.8.1 (http://www.scipy.org/)

All softwares above are assumed to be installed in your searching path. Ask your admistrator for assistance if necessary. 

Getting Soure Code
------------------
	git clone git://github.com/cauyrd/ScanIndel.git
	cd ScanIndel
Running ScanIndel
-----------------
### command-line usage
	python ScanIndel.py -i sample.txt -p config.txt [options]
#### Options:
	 -F							:setting min-alternate-fraction for FreeBayes (default 0.2)
	 -o							:setting output directory (default ./)
	 -C							:setting min-alternate-count for FreeBayes (default 2)
	 -s							:softclipping percentage triggering BLAT re-alignment (default 0.2)
	 -t							:setting -t for FreeBayes to provide a BED-format file limiting the analysis to these regions
	 --min_percent_hq			:min percentage of high quality base in soft clipping reads (default 0.8)
	 --lowqual_cutoff			:low quality cutoff value (default 20)
	 --mapq_cutoff				:low mapping quality cutoff (default 1)
	 --blat_ident_pct_cutoff	:Blat sequence identity cutoff (default 0.8)
	 --gfServer_port			:gfServer service port number, changing this value to allow multiple ScanIndel running at a single machine (default 50000)
	 --hetero_factor			:The factor about the indel heterogenirity and heterozygosity (default 0.1)
	 --bam 						:the input file is BAM format
	 --rmdup					:exccute duplicate removal step before realignment
	 -h --help					:produce this menu
	 -v --version				:show version of this tool
#### Input:
	sample.txt    				:this file contains the listed samples to be analyzed (one per line), the input can be raw read FastQ file or aligned BAM file and use --bam when running (default name is sample.txt)
	config.txt    				:this file contains the path of reference file for each BWA, BLAT and Freebayes (default name is config.txt)
#### Output:
The output files include the VCF file for detected variant and BAM files for BWA-MEM and BLAT mapping.

	*.reads.bam					:BAM file for read after blat alignment.
	*.contigs.bam 				:BAM file for de novo assembled contigs after BWA and BLAT mapping.
	*.mapping.indel.vcf			:VCF file includes putative INDELs from softclipping read re-alignment.
	*.assembly.indel.vcf		:VCF file includes putative INDELs from de novo assembly.
	*.merged.indel.vcf			:VCF file that include all putative INDELs by merging the results from *mapping.indel.vcf and *.assembly.indel.vcf
#### Example:
To run a test example, please go to folder 'example' and follow the README file to run test data.
#### Citing ScanIndel:
Yang, Rendong, et al. "ScanIndel: a hybrid framework for indel detection via gapped alignment, split reads and de novo assembly." *Genome medicine* 7.1 (2015): 1-12.
