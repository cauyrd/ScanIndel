#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: ScanIndel.py
#
#        USAGE: ./ScanIndel.py -i sample.txt -c config.txt [opts]
#
#  DESCRIPTION: Indel detection for targeted NGS data
#
#      OPTIONS: ---
# REQUIREMENTS: bwa, samtools, blat, freebayes, bedtools
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.1
#      CREATED: Fri May  2 14:37:20 CDT 2014
#     REVISION: Wed Jul 23 14:11:15 CDT 2014
#===============================================================================
import sys
import os
import re
from time import time, strftime
import getopt

def read_config_file(filename):
	path = {}
	ifp = open(filename)
	for line in ifp:
		if line[0] == '#' or line == '\n':
			continue
		item = line.rstrip().split('=')
		path[item[0]]=item[1]
	ifp.close()
	return path

def read_sample_file(filename):
	input = {}
	ifp = open(filename)
	for line in ifp:
		if line[0] == '#' or line == '\n':
			continue
		item = line.rstrip().split()
		input[item[0]] = ' '.join(item[1:])
	ifp.close()
	return input

def parse_vaf(term):
	pattern = re.compile('.+AO=(\d+);.+DP=(\d+);.+')
	item = re.match(pattern,term)
	try:
 		ratio = float(item.group(1))/float(item.group(2))	
	except:
		print 'multiple event detected!'
		return 0
	return ratio

def usage():
	"""helping information"""
	print 'Usage:'
	print ' python ScanIndel.py -c config.txt -i sample.txt [opts]'
	print 'Opts:'
	print ' -f  :min-alternate-fraction for FreeBayes (default 0.2)'
	print ' -s  :softclipping percentage triggering BLAT re-alignment (default 0.2)'
	print ' -l  :minmal length of indels to be included in final output (default 4)'
	print ' -v  :minmal alternate fraction for other type of variants [snp, mnp, complex] except indels (default 0.1)'
	print ' -h  :produce this menu'

def main():

	# parameters parsing
	sample_file = 'sample.txt'
	config_file = 'config.txt'
	freebayes_f = 0.2
	softclip_ratio = 0.2
	indel_len = 4 
	vaf_cutoff = 0.1
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:c:f:s:l:v:h')
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	for o,a in opts:
		if o == '-i':
			sample_file = a
		elif o == '-c':
			config_file = a
		elif o == '-f':
			freebayes_f = float(a)
		elif o == '-s':
			softclip_ratio = float(a)
		elif o == '-l':
			indel_len = int(a)
		elif o == '-v':
			vaf_cutoff = float(a)
		elif o == '-h':
			usage()
			sys.exit(0)
		else:
			assert False, "unhandled option"

	# timing progrom start
	start = time()
	print 'ScanIndel starts running: '+strftime("%Y-%m-%d %H:%M:%S")
	sample = read_sample_file(sample_file)
	reference = read_config_file(config_file)	

	path = os.path.dirname(os.path.realpath(__file__))
	
	print "Start up BLAT server"
	cwd = os.getcwd()
	os.chdir(reference['blat'])
	os.system('gfServer start localhost 50000 *.2bit &')
	os.chdir(cwd)

	for each in sample:
		print "Analyzing sample:", each
		print "STEP1: BWA-MEM start aligning raw reads..." 
		os.system("bwa mem -t8 "+reference['bwa']+" "+sample[each]+" >"+each+".sam")
		os.system("samtools view -bS "+each+".sam >"+each+".bam")
		os.system("samtools sort "+each+".bam "+each+".sorted")
		os.system("samtools index "+each+".sorted.bam")
		os.remove(each+".sam")
		os.remove(each+".bam")
		
		print "STEP2: BLAT start re-aligning soft-clipped reads..."
		os.system('python '+path+'/script/bwa2blat.py '+reference['blat']+' '+each+'.sorted.bam '+str(softclip_ratio))

		print "STEP3: Freebayes start variant calling..."
		os.system('freebayes --pooled-discrete -F '+str(freebayes_f)+' -f '+reference['freebayes']+' '+each+'.sorted.bam.blat.bam > '+each+'.raw.vcf')
		os.system('vcffilter -f "( LEN > '+str(indel_len-1)+' & TYPE = ins ) | ( LEN > '+str(indel_len-1)+' & TYPE = del )" '+each+'.raw.vcf >'+each+'.indels.raw.vcf')
		os.system('vcffilter -f "! ( TYPE = ins ) & ! ( TYPE = del )" '+each+'.raw.vcf >'+each+'.others.raw.vcf')
		os.system('intersectBed -a '+each+'.indels.raw.vcf -b '+path+'/bedfile/hg19_coding_exon.bed -wa -header -u >'+each+'.indels.exon.vcf')
		os.system('intersectBed -a '+each+'.others.raw.vcf -b '+path+'/bedfile/hg19_coding_exon.bed -wa -header -u >'+each+'.others.exon.raw.vcf')
		os.remove(each+'.raw.vcf')
		os.remove(each+'.indels.raw.vcf')
		os.remove(each+'.others.raw.vcf')
		
		# post filtering SNP,MNP and complex variants based on calculated VAF.
		ifp = open(each+'.others.exon.raw.vcf')
		ofp = open(each+'.others.exon.vcf','w')
		for line in ifp:
			if line[0] == '#':
				print >> ofp, line.rstrip()
			else:
				if parse_vaf(line.rstrip()) >= vaf_cutoff:
					print >> ofp, line.rstrip()
		ifp.close()
		ofp.close()
		os.remove(each+'.others.exon.raw.vcf')

	print "Cleanup sever"
	os.system('gfServer stop localhost 50000')
	print "ScanIndel running done: "+strftime("%Y-%m-%d %H:%M:%S")
	end = time()
	print 'Analyzing your data takes: '+str(end-start)+' seconds.'

if __name__ == '__main__':
	main()
