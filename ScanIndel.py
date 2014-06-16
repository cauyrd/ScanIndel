#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: ScanIndel.py
#
#        USAGE: ./ScanIndel.py sample.txt(optional) config.txt(optional)
#
#  DESCRIPTION: variant detection for oncogenes 
#
#      OPTIONS: ---
# REQUIREMENTS: bwa, samtools, blat, freebayes, bamtools
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: Fri May  2 14:37:20 CDT 2014
#     REVISION: ---
#===============================================================================
import sys
import os
from time import time, strftime

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

def main():
	# timing progrom start
	start = time()
	print 'BBF starts running: '+strftime("%Y-%m-%d %H:%M:%S")
	print "Data input"
	if len(sys.argv) > 1:
		sample_file = sys.argv[1]
		config_file = sys.argv[2]
	else:
		sample_file = "sample.txt"
		config_file = "config.txt"
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
		print "Mapping with BWA-mem" 
		os.system("bwa mem -t8 "+reference['bwa']+" "+sample[each]+" >"+each+".sam")
		os.system("samtools view -bS "+each+".sam >"+each+".bam")
		os.system("samtools sort "+each+".bam "+each+".sorted")
		os.system("samtools index "+each+".sorted.bam")
		os.remove(each+".sam")
		os.remove(each+".bam")
		
		print "Blat remapping softclipped reads"
		os.system('python '+path+'/script/bwa2blat.py '+reference['blat']+' '+each+'.sorted.bam')

		print "Freebayes variant calling"
		os.system('freebayes -F 0.005 -f '+reference['freebayes']+' '+each+'.sorted.bam.blat.bam | vcffilter -f "( LEN > 4 & TYPE = ins ) | ( LEN > 4 & TYPE = del )" >'+each+'.scanindel.raw.vcf')
		os.system('intersectBed -a '+each+'.scanindel.raw.vcf -b '+path+'/bedfile/hg19_coding_exon.bed -wa -header >'+each+'.scanindel.exon.vcf')
		os.remove(each+'.scanindel.raw.vcf')
	print "Cleanup sever"
	os.system('gfServer stop localhost 50000')
	print "BBF finish running: "+strftime("%Y-%m-%d %H:%M:%S")
	end = time()
	print 'Programing runing: '+str(end-start)+' seconds.'

if __name__ == '__main__':
	main()
