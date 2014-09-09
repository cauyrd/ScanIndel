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
# REQUIREMENTS: See README.md file
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.3.1
#      CREATED: Fri May  2 14:37:20 CDT 2014
#     REVISION: Wed Jul 23 14:11:15 CDT 2014
#               Thu Jul 31 16:49:09 CDT 2014
#               Thu Aug  7 13:20:16 CDT 2014
#               Fri Aug  8 15:59:04 CDT 2014
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

def usage():
	"""helping information"""
	print 'Usage:'
	print ' python ScanIndel.py -c config.txt -i sample.txt [opts]'
	print 'Opts:'
	print ' -f  :min-alternate-fraction for FreeBayes (default 0.2)'
	print ' -s  :softclipping percentage triggering BLAT re-alignment (default 0.2)'
	print ' -l  :minimal length of indels to be included in final output (default 4)'
	print ' -d  :minimal sequencing depth to identify variants (default 20)'
	print ' -t  :limit analysis to targets listed in a provided BED-format file'
	print ' --bam  :the input file is BAM format'
	print ' -h --help :produce this menu'
	print 'author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'
	print 'version: v1.3.1'

def main():

	# parameters parsing
	sample_file = 'sample.txt'
	config_file = 'config.txt'
	freebayes_f = 0.2
	softclip_ratio = 0.2
	indel_len = 4 
	depth = 20
	bam = False
	bedfile = ''
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:c:f:s:l:d:t:h', ['bam','help'])
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
		elif o == '-d':
			depth = int(a)
		elif o == '-t':
			bedfile = a
		elif o == '--bam':
			bam = True
		elif o in ('-h', '--help'):
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
		if not bam:
			print "BWA-MEM start aligning raw reads..." 
			os.system("bwa mem -t8 "+reference['bwa']+" "+sample[each]+" >"+each+".sam")
			os.system("samtools view -bS "+each+".sam >"+each+".bam")
			os.system("samtools sort "+each+".bam "+each+".sorted")
			os.system("samtools index "+each+".sorted.bam")
			os.remove(each+".sam")
			os.remove(each+".bam")
			blat_input = each+'.sorted.bam'
		else:
			blat_input = sample[each]
			os.system('samtools index '+blat_input)

		print "BLAT start re-aligning soft-clipped reads..."
		os.system('python '+path+'/script/bwa2blat.py '+reference['blat']+' '+blat_input+' '+str(softclip_ratio)+' > '+each+'.softclip_realigned.bam')

		print "Freebayes start variant calling..."
		os.system('freebayes -F '+str(freebayes_f)+' -f '+reference['freebayes']+' '+each+'.softclip_realigned.bam | vcfallelicprimitives -k -g | vt normalize -r '+reference['freebayes']+' - > '+each+'.raw.vcf')
		os.system('vcffilter -f "( LEN > '+str(indel_len-1)+' & TYPE = ins & DP > '+str(depth-1)+' ) | ( LEN > '+str(indel_len-1)+' & TYPE = del & DP > '+str(depth-1)+' ) | ( LEN > '+str(indel_len-1)+' & TYPE = complex & DP > '+str(depth-1)+' )" '+each+'.raw.vcf >'+each+'.indel.vcf')
		os.system('vcffilter -f "! ( TYPE = ins ) & ! ( TYPE = del ) & ! ( TYPE = complex ) & DP > '+str(depth-1)+'" '+each+'.raw.vcf >'+each+'.snp.vcf')
		os.remove(each+'.raw.vcf')

		# filtering variants in regions listed in a BED file if provided by -t parameter.
		if bedfile:
			os.system('intersectBed -a '+each+'.indel.vcf -b '+bedfile+' -wa -header -u >'+each+'.indel.filter.vcf')
			os.system('intersectBed -a '+each+'.snp.vcf -b '+bedfile+' -wa -header -u >'+each+'.snp.filter.vcf')
			os.system('mv '+each+'.indel.filter.vcf '+each+'.indel.vcf')
			os.system('mv '+each+'.snp.filter.vcf '+each+'.snp.vcf')

	print "Cleanup sever"
	os.system('gfServer stop localhost 50000')
	print "ScanIndel running done: "+strftime("%Y-%m-%d %H:%M:%S")
	end = time()
	print 'Analyzing your data takes: '+str(end-start)+' seconds.'

if __name__ == '__main__':
	main()
