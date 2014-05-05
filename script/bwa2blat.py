#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: bwa2blat.py
#
#        USAGE: ./bwa2blat.py blatDir input.bam output.bam(optional)
#
#  DESCRIPTION: the program is used to revise the CIGAR string with soft clip in BWA using BLAT mapping 
#
#      OPTIONS: ---
# REQUIREMENTS: psl2sam.pl is used
#         BUGS: ---
#        NOTES: need to set the path for BLAT, psl2sam.pl before running
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: Wed Apr 30 13:27:17 CDT 2014
#     REVISION: ---
#===============================================================================
import sys
import os
import time
import re
import pysam

path = os.path.dirname(os.path.realpath(__file__))
bwa_bam = pysam.Samfile(sys.argv[2],'rb')
if len(sys.argv) > 3:
	blat_bam = pysam.Samfile(sys.argv[3],'wb', template=bwa_bam)
else:
	blat_bam = pysam.Samfile(sys.argv[2]+'.blat.bam', 'wb', template=bwa_bam)

for read in bwa_bam.fetch():
	rlen = read.qlen
	if not read.cigar:
		continue
	for each in read.cigar:
		if each[0] == 4 and each[1]/float(rlen) >= 0.2:

			# generating fasta file
			fa = open('temp.fa','w')
			print >> fa,'>'+read.qname
			print >> fa, read.seq
			fa.close()

			# run blat
			os.system('gfClient localhost 50000 '+sys.argv[1]+' temp.fa temp.psl >/dev/null 2>&1 ')
			os.system('perl '+path+'/psl2sam.pl temp.psl >temp.sam')

			# matching genomic coordinate
			sam = open('temp.sam')
			for line in sam:
				blatitem = line.rstrip().split()
				if blatitem[2] == bwa_bam.getrname(read.tid) and abs(int(blatitem[3] ) - read.pos) < 2:
					read.cigarstring = re.sub('H','S',blatitem[5])
					break
			sam.close()
			break
	blat_bam.write(read)
os.system('rm temp.*')
bwa_bam.close()
blat_bam.close()
