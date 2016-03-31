#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Dustin E Schones and Keji Zhao
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010


import pysam
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import Utility

grep = "grep";
cat = "cat";
#grep = "/bin/grep";
#cat = "/bin/cat";

plus = re.compile("\+");
minus = re.compile("\-");

def separateByChromBamToBed(chroms, bam, extension, requiredFlag= 0, filterFlag= 0, mapq= 0):
    """
    """
    chromFileDict= {}
    for chrom in chroms:
        ## Prepare output files
        tmpFile = chrom + extension
        chromFileDict[chrom]= open(tmpFile, 'w')

    inBam= pysam.AlignmentFile(bam)
    for aln in inBam:
        
        if aln.mapping_quality < mapq:
            continue
        if (aln.flag & requiredFlag) != requiredFlag:
            continue
        if (aln.flag & filterFlag) != 0:
            continue
        if aln.reference_name not in chromFileDict:
            continue # This can happen if some chroms have been removed from GenomeData
        # Segregate reads
        chrom=  aln.reference_name
        if aln.is_reverse:
            strand= '-'
        else:
            strand= '+'
        bedline= '\t'.join([chrom, str(aln.reference_start), str(aln.reference_end + 1), aln.query_name, str(aln.mapping_quality), strand])
        chromFileDict[chrom].write(bedline + '\n')
    inBam.close()
    
    for chrom in chromFileDict:
        chromFileDict[chrom].close()


def combineAllGraphFilesBedToBam(chroms, extension, template_bam, final_out):
    """
    Combine the seperately processed chromosomes, return the output file name
    NB: Some attributes of AlignedSegment cannot be set. 
    See also http://pysam.readthedocs.org/en/latest/usage.html#creating-bam-cram-sam-files-from-scratch
    """
    template= pysam.AlignmentFile(template_bam);
    outfile = pysam.AlignmentFile(final_out, 'wb', template= template);
    template.close()
    for chrom in chroms:
        inbed= open(chrom + extension)
        tid= outfile.get_tid(chrom)
        for line in inbed:
            line= line.strip().split('\t')
            aln= pysam.AlignedSegment() 
            aln.query_name= line[3]
            aln.reference_id= tid
            aln.reference_start= int(line[1])
            aln.cigarstring= str(int(line[2]) - int(line[1]) - 1) + 'M'
            aln.is_reverse= False
            if line[5] == '-':
                aln.is_reverse= True
            outfile.write(aln)
        inbed.close()            
    outfile.close();
    return final_out

def cleanup(chroms, extension):
    for chrom in chroms:
        file = chrom + extension;
        try:
            if os.remove('%s' %
                         (file)): raise
        except: 
            sys.stderr.write("")
#            sys.stderr.write("clean up failed\n");

def combineAllGraphFiles(chroms, extension, final_out):
    """
    Combine the seperately processed chromosomes, return the output file name
    TODO: Do not use a system call to concatenate
    """
    outfile = open(final_out,'w');
    outfile.close();
    
    for chrom in chroms:
        file = chrom + extension;
        if Utility.fileExists(file):
            try:
                if os.system('%s %s >> %s' %
                    (cat, file, final_out)): raise
            except: 
                sys.stderr.write( "")
    #            sys.stderr.write("cat failed\n")
        else:
            print file, " file does not exist."
    return final_out

# DEPRECTED
# =========
#def separateByChrom(chroms, bedfile, extension):
#    """This is a faster version of the original function to split reads by chromosome
#    """
#    chromFileDict= {}
#    for chrom in chroms:
#        ## Prepare output files
#        tmpFile = chrom + extension
#        chromFileDict[chrom]= open(tmpFile, 'w')
#
#    fin= open(bedfile)    
#    for line in fin:
#        # Segregate reads
#        chrom= line.strip().split('\t')[0]
#        chromFileDict[chrom].write(line)
#
#    for chrom in chromFileDict:
#        chromFileDict[chrom].close()
    
#def separateByChrom_deprecated(chroms, file, extension):
#    """
#    It is ok if the chroms do not include all those existing in the file.
#    """
#    for chrom in chroms:
#        match = chrom + "[[:space:]]";
#        tmpFile = chrom + extension;
#        cmd= '%s %s %s > %s' %(grep, match, file, tmpFile);
#        print cmd
#        try:
#            if os.system(cmd): raise
#        except:
#            sys.stderr.write( "Warning: " + str(chrom) + " reads do not exist in " + str(file) + "\n");

#def separateByChromBam(chroms, bam):
#    """This is a faster version of the original function to split reads by chromosome
#    """
#    inBam= pysam.AlignmentFile(bam)
#
#    chromFileDict= {}
#    for chrom in chroms:
#        ## Prepare output files
#        tmpFile = chrom + '.bam'
#        chromFileDict[chrom]= pysam.AlignmentFile(tmpFile, 'wb', template= inBam)
#
#    for aln in inBam:
#        # Segregate reads
#        chrom=  aln.reference_name
#        chromFileDict[chrom].write(aln)
#    inBam.close()
#    
#    for chrom in chromFileDict:
#        chromFileDict[chrom].close()


