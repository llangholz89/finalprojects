#! usr/bin/env python

#####################################################
'''Creates tables for each chromosome that contain:
1) Chromosome name 2) gene name 3) gene length 4) GC%'''
#####################################################

from __future__ import absolute_import, division, print_function
from collections import Counter
import pybedtools
import gzip
import ipdb

#---------------------
# hg19 genes bed file:
test_file='/Users/loganlangholz/Documents/class_files/Graduate/Spring_2016/Genome_analysis/Projects/data-sets/bed/test.genes.hg19.bed.gz'
filename='/Users/loganlangholz/Documents/class_files/Graduate/Spring_2016/Genome_analysis/Projects/data-sets/bed/genes.hg19.bed.gz'

test_genes = pybedtools.BedTool(test_file)
genes = pybedtools.BedTool(filename)

#-----------------------------------------------------

# fxn to get GC content from chr and interval of record in bed file:
def GC_content(chrom, start, end):

    fasta_file = '/Volumes/My Book/Genome_Analysis-Final_Project/%s.fa' %chrom
    fasta = pybedtools.BedTool(fasta_file)

    bed = pybedtools.BedTool('%s\t%s\t%s' %(chrom, start, end), \
    from_string = True)

    nuc_content = bed.nucleotide_content(fasta)
    GC = float(nuc_content[0][4]) #GC% is index4 of the 1st and only feature

    return GC

#-------------------------------------------------------
#Create tsv file with Chrom, Gene Name, Gene Length:

fw = open('tbl_chr-genes.tsv', 'w')

for record in genes:

    chrom = str(record.chrom)

    length = record.length
    name = str(record.name)
    start = int(record.start)
    end = int(record.end)
    GC_percent = GC_content(chrom, record.start, record.end)

    fw.write('%s\t%d\t%d\t%s\t%d\t%f\n' %(chrom, start, end, name, length, GC_percent))

fw.close()
#---------------------
#Create a bed file with chr, start, end, name, and GC content

### Note: hg19.chr*.fa files are located on external hard drive ###
'''
for i in range(0,len(chroms)):
    fasta = '/Volumes/My Book/Genome_Analysis-Final_Project/%s.fa.gz' % chroms[i]
    import ipdb; ipdb.set_trace()
    #fa = genes.sequence(fi = fasta) #this will run the whole file,only
    #want it by chromosome...

    #fw = open('%s.tsv' %chroms[i], 'w')
'''




