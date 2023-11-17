#!/usr/bin/env python3
##add_3prime_buffer_to_gff.py
##3-12-20
##Groves Dixon

#import modules
import argparse
from sys import exit
import pandas as pd

##################################
############ FUNCTIONS ###########
##################################

def modify_gff(input_file, output_file, buffer_length):
  """Function to read in and modify the gff by adding
  buffer_length to the 3' end based on strandedness, then output the
  modified gff"""
  input_file = 'Amil.coding.gff3'
  output_file = 'Amil.coding_genes_3primePlus700.gff3'
  buffer_length = 700
  gff = pd.read_csv(input_file, sep='\t', names = ['chr', 'ref', 'feature', 'start', 'end', 'x', 'strand', 'xx', 'description'])
  gene_sub = gff.loc[gff['feature']=='gene']
  forward = gene_sub['strand']=='+'
  #make modifications for forward strand genes
  gene_sub.loc[forward, 'modified_start'] = gene_sub['start']
  gene_sub.loc[forward, 'modified_end'] = gene_sub['end'] + buffer_length
  #make modifications for reverse strand genes
  gene_sub.loc[~forward, 'modified_start'] = gene_sub['start'] - buffer_length
  gene_sub.loc[~forward, 'modified_end'] = gene_sub['end']
  gene_sub.loc[gene_sub['modified_start'] < 0, 'modified_start'] = 0 #ensure we didn't end up with negative numbers

  #replace columns
  mod_gff = gene_sub.copy()
  mod_gff.loc[:, 'start'] = mod_gff['modified_start'].astype(int)
  mod_gff.loc[:, 'end'] = mod_gff['modified_end'].astype(int)
  mod_gff = mod_gff.drop(['modified_start', 'modified_end'], 1)

  #output
  mod_gff.to_csv(output_file, sep='\t', index=False, header = False)


##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
  ##SET UP ARGUMENT PARSING
  Description = '''
  Description:
  Modify a GFF to add a buffer to the 3' end of each gene. This can serve as extra reassurance you capture
  RNA-seq reads from the 3' UTR of transcripts (especially important for tag-seq).
  '''

  parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
  parser.add_argument('-i', required = True, dest = 'table', help = 'two column table with old seq names thrn new ones')
  parser.add_argument('-o', required = True, dest = 'output_file', help = 'The the output file name')
  parser.add_argument('-b', required = True, dest = 'buffer_length', help = 'The length of the buffer to add in bp')
  

  #--- PARSE ARGUMENTS ---#
  args = parser.parse_args()
  input_file = args.table
  output_file = args.output_file
  buffer_length = args.buffer_length

  #----- RUN FUNCTIONS -----#
  modify_gff(input_file, output_file, buffer_length)
  
  
  
  
