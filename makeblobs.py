#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   		: makeblobs.py
Version 	: 0.1
Author 		: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs 		: ?
To do 		: Add span in stats file, remove double entry in statsfile if only one blast file  
"""

from __future__ import division
import sys 
import argparse
import os
import BlobCollection

from ArgUtil import integer_cline_list
from ArgUtil import param_info

class InputObject():
	def __init__(self, args):
		self.assembly_type, self.assembly_file = self.getAssembly(args)
		self.exclude_assembly_cov = args.exclude_assembly_cov
		self.mapping_files = self.getDictOfMappingFiles(args)
		self.blast_files, self.blast_order = self.getDictOfBLASTFiles(args)
		self.outfiles = self.getDictOfOutfiles(args.o)
		self.tax_label = args.tax_label
		self.target_taxa = args.target_taxa
		self.taxrule = args.taxrule 
		self.rank = args.rank

	def printParameters(self):
		print "\n\tAssembly file :\n\t\t- {} (type = '{}')".format(self.assembly_file, self.assembly_type)
		print "\tMapping file(s) :"
		for lib_name, mapping_file in self.mapping_files.items():
			print "\t\t- {} : {}".format(lib_name, mapping_file)
		if not self.assembly_type == 'unknown':
			if self.exclude_assembly_cov:
				print "\t\t- Exclude assembly coverage : {} ".format(self.assembly_file)
			else:
				print "\t\t- Include assembly coverage : {} ".format(self.assembly_file)
		print "\tBLAST file(s) :"
		for blast_name, blast_file in self.blast_files.items():
			print "\t\t- {} : {}".format(blast_name, blast_file)
		if len(self.blast_order) > 1:
			print "\tTaxonomification rule : \"{}\"".format(self.taxrule)
			if self.taxrule == 'A':
				print "\t\t - Taxid with maximal sum of scores is chosen"
			elif self.taxrule == 'B':
				print "\t\t - Confidence in taxid decreases with order in BLAST libs :\n\t\t {}".format(str(self.blast_order))
			else:
				pass
		if self.rank is not None:
			print "\tWhere possible, taxonomy will be inferred for the following taxonomic rank:\n\t\t - \"{}\"".format(self.rank.capitalize())
		print('\tTaxon label type: "{}"'.format(self.tax_label))
		if self.target_taxa is not None:
			print('\tTarget taxa: {}'.format(', '.join([ str(t) for t in self.target_taxa ])))
		print "\tOutfiles :\n\t\t - Blob-file : {}\n\t\t - Stats-file : {}\n".format(self.outfiles['blobs'], self.outfiles['stats']) 


	def getAssembly(self, args):
		assembly = {'unknown' : args.a, 
				'spades' : args.spades,
				'velvet' : args.velvet,
				'abyss' : args.abyss
				}
		assembly = {k: v for k, v in assembly.items() if v}	
		try:
			# get assembly type, assembly file and check whether file exists
			[(assembly_type, assembly_file)] = assembly.items()
			if os.path.exists(assembly_file):
				pass
			else:
				sys.exit("[ERROR] : Assembly file {} does not exist.".format(assembly_file))		
		except ValueError, e:
			# throw ValueError if there are too many (>1) elements to unpack from assembly   
			sys.exit("[ERROR] : Please specify an assembly file.")	
		if assembly_type == 'unknown':
			self.exclude_assembly_cov = True
		return assembly_type, assembly_file

	def getDictOfMappingFiles(self, args):
		mappings = {}
		files = {'CAS' : args.cas, 'BAM' : args.bam, 'SAM' : args.sam, 'COV' : args.cov}
		files = {k: v for k, v in files.items() if v}	
		count = {}
		set_of_files = set()
		for mapping_type, mapping_files in files.items():
			for mapping_file in mapping_files:
				if os.path.exists(mapping_file):
					if mapping_file in set_of_files:
						sys.exit("[ERROR] : Mapping file {} listed more than once.".format(mapping_file))
					set_of_files.add(mapping_file)
					count[mapping_type] = count.get(mapping_type, 0) + 1
					mappings[mapping_type + "_" + str(count[mapping_type])] = mapping_file
				else:
					sys.exit("[ERROR] : Mapping file %s does not exist." %mapping_file)	
		if not mappings and (self.exclude_assembly_cov or self.assembly_type == 'unknown'):
			sys.exit("[ERROR] : Please specify at least one mapping file.")	
		return mappings

	def getDictOfBLASTFiles(self, args):
		blasts = {}
		order = []  
		blast_files = args.blast
		set_of_files = set()
		blast_count = 0
		for blast_file in blast_files:
			if os.path.exists(blast_file):
				if blast_file in set_of_files:
					sys.exit("[ERROR] : BLAST file %s listed more than once." %blast_file)		
				set_of_files.add(blast_file)
				blast_count += 1
				blast_lib = "BLAST_" + str(blast_count)
				blasts[blast_lib] = blast_file
				order.append(blast_lib)
			else:
				sys.exit("[ERROR] : BLAST file %s does not exist." %blast_file)	
		if not blasts:
			sys.exit("[ERROR] : Please specify at least one BLAST file.")	
		return blasts, order

	def getDictOfOutfiles(self, outprefix):
		outfiles = {}
		if outprefix:
			outfiles['blobs'] = outprefix + ".blobplot.txt"
			outfiles['stats'] = outprefix + ".blobplot.stats.txt"
		else:
			outfiles['blobs'] = "blobplot.txt"
			outfiles['stats'] = "blobplot.stats.txt"
		return outfiles

def getInput():
	parser = argparse.ArgumentParser(
		prog='makeblobs.py',
		usage = '%(prog)s -a <ASSEMBLY> -cas <CAS> -blast <BLAST> -o <OUTPUT> [-h]',
		add_help=True)
	# only ONE assembly file
	parser.add_argument('-a', metavar = 'ASSEMBLY_FASTA', default='', help='Assembly file')
	parser.add_argument('-spades', metavar = 'SPADES_FASTA', default='', help='SPADES assembly file')
	parser.add_argument('-velvet', metavar = 'VELVET_FASTA', default='', help='VELVET assembly file')
	parser.add_argument('-abyss', metavar = 'ABYSS_FASTA', default='', help='ABYSS assembly file')
	parser.add_argument('-exclude_assembly_cov', action='store_true' , default=False, help='Exclude coverage from assembly file')
	# multiple mapping files
	parser.add_argument('-cov', metavar = 'COV_FILE', default=[], nargs='+', help='COV (mapping) file')
	parser.add_argument('-bam', metavar = 'BAM_FILE', default=[], nargs='+', help='BAM (mapping) file')
	parser.add_argument('-sam', metavar = 'SAM_FILE', default=[], nargs='+', help='SAM (mapping) file')
	parser.add_argument('-cas', metavar = 'CAS_FILE', default=[], nargs='+', help='CAS (mapping) file')
	# multiple BLAST files
	#parser.add_argument('-tax', action='tax', default=default_rank, choices=RANKS, help='Select target taxonomic rank {!r}. Default = {!r}'.format(RANKS, default_rank))
	parser.add_argument('-blast', metavar = 'BLAST_FILE', default=[], nargs='+', help='BLAST file') 
	parser.add_argument('-rank', default=param_info['rank']['default'], choices=param_info['rank']['choices'], help=param_info['rank']['help']) 
	parser.add_argument('-tax_label', default=param_info['tax_label']['default'], choices=param_info['tax_label']['choices'], help=param_info['tax_label']['help'])
	parser.add_argument('-target_taxa', type=integer_cline_list, default=param_info['target_taxa']['default'], help=param_info['target_taxa']['help'])
	parser.add_argument('-taxrule', metavar = 'A or B', default='A', help='Tax-rule on how to deal with multiple BLAST libs. A : "higher bitscore wins", B : "Decreasing trust in BLAST libs"') 
	parser.add_argument('-o', metavar = 'OUTPUT_PREFIX', default='', help='Output prefix') 
	# Version number
	parser.add_argument('-v', action='version', version='%(prog)s version 0.1')
	args = parser.parse_args()
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	return InputObject(args)

def parseInput(parameters):
	'''
	1. Create a BlobCollection Object with chosen preferred rank
	'''
	data = BlobCollection.BlobCollection(rank_level=parameters.rank, tax_label=parameters.tax_label)
	'''
	2. Save outfile-dict as part of BlobCollection object
	'''
	data.outfiles = parameters.outfiles
	'''
	3. Parse contigs of assembly into Blob objects in BlobCollection object
	'''
	data.getBlobsFromAssembly(parameters.assembly_file, parameters.assembly_type, parameters.exclude_assembly_cov)
	'''
	4a. Parse coverage library files into Blob objects in BlobCollection object
	'''
	data.getCovForBlobs(parameters.mapping_files)
	'''
	4b. Print COVs to files if mappings are BAM/SAM/CAS; to avoid having to parse them again 
	'''
	data.printCOVToFiles(parameters.mapping_files)
	'''
	5. Parse BLAST files into Blob objects in BlobCollection object and group into taxonomic bins by translating taxids to taxonomic group
	'''
	data.getTaxForBlobs(parameters.blast_files, target_taxa=parameters.target_taxa)
	'''
	6. Infer a "taxonomy" from the taxonomic bins based on taxrule for each of the blobs
	'''
	data.getConsensusTaxForBlobs(parameters.taxrule, parameters.blast_order)
	'''
	7. Return BlobCollection object
	'''
	return data

if __name__ == "__main__":
	
	__version__ = 0.1

	# Put all parameters of the run into InputObject
	parameters = getInput()
	
	# Print InputObject so that user gets feedback about whats happening
	parameters.printParameters()
	
	# Parse files according to InputObject into BlobCollection object 
	data = parseInput(parameters)

	# Do stats on BlobCollection object so that the stats file can be printed
	data.getStats()

	# Write Output to the two files *blobplot.txt and *stats.txt
	data.writeOutput(__version__)

