#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''BLAST blob contig assembly against NCBI 'nt' database.'''

from Bio.Blast.Applications import NcbiblastnCommandline
import os
from shutil import copyfile
import sys
from tempfile import NamedTemporaryFile

from ArgUtil import InputHandler
from ArgUtil import printArgs
from NcbiUtil import filterBestBlastHits
from NcbiUtil import NcbiBlastdb
from NcbiUtil import NcbiGilist
from NcbiUtil import validateBlastFields
import Util
from Util import update

################################################################################

def blastScreen(args):
	'''BLAST sequences against the NCBI 'nt' database.'''

	Util.checkBlastnVersion()

	# Get arguments.
	assembly_file = args.assembly_file
	blast_fields = args.blast_fields
	blast_file = args.blast_file
	blast_task = args.blast_task
	blastdb = args.blastdb
	num_threads = str(args.num_threads)
	
	validateBlastFields(blast_fields)
	
	# Set BLAST format string.
	blast_format = "'6 {}'".format(' '.join(blast_fields))
	
	# Set BLAST database path.
	db = os.path.join(blastdb, 'nt')
	
	update("Starting BLAST screen of NCBI 'nt' database")
	
	update("Checking NCBI BLAST 'nt' database in directory: {!r}".format(blastdb))
	if not NcbiBlastdb.present(blastdb, 'nt'):
		NcbiBlastdb('nt').download(blastdb)
		
	update("Creating temp directory")
	with Util.tempDirectory() as twd:
		update("Created temp directory: {!r}".format(twd))
		
		full_blast = NamedTemporaryFile(dir=twd, delete=False)
			
		update("Running BLAST search against NCBI 'nt' database")
		blast_cline = NcbiblastnCommandline(task=blast_task, 
			num_threads=num_threads, query=assembly_file, db=db, culling_limit=2,
			outfmt=blast_format, evalue=1e-5, out=full_blast.name)		
		blast_cline(stdout=False, stderr=False)
		
		filt_blast = NamedTemporaryFile(dir=twd, delete=False)
		
		update("Filtering BLAST results to retain hit with best E-value and bitscore for each contig")
		filterBestBlastHits(full_blast.name, filt_blast.name, 
			blast_fields=blast_fields)
		
		update("Writing BLAST screen results to {!r}".format(blast_file))
		copyfile(filt_blast.name, blast_file)

	update("Completed BLAST screen of NCBI 'nt' database")

def blastTaxa(args):
	'''BLAST sequences against specified taxa.'''
	
	Util.checkBlastnVersion()
	
	# Get arguments.
	assembly_file = args.assembly_file
	blast_fields = args.blast_fields
	blast_file = args.blast_file
	blast_task = args.blast_task
	blastdb = args.blastdb
	email = args.email
	num_threads = args.num_threads
	taxa = args.target_taxa
	
	validateBlastFields(blast_fields)
	
	# Set BLAST format string.
	blast_format = "'6 {}'".format(' '.join(blast_fields))
	
	# Set BLAST database path.
	db = os.path.join(blastdb, 'nt')
	
	update("Starting BLAST search of specified taxa in NCBI 'nt' database")

	update("Checking NCBI BLAST 'nt' database in directory: {!r}".format(blastdb))
	if not NcbiBlastdb.present(blastdb, 'nt'):
		NcbiBlastdb('nt').download(blastdb)

	update("Creating temp directory")
	with Util.tempDirectory() as twd:
		update("Created temp directory: {!r}".format(twd))
		
		update("Creating temp GI list for taxa: {}".format(
			', '.join([ str(t) for t in taxa ])))
		gilist_temp = NamedTemporaryFile(dir=twd, suffix='.n.gil', delete=False)
		gilist = NcbiGilist()
		for t in taxa:
			gilist += NcbiGilist.fromTaxid(t, 'nucl', email)
		gilist.write(gilist_temp.name)
		
		full_blast = NamedTemporaryFile(dir=twd, delete=False)
				
		update("Running BLAST search against specified taxa")
		blast_cline = NcbiblastnCommandline(task=blast_task, outfmt=blast_format, 
			db=db, query=assembly_file, gilist=gilist_temp.name, evalue=1e-5, culling_limit=2,
			num_threads=num_threads, out=full_blast.name)
		blast_cline(stdout=False, stderr=False)
		
		filt_blast = NamedTemporaryFile(dir=twd, delete=False)
		
		update("Filtering BLAST results to retain hit with best E-value and bitscore for each contig")
		filterBestBlastHits(full_blast.name, filt_blast.name, 
			blast_fields=blast_fields)
		
		update("Writing BLAST results to {!r}".format(blast_file))
		copyfile(filt_blast.name, blast_file)
		
	update("Completed BLAST search of specified taxa in NCBI 'nt' database")

################################################################################

supported_params = ('assembly', 'BLAST options', 'email', 'num_threads')
required_params = ('assembly', 'blast_file', 'blastdb') 

def main(argv):
   
   	Util.checkPythonVersion()
   
	params = InputHandler(supported=supported_params, required=required_params)
   
	argparser = params.initArgparser()
   
   	argparser.description = __doc__
   
	args = argparser.parse_args(argv[1:])
	
	update("Processing arguments")
	
	args = params.processArgs(args)
	
	printArgs(args)
	
	blastScreen(args)
	
	update("Done.")
	
################################################################################
	
if __name__ == '__main__':
	main(sys.argv)
		
################################################################################
