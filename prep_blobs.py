#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Prepare sequencing read data for the blobtools-light pipeline.

If a BLAST database directory is given, this is checked to determine if the NCBI 
'nt' database is present in the given directory. If not found, the 'nt' database 
is downloaded from the NCBI site to the given directory.

A blob contig assembly is created if both read files and an assembly file are 
specified, and if the specified assembly file does not exist or is older than 
the given read files.

Mapping of reads to blob contig assembly is attempted if read files, an assembly 
file, and a mapping file are specified.

If a BLAST result file is specified, the blob contig file is screened using 
BLAST against the NCBI 'nt' database.
'''

import os
import sys

from ArgUtil import InputHandler
from ArgUtil import printArgs
from NcbiUtil import NcbiBlastdb
import Util
from Util import update

from assemble_blobs import assembleBlobs
from blast_blobs import blastScreen
from map_blobs import mapBlobs

from assemble_blobs import supported_params as assembly_params
from blast_blobs import supported_params as blast_params
from map_blobs import supported_params as mapping_params

################################################################################

def prepBlobs(args):
	'''Prepare reads for the blobtools-light pipeline.'''
	
	# Get arguments.
	assembly_file = args.assembly_file
	blast_file = args.blast_file
	blastdb = args.blastdb
	forward_reads = args.forward_reads
	input_reads = args.input_reads
	main_input_reads = args.main_input_reads
	mapping = args.mapping_file
	orphan_reads = args.orphan_reads
	reverse_reads = args.reverse_reads

	# Verify that every specified read file exists.
	if len( filter(os.path.exists, input_reads) ) != len(input_reads):
		raise ValueError("one or more input read files not found") 
	
	update("Starting blobtools-light prep pipeline")
	
	if len(input_reads) > 0:
		update("Validating input read sequence IDs")
		Util.validateFastqSeqids(main_input_reads, orphan_reads=orphan_reads)
	
	if blastdb is not None:
		update("Checking NCBI BLAST 'nt' database in directory: {!r}".format(blastdb))
		if not NcbiBlastdb.present(blastdb, 'nt'):
			update("Downloading NCBI BLAST 'nt' database")
			NcbiBlastdb('nt').download(blastdb)
	
	# If blob contig assembly specified, perform actions to create and use it.
	if assembly_file is not None:
	
		# If reads specified and blob contig assembly either doesn't
		# exist or is older than reads, create blob contig assembly.
		if len(input_reads) > 0 and ( not os.path.exists(assembly_file) or 
			max( map(os.path.getmtime, input_reads) ) > os.path.getmtime(assembly_file) ):
			assembleBlobs(args)

		# If mapping specified, map reads to blob contig assembly.
		if mapping is not None:
			if len(input_reads) == 0:
				raise ValueError("cannot map reads without read data")
			mapBlobs(args)
		
		# If BLAST file specified, run BLAST search of blob 
		# contig assembly against the NCBI 'nt' database.
		if blast_file is not None:
			if not os.path.exists(assembly_file):
				raise ValueError("cannot run BLAST without assembly file")
			if blastdb is None:
				raise ValueError("cannot run BLAST without a database")
			blastScreen(args)

	update("Completed blobtools-light prep pipeline")

################################################################################

supported_params = tuple(set(assembly_params + blast_params + mapping_params))
required_params = ()

def main(argv):

	Util.checkPythonVersion()
	
	params = InputHandler(supported=supported_params, required=required_params)
	
	argparser = params.initArgparser()
	
	argparser.description = __doc__
	
	args = argparser.parse_args(argv[1:])	
	
	update("Processing arguments")
	
	args = params.processArgs(args)	
	
	printArgs(args)
	
	prepBlobs(args)  
	
	update("Done.") 
	
################################################################################
	
if __name__ == '__main__':
	main(sys.argv)

################################################################################
