#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Map sequencing reads to blob contig assembly.'''

import inspect
import os
import re
from shutil import copyfile
from subprocess import check_call
import sys

from ArgUtil import InputHandler
from ArgUtil import printArgs
import Util
from Util import mapping_types
from Util import update

################################################################################

def mapBlobs(args):
	'''Map reads to blob contig assembly.'''	
	
	if args.mapping_type not in mapping_types:
		raise ValueError("unsupported mapping type: {!r}".format(mapping_type))
	
	mapBlobsBowtie2(args)

def mapBlobsBowtie2(args):
	'''Map reads to blob contig assembly using Bowtie2.''' 
	
	Util.checkBowtie2Version()
	
	# Get arguments.
	assembly_file = args.assembly_file
	forward_reads = args.forward_reads
	mapping_file = args.mapping_file
	mapping_type = args.mapping_type
	num_threads = str(args.num_threads)
	orphan_reads = args.orphan_reads
	reverse_reads = args.reverse_reads
	
	phred_option = ['--phred64'] if args.phred64 == True else ['--phred33']
		 
	update("Starting Bowtie2 read mapping")

	update("Creating temp directory")
	with Util.tempDirectory() as twd:
		update("Created temp directory: {!r}".format(twd))
		
		update("Moving to temp directory")
		with Util.tempMove(twd):
			
			update("Checking for index files of assembly {!r}".format(assembly_file))
			
			# Get assembly prefix path: assembly file path without extension.
			assembly_prefix, ext = os.path.splitext(assembly_file)
			 
			# Get assembly file stem: assembly file basename without extension.
			assembly_stem = os.path.basename(assembly_prefix)
			
			# Get directory containing assembly file.
			assembly_dir = os.path.dirname(assembly_file)
			
			# Create regexes for forward and reverse index files.
			findex_pattern = re.compile("^{}[.](\d+)[.]bt2$".format(assembly_stem))
			rindex_pattern = re.compile("^{}[.]rev[.](\d+)[.]bt2$".format(assembly_stem))
			
			# Init list of index files and index file numbers.
			findex_nums, rindex_nums = ([], [])
			index_paths = list()
			
			# Get filename of each index file in assembly directory, 
			# as well as its index file number.
			for filename in os.listdir(assembly_dir):
				
				fmatch = findex_pattern.match(filename)
				if fmatch is not None:
					index_paths.append(os.path.join(assembly_dir, filename))
					findex_nums.append(int(fmatch.group(1)))
					continue
					
				rmatch = rindex_pattern.match(filename)  
				if rmatch is not None:
					index_paths.append(os.path.join(assembly_dir, filename))
					rindex_nums.append(int(rmatch.group(1)))
					continue		 
			
			# Assume indexing not necessary.
			indexing = False
			
			# If either forward or reverse index file numbers do not form
			# a consecutive list beginning at 1, indexing is necessary.
			for nums in (findex_nums, rindex_nums):
				nums = sorted(nums)
				try:
					assert nums[0] == 1
					assert all( y - x == 1 for x, y in zip(nums[:-1], nums[1:]) )
				except (AssertionError, IndexError):
					indexing = True
  
			# If any index files are older than the 
			# assembly file, indexing is necessary.
			assembly_mtime = os.path.getmtime(assembly_file)
			index_mtimes = map(os.path.getmtime, index_paths)
			if any( mtime < assembly_mtime for mtime in index_mtimes ):
				indexing = True
  
  			# Open NULL output handle.
  			with open(os.devnull, 'w') as nul:
  			
				if indexing:
					update("Index files not found, creating temp index")
					index_prefix = os.path.basename(assembly_stem)
					check_call(['bowtie2-build', '-q', '--seed', '$RANDOM', 
						assembly_file, index_prefix], stdout=nul, stderr=nul)
				else:
					update("Index files found, using existing index")
					index_prefix = assembly_prefix  
			
				# Map reads to blob contig assembly.
				sam_file = index_prefix + '.sam'
				with open(sam_file, 'w') as fsam:
				
					if reverse_reads is None:
						update("Mapping PE reads to assembly")
						check_call(['bowtie2', '--threads', num_threads, '--seed', 
						'$RANDOM'] + phred_option + ['--very-sensitive-local', 
						'-1', forward_reads, '-2', reverse_reads, '-x', index_prefix], 
						stdout=fsam, stderr=nul)
					else:
						update("Mapping SE reads to assembly")
						check_call(['bowtie2', '--threads', num_threads,  
							'--seed', '$RANDOM'] + phred_option + 
							['--very-sensitive-local', '-U', forward_reads, '-x', 
							index_prefix], stdout=fsam, stderr=nul)
			
				if orphan_reads is not None:
					update("Mapping orphan reads to assembly")
					with open(sam_file, 'a') as fsam:
						check_call(['bowtie2', '--threads', num_threads, 
							'--seed', '$RANDOM'] + phred_option + 
							['--very-sensitive-local', '--no-head', '-U', 
							orphan_reads, '-x', index_prefix], stdout=fsam, 
							stderr=nul)

			update("Writing {} file to {!r}".format(mapping_type, mapping_file))
			if mapping_type == 'SAM':
				copyfile(sam_file, mapping_file)
			elif mapping_type == 'BAM':
				Util.sam2bam(sam_file, mapping_file)

	update("Completed Bowtie2 read mapping")

################################################################################

supported_params = ('alignment options', 'assembly', 'input reads', 
	'num_threads', 'Phred options')
required_params = ('assembly', 'input reads', 'mapping')

def main(argv):
	
	Util.checkPythonVersion()

	params = InputHandler(supported=supported_params, required=required_params)
   
	argparser = params.initArgparser()
	
	argparser.description = __doc__
	
	args = argparser.parse_args(argv[1:])
	
	update("Processing arguments")
	
	args = params.processArgs(args)
	
	printArgs(args) 
	   
	mapBlobs(args)   

	update("Done.")

################################################################################
	
if __name__ == '__main__':
	main(sys.argv)

################################################################################
