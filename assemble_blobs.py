#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Assemble raw sequencing reads into blob contig assembly.

Each contig ID is generated from the CRC64 checksum hex digest of its 
corresponding contig sequence. Using this method, contigs with identical 
sequence will have the same contig ID, and duplicate sequences can be
identified more efficiently.
'''

from Bio import SeqIO
import os
import re
from shutil import copyfile
from subprocess import check_call
import sys

from ArgUtil import InputHandler
from ArgUtil import printArgs
import Util
from Util import update

################################################################################

contig_prefix = 'CRC64JONES-'

################################################################################

def assembleBlobs(args):
	'''Assemble reads into blob contigs.'''
	
	assembly_type = args.assembly_type
	
	if assembly_type == 'SGA':
		assembleBlobsSGA(args)
	elif assembly_type == 'unknown':
		assembleBlobsSGA(args) 
	else:
		raise ValueError("unsupported assembly type: {!r}".format(assembly_type))

def assembleBlobsSGA(args):
	'''Assemble reads into blob contigs using SGA.'''
	
	# TODO: choose preprocessing index type based on read length
	
	Util.checkSGAVersion()
	
	# Set preprocessing quality threshold.
	quality_trim = '20'
	
	# Set error-correction kmer length.
	kmer_length = '41'
	
	# Get arguments.
	assembly_file = args.assembly_file
	forward_reads = args.forward_reads
	num_threads = str(args.num_threads)
	orphan_reads = args.orphan_reads
	reverse_reads = args.reverse_reads
 
	phred_option = ['--phred64'] if args.phred64 else []
	len_cutoff = min(args.len_cutoffs)
	
	update("Starting SGA contig assembly")
	
	update("Creating temp directory")
	with Util.tempDirectory() as twd:
		update("Created temp directory: {!r}".format(twd))
		
		update("Moving to temp directory")
		with Util.tempMove(twd):

			# Open NULL output handle.
			with open(os.devnull, 'w') as nul:
			
				stem = 'temp'
			
				prep_file = stem + '.fq'
				if reverse_reads is not None:
				
					update("Preprocessing PE reads")
					orphaned_file = 'orphaned.fq'
					check_call(['sga', 'preprocess'] + phred_option + [ '-p', '1', 
						'--pe-orphans=' + orphaned_file, '-q', quality_trim, 
						'-o', prep_file, forward_reads, reverse_reads ], 
						stdout=nul, stderr=nul)
					
					update("Appending orphaned reads to preprocessed PE reads")
					Util.appendTextFile('orphaned.fq', prep_file)
		
				else:
			
					update("Preprocessing SE reads")
					check_call(['sga', 'preprocess'] + phred_option + 
						[ '-q', quality_trim, '-o', prep_file, forward_reads ], 
						stdout=nul, stderr=nul)
				
				if orphan_reads is not None:
			
					update("Preprocessing existing orphan reads")
					orphans_file = 'orphans.fq'
					check_call(['sga', 'preprocess'] + phred_option + 
						[ '-q', quality_trim, '-o', orphans_file, orphan_reads ], 
						stdout=nul, stderr=nul)
					
					update("Appending existing orphan reads to preprocessed reads")
					Util.appendTextFile(orphans_file, prep_file)
			
				update("Read preprocessing complete")
						   
				update("Indexing preprocessed reads")
				# (Outputs index files with extensions: '.bwt' '.sai')
				check_call(['sga', 'index', '-t', num_threads, '-a', 'ropebwt', 
					'--no-reverse', prep_file], stdout=nul, stderr=nul) 

				update("Error-correcting preprocessed reads")
				ec_file = stem + '.ec.fq'
				check_call(['sga', 'correct', '-t', num_threads, '-k', kmer_length, 
					prep_file, '-o', ec_file], stdout=nul, stderr=nul)
			
				update("Indexing error-corrected reads")
				# (Outputs files with extensions: '.bwt' '.sai' '.rbwt' '.rsai')
				check_call(['sga', 'index', '-t', num_threads, '-a', 'ropebwt', 
					ec_file], stdout=nul, stderr=nul)

				update("Removing duplicates")
				nr_file = stem + '.nr.fa'
				check_call(['sga', 'rmdup', '-t', num_threads, ec_file, 
					'-o', nr_file], stdout=nul, stderr=nul)
			
				update("Finding overlaps")
				asqg_file = stem + '.nr.asqg.gz'
				check_call(['sga', 'overlap', '-t', num_threads, nr_file], 
					stdout=nul, stderr=nul)
			
				update("Assembling reads")
				contig_file = stem + '-contigs.fa'
				variant_file = stem + '-variants.fa'
				check_call(['sga', 'assemble', asqg_file, '-o', stem], 
					stdout=nul, stderr=nul)
		
				update("Combining SGA output into one assembly file")
				# Filtering by minimum sequence length and setting contig IDs.
				cat_file = stem + '.cat.fa'
				with open(cat_file, 'w') as fout:
					for seq_file in (contig_file, variant_file):
						with open(seq_file, 'r') as fin:
							for record in SeqIO.parse(fin, 'fasta'):
								if len(record.seq) >= len_cutoff:
									record.id = genContigID(record)
									record.name = record.description = ''
									SeqIO.write(record, fout, 'fasta')
				
				update("Filtering duplicate contigs")
				# Filtering by duplicate contig IDs and sequences.
				uniq_file = stem + '.uniq.fa'
				filterDuplicateContigs(cat_file, uniq_file)
			
			update("Writing SGA contig assembly to {!r}".format(assembly_file))
			copyfile(uniq_file, assembly_file)

	update("Completed SGA contig assembly")

def filterDuplicateContigs(input_file, output_file):
	'''Filter duplicate contigs in a suitable-formatted FASTA file.
	
	This function requires that contig IDs have the prefix 'BLOBSUM-', followed 
	by the CRC64 hex digest of the contig sequence. The CRC function used is the
	predefined function 'crc-64-jones' implemented in crcmod, as described by 
	Jones (2002).
	
	References
		
	Bairoch, Apweiler (2000) The SWISS-PROT protein sequence database and its 
	supplement TrEMBL in 2000. Nucleic Acids Research. 28(1):45-8. [PMID:10592178]
	
	Jones, David (2002). An improved 64-bit cyclic redundancy check for protein 
	sequences. University College London. 
	
	Press, Flannery, Teukolsky, Vetterling (1993) Cyclic redundancy and other 
	checksums. Numerical recipes in C (Second Edition), New York: Cambridge 
	University Press. [ISBN:0-521-43108-5]
	'''
	
	pattern = re.compile("^{}[A-F0-9]+$".format(contig_prefix))
	
	# Get contig ID counts.
	input_ids = dict()
	with open(input_file, 'r') as fin:
		for record in SeqIO.parse(fin, 'fasta'):
		
			contig_id = record.id
			
			if pattern.match(contig_id) is None:
				raise ValueError("unknown contig ID: {}".format(contig_id))
			
			input_ids.setdefault(contig_id, 0)
			input_ids[contig_id] += 1

	# Compile set of duplicate contig IDs.
	duplicate_ids = set([ k for k in input_ids if input_ids[k] > 1 ])
	
	# Init duplicate sequences.
	duplicate_seqs = dict()
	
	# Init set of output contig IDs.
	output_ids = set()
		
	# If duplicate IDs found, filter duplicates..
	if len(duplicate_ids) > 0:
		
		with open(output_file, 'w') as fout:
			with open(input_file, 'r') as fin:
				for record in SeqIO.parse(fin, 'fasta'):
					contig_id = record.id
					if contig_id in duplicate_ids:
						if contig_id in output_ids:
							if record.seq != duplicate_seqs[contig_id]:
								raise ValueError("contig ID hash collision (crc-64-jones): {}".format(contig_id))
							continue
						duplicate_seqs[contig_id] = record.seq
					SeqIO.write(record, fout, 'fasta')
					output_ids.add(contig_id)

	# ..otherwise simply copy input file to output file.
	else:
	
		copyfile(input_file, output_file)

def genContigID(record):
	'''Generate a contig ID from the CRC64 checksum hex digest of the contig sequence.'''
	return '{}{}'.format( contig_prefix, Util.stringCRC64(str(record.seq)) )

################################################################################
supported_params = ('assembly options', 'input reads', 'num_threads', 
	'Phred options')
required_params = ('assembly', 'forward_reads')

def main(argv):

	Util.checkPythonVersion()
	
	params = InputHandler(supported=supported_params, required=required_params)
	
	argparser = params.initArgparser()
	
	argparser.description = __doc__
	
	args = argparser.parse_args(argv[1:])	
	
	update("Processing arguments")
	
	args = params.processArgs(args)	
	
	printArgs(args)
	
	assembleBlobs(args)   
	
	update("Done.")

################################################################################
	
if __name__ == '__main__':
	main(sys.argv)

################################################################################
