#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Utilities for the blobtools-light pipeline.'''

from __future__ import division
from binascii import hexlify
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from collections import deque
from collections import OrderedDict
from contextlib import contextmanager
from crcmod.predefined import mkPredefinedCrcFun
import csv
from ete2 import NCBITaxa
from gzip import GzipFile
from hashlib import md5
import io
from itertools import izip_longest
from math import log
from matplotlib.colors import cnames
from MiscFunctions import keyWithMaxVal
import os
import pysam
from pysam import Samfile
import re
from shutil import rmtree
from subprocess import PIPE 
from subprocess import Popen
import sys
from tempfile import mkdtemp
from warnings import warn

################################################################################

assembly_types = ('unknown', 'SGA')

crc64 = mkPredefinedCrcFun('crc-64-jones')

filter_types = ('include', 'exclude')
default_filter_type = 'exclude'

mapping_types = ('SAM', 'BAM')

tax_label_types = ('scientific name', 'Taxonomy ID')
default_tax_label_type = 'scientific name'

RANKS = ('species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom')
default_rank = 'phylum'

radix_vitae = 1

pseudo_domains = set([
	12908,  # unclassified sequences
	28384   # other sequences
])

special_taxa = set( list(pseudo_domains) + [
	 2387,  # transposons
	 2673,  # insertion sequences
 	29278,  # vectors
 	36549,  # plasmids
 	81077   # artificial sequences
])

imprecise_sciname_prefixes = set([
	'environmental samples',
	'miscellaneous',
	'unclassified',
	'uncultured',
	'unidentified',
	'unknown',
	'unspecified'
])

special_plot_colors = { 
	'no-hit': cnames['lightgray'], 
	'ambig-hit': cnames['darkseagreen'], 
	'other-taxa': cnames['darkslategray']
}

################################################################################

def appendTextFile(src, dest):
	'''Append the specified source text file to the given destination text file.'''
	with open(src, 'r') as fin:		
		with open(dest, 'a') as fout:
			for line in fin:
				fout.write(line)

def bam2cov(bam_file, assembly_file, cov_file):
	'''Convert BAM file to blobtools-light COV format.'''

	checkSamtoolsVersion()

	# Create a temporary sorted and indexed copy of the input BAM file.
	with NamedTemporaryFile(mode='wb', dir=twd, suffix='.bam', delete=False) as sort_bam:
		with open(sort_bam, 'wb') as fbam:  
			for data in pysam.sort('-O', 'bam', '-T', 'temp', bam_file):
				fbam.write(data)
		pysam.index(sort_bam)

	assembly_contigs = set()
	
	with pysam.Samfile(sort_bam, 'rb') as fbam:

		with open(assembly_file, 'r') as fa:
		
			with open(cov_file, 'w') as fcov:
			
				for record in SeqIO.parse(fa, 'fasta'):
					
					contig_id = record.description
					
					assembly_contigs.add(contig_id)
					
					seq = record.seq.upper()
					
					corrected_length = getCorrectedNucSeqLength(seq)
					
					if corrected_length > 0 and contig_id in fbam.references:
						
						total_depth = 0
						
						for column in fbam.pileup(contig_id):
							
							i = column.reference_pos
							
							if seq[i] in unambiguous_dna.letters:
								total_depth += column.nsegments
							
						cov = total_depth / corrected_length
						
					else:
						
						cov = 0
					
					fcov.write("{}\t{}\n".format(contig_id, cov))
		
		# Verify that every contig in the BAM file exists in the assembly.
		for contig_id in fbam.references:
			if contig_id not in assembly_contigs:
				raise RuntimeError("SAM/BAM sequence header {!r} does not seem to be part "
					"of the assembly in file {!r}".format(contig_id, assembly_file))

def checkPythonVersion():
	'''Check that Python is version 2.7+.'''
	if sys.version_info[0] != 2 or sys.version_info[1] < 7:
		raise RuntimeError("Python version 2.7+ required")
		
def checkSoftwareVersion(prog, command, pattern, minimum_version):
	'''Check that software is minimum required version or higher.'''
	
	if not isinstance(prog, basestring):
		raise TypeError("program name must be of string type")
	
	try:
		assert hasattr(command, '__iter__') 
		assert not isinstance(command, basestring)
		assert all( isinstance(x, basestring) for x in command )
	except AssertionError:
		raise TypeError("program command must be an iterable of strings")
	
	try:
		assert pattern.__class__.__name__ == 'SRE_Pattern'
	except (AssertionError, AttributeError):
		raise TypeError("version pattern must be of type 'SRE_Pattern'")

	try:
		assert isinstance(minimum_version, basestring)
		min_ver = [ int(x) for x in minimum_version.split('.') ]
	except (AssertionError, ValueError):
		raise TypeError("version number must be a period-delimited string of integers")

	try:
		p = Popen(command, stdin=PIPE, stdout=PIPE)
	except (IOError, OSError):
		raise RuntimeError("{} not found".format(prog))	

	out, err = p.communicate()

	m = pattern.search(out)

	try:
		real_ver = [ int(x) for x in m.groups() ]
	except AttributeError:
		raise RuntimeError("{} version output could not be parsed".format(prog))
	except ValueError:
		raise RuntimeError("{} version is not numeric".format(prog))
	
	if len(real_ver) != len(min_ver):
		raise RuntimeError("{} version numbers cannot be compared".format(prog))
	
	ver_diff = [ x2 - x1 for x1, x2 in zip(min_ver, real_ver) ]

	for i in range( len(ver_diff) ):
		if all( ver_diff[j] == 0 for j in range(i) ) and ver_diff[i] < 0:
			raise RuntimeError("{} version {} or higher required".format(prog, minimum_version))

def checkBlastnVersion(minimum_version='2.2.29'):
	'''Check BLASTN meets minimum version requirement.'''
	prog = 'BLASTN'
	command = ["blastn", "-version"]
	pattern = re.compile("blastn: (\d+)[.](\d+)[.](\d+)[+]")
	checkSoftwareVersion(prog, command, pattern, minimum_version)   

def checkBowtie2Version(minimum_version='2.2.5'):
	'''Check Bowtie2 meets minimum version requirement.'''
	prog = 'Bowtie2'
	command = ["bowtie2", "--version"]
	pattern = re.compile("bowtie2-align.+ version (\d+)[.](\d+)[.](\d+)")
	checkSoftwareVersion(prog, command, pattern, minimum_version)

def checkSamtoolsVersion(minimum_version='1.2'):
	'''Check SAMtools meets minimum version requirement.'''
	prog = 'SAMtools'
	command = ["samtools", "--version"]
	pattern = re.compile("samtools (\d+)[.](\d+)")
	checkSoftwareVersion(prog, command, pattern, minimum_version)

def checkSGAVersion(minimum_version='0.10.13'):
	'''Check SGA meets minimum version requirement.'''
	prog = 'SGA'
	command = ["sga", "--version"]
	pattern = re.compile("String Graph Assembler [(]sga[)] Version (\d+)[.](\d+)[.](\d+)")
	checkSoftwareVersion(prog, command, pattern, minimum_version)

def formatBaseCount(n):
	'''Format an integer nucleotide count as a string for display.'''
	
	metric_suffixes = ('nt', 'kb', 'Mb', 'Gb')
	
	try:
		n = int(n)
	except (TypeError, ValueError):
		raise RuntimeError("invalid nucleotide count: {!r}".format(n))
	
	# Set logarithm base.
	base = 1000
	
	# Get exponent of nucleotide count.
	exp = int( log(n, base) )
	
	# Enforce upper bound on exponent.
	exp = min(exp, len(metric_suffixes))
	
	n = n / base ** exp
	
	if exp == 0:
		x = "{:d}{}".format(int(n), metric_suffixes[exp])
	else:
		x = "{:.2f}{}".format(n, metric_suffixes[exp])
		
	return x

def formatTaxon(taxid, tax_label='scientific name', taxonomy=None):
	'''Format NCBI Taxonomy ID as string for display.'''
	
	# Load NCBI Taxonomy data if necessary.
	if tax_label == 'scientific name' and taxonomy is None:
		taxonomy = NCBITaxa()
	
	if taxid is None:
		taxon = 'no-hit'
	elif taxid == radix_vitae:
		taxon = 'ambig-hit'
	else:
		try:
			taxid = int(taxid)
		except (TypeError, ValueError):
			raise RuntimeError("invalid Taxonomy ID: {!r}".format(taxid))

		if tax_label == 'scientific name':
			taxon = taxonomy.translate_to_names([taxid])[0]
			if isinstance(taxon, int):
				raise RuntimeError("scientific name not found - "
					"you may need to update local NCBI Taxonomy database")
		elif tax_label == 'Taxonomy ID':
			taxon = str(taxid)
				
	return taxon

def getConsensusTaxonomyID(taxids, threshold=0.5, taxid_weights=None, taxonomy=None):
	'''Get a consensus Taxonomy ID, given a list of taxids.'''

	# Ensure taxids is a list of integers.
	try:
		taxids = [ int(t) for t in taxids ]
	except (TypeError, ValueError):
		raise ValueError("invalid taxids: {!r}".format(taxids))

	# Check threshold is in valid range.
	min_threshold, max_threshold = 0.5, 1.0
	if threshold < min_threshold or threshold > max_threshold:
		raise ValueError("specified threshold ({!r}) out of range ({}-{})".format(
			threshold, min_threshold, max_threshold))

	# If taxid weights defined, validate these against taxids..
	if taxid_weights is not None:
		if set(taxid_weights) != set(taxids):
			raise ValueError("taxid weight info mismatch")
	# ..otherwise set taxid weights to equal weighting.
	else:
		taxid_weights = { t: 1 for t in taxids }

	# Calculate total weight of taxids.
	total_weight = sum( taxid_weights.values() )

	# Get threshold weight.
	threshold_weight = threshold * total_weight

	# Get taxid with highest weight.
	argmax_taxid = keyWithMaxVal(taxid_weights)

	# If any single taxid has weight greater than threshold, choose that taxid..
	if taxid_weights[argmax_taxid] > threshold_weight:
						
		taxid = argmax_taxid
	
	# ..otherwise resolve taxids using NCBI Taxonomy information.
	else:
	
		# Load NCBI Taxonomy data if necessary.
		if taxonomy is None:
			taxonomy = NCBITaxa()
			
		# Get minimal topology of the specified taxids.
		tree = taxonomy.get_topology(taxids)
		
		# Check that NCBI topology includes all Taxonomy IDs.
		tree_taxa = [ node.taxid for node in tree.traverse() ]
		if any( taxid not in tree_taxa for taxid in taxids ):
			raise RuntimeError("partial taxonomy tree - "
				"you may need to update local NCBI Taxonomy database")
		
		# Create deque of leaf nodes for leaf-to-root weight calculations.
		# Nodes will be pushed on left of deque and popped from the right.
		nodes = deque([ x for x in tree.iter_leaves() ])
						
		# Init cumulative taxid weight.
		cum_weight = dict()
						
		# Init rotated node count.
		rotations = 0

		while len(nodes) > 0:
							
			# Check we haven't rotated through all available nodes.
			if rotations == len(nodes):
				tree_str = tree.write(format=8, format_root_node=True)
				raise RuntimeError("failed to process NCBI Taxonomy tree: {}".format(tree_str))
						
			try:
				# Pop node from right of deque.
				node = nodes.pop()
		
				# Get taxid of this node.
				t = node.taxid
		
				# Init cumulative weight for this taxid.
				w = taxid_weights[t] if t in taxid_weights else 0
			
				# If node has children, add the sum of their 
				# cumulative weights to that of this taxid.
				if not node.is_leaf():
					w += sum( cum_weight[child.taxid]  # KeyError possible
						for child in node.get_children() )
		
				# Set cumulative weight for this taxid.
				cum_weight[t] = w
			
				# Append parent node to left of deque.
				if not node.is_root():
					nodes.appendleft(node.up)
		
				# Reset rotated node count.
				rotations = 0
		
			except KeyError:
				# Cumulative weight has not been calculated
				# for all children of current node, so push
				# node back onto left of deque and try again
				# after processing all other nodes in deque.
				nodes.appendleft(node)
				rotations += 1
						
		# Traverse tree from root, set consensus Taxonomy ID
		# where cumulative weight of tree exceeds threshold.
		for node in tree.traverse("levelorder"):
			if cum_weight[node.taxid] > threshold_weight:
				taxid = node.taxid

	return taxid

def getCorrectedNucSeqLength(seq):
	'''Get corrected length of nucleotide sequence.'''
	return sum( seq.upper().count(x) for x in unambiguous_dna.letters )

def getFastaSeqids(fasta_file):
	'''Get list of sequence IDs from a FASTA file.'''   
	with open(fasta_file, 'r') as fa:
		seqids = [ record.id for record in SeqIO.parse(fa, 'fasta') ]
	return seqids

def fileMD5(filepath):
	'''Get hex digest of MD5 checksum for the given file.'''
	with open(filepath, 'rb') as f:
		data = f.read()
	md5sum = md5(data)
	return md5sum.hexdigest()

def openFastq(filepath, mode='r'):
	'''Open a FASTQ file (optionally GZIP-compressed).'''
					 
	text_modes = ('r', 'w', 'a', 'rU')
	gzip_modes = ('r', 'w', 'a', 'rb', 'wb', 'ab')
	modes = text_modes + gzip_modes
	valid_modes = [ m for i, m in enumerate(modes) if m not in modes[:i] ]
	gzip_magic = '1f8b' 
		
	if mode not in valid_modes:
		raise ValueError("invalid file mode: {!r}".format(mode))
	
	# Assume no GZIP compression/decompression.	   
	gzipping = False
					
	if mode.startswith('r'):
		
		with io.open(filepath, 'rb') as fh:
			sample = fh.peek()

		magic = hexlify(sample[:2])
		method = ord(sample[2:3])
		
		if magic == gzip_magic:
			if method == 8:
				gzipping = True
			else:
				raise ValueError("input compressed with unknown GZIP method")
			
	elif filepath.endswith('.gz'):
		gzipping = True
	
	if gzipping:

		if mode not in gzip_modes:
			raise ValueError("file mode {!r} should not be used for GZIP-compressed content".format(mode))
		fh = GzipFile(filepath, mode=mode)

	else:

		if mode not in text_modes:
			raise ValueError("file mode {!r} should not be used for plain text".format(mode))
		fh = open(filepath, mode=mode)

	return fh	

def readBlobTable(filepath, cov_libs=[]):
	'''Read blobtools-light blob table file.'''
	
	blob_table_fields = ('contig_id', 'length', 'gc', 'cov', 'taxonomy')
	min_cov = 0.1
	
	data = OrderedDict()
	
	with open(filepath, 'r') as fh:

		table = csv.DictReader( filter(lambda line: not line.startswith('#'), fh), 
			fieldnames=blob_table_fields, dialect='excel-tab')	
	
		for row in table:
		
			b = row['contig_id']
			
			if b in data:
				raise ValueError("duplicate blob ID {!r} in file {!r}".format(b, filepath))
			
			cov_dict = dict()
			for x in row['cov'].split(';'):
				k, v = x.split('=')
				cov_dict[k] = float(v)
			
			if len(cov_libs) == 0:
				cov_libs = cov_dict.keys()
			
			sum_cov = 0.0
			for cov_lib in cov_libs:
				try:
					sum_cov += cov_dict[cov_lib]
				except KeyError:
					raise RuntimeError("coverage not found for blob contig {!r} in coverage library {!r}".format(b, cov_lib))   
			
			blast_dict = dict()
			for x in row['taxonomy'].split(';'):
				k, v = x.split('=')
				blast_dict[k] = v.split(':')[0]
			
			data[b] = {
				'contig_id': b,
				'length': int(row['length']),
				'gc': float(row['gc']),
				'cov': max(sum_cov, min_cov),
				'tax': blast_dict['tax']
			}
			
	return data

def readMappingInfo(filepath):
	'''Read mapping info from the specified file.'''

	checkSamtoolsVersion()

	if filepath.upper().endswith('.SAM'):
		mode = 'r'
	elif filepath.upper().endswith('.BAM'):
		mode = 'rb'
	elif filepath.upper().endswith('.CAS'):
		raise NotImplementedError("cannot currently parse CAS format files")
	elif filepath.upper().endswith('.COV'):
		raise RuntimeError("cannot get mapping info from COV format files")
	else:
		raise ValueError("unknown format of mapping file: {!r}".format(filepath))

	reads = set()
	read2ref = dict()
	
	with pysam.Samfile(filepath, mode=mode) as fh:
		
		for record in fh:

			r = record.qname
			reads.add(r)

			if record.is_unmapped:
				continue

			ref = fh.references[ record.rname ]
			
			mapq = record.mapping_quality
			if mapq == 255: mapq = -1
			
			fwd, rev = range(2)
			
			if record.is_paired:
			
				read2ref.setdefault(r, [None, None])
				i = fwd if record.is_read1 else rev
			
			else:
			
				read2ref.setdefault(r, [None])
				i = fwd
				
			if read2ref[r][i] is not None:
				
				prev_mapq = read2ref[r][i].values()[0]
				if prev_mapq == 255: prev_mapq = -1
				
				if mapq > prev_mapq:
					read2ref[r][i] = { ref: mapq }
				elif mapq == prev_mapq:
					read2ref[r][i][ref] = mapq
				else:
					continue
			else:
				read2ref[r][i] = { ref: mapq }
				
		refs = fh.references
		
	return sorted(reads), read2ref, refs

def sam2bam(sam_file, bam_file):
	'''Convert SAM file to BAM format.'''
	checkSamtoolsVersion()
	with open(bam_file, 'wb') as fbam:
		for data in pysam.view('-b', sam_file):
			fbam.write(data)
				   
def sam2cov(sam_file, assembly_file, cov_file):
	'''Convert SAM file to blobtools COV format.'''
	with tempDirectory() as twd:
		with NamedTemporaryFile(mode='wb', dir=twd, suffix='.bam') as fbam:
			sam2bam(sam_file, fbam.name)
			bam2cov(fbam.name, assembly_file, cov_file)

def stringMD5(string):
	'''Get uppercase hex digest of MD5 checksum for the given string.'''
	checksum = md5(string)
	hexdigest = checksum.hexdigest()	
	return hexdigest.upper()

def stringCRC64(string):
	'''Get uppercase hex digest of CRC64 checksum for the given string.
	
	The CRC function used is the predefined function 'crc-64-jones' implemented 
	in crcmod, as described by Jones (2002).
	
	References
		
	Bairoch, Apweiler (2000) The SWISS-PROT protein sequence database and its 
	supplement TrEMBL in 2000. Nucleic Acids Research. 28(1):45-8. [PMID:10592178]
	
	Jones, David (2002). An improved 64-bit cyclic redundancy check for protein 
	sequences. University College London. 
	
	Press, Flannery, Teukolsky, Vetterling (1993) Cyclic redundancy and other 
	checksums. Numerical recipes in C (Second Edition), New York: Cambridge 
	University Press. [ISBN:0-521-43108-5]
	'''
	checksum = crc64(string)
	hexdigest = str(hex(checksum))
	if hexdigest.startswith('0x'):
		hexdigest = hexdigest[len('0x'):]
	hexdigest = hexdigest.rstrip('L')
	return hexdigest.upper()

@contextmanager
def tempDirectory(suffix='', prefix='tmp', name=None, dir=None):
	'''Create temporary directory that will delete itself on exit.'''
	
	# If a temp directory name was specified, ensure it exists..
	if name is not None:
		
		# Verify temp directory name is a valid pathname component.
		if os.path.split(name)[0] != '':
			raise ValueError("temp directory name must be a valid pathname component")
		
		# Set temp directory name.
		twd = name
		
		# Prepend directory if specified.
		if dir is not None:
			os.path.join(dir, twd)
		
		# Get canonical path of temp directory.
		twd = os.path.realpath(twd)
		
		# Ensure temp directory exists.
		if not os.path.exists(twd):
			os.makedirs(twd)
	
	# ..otherwise, create temp directory in usual way.
	else:
		twd = mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
	
	try:
		yield twd
	finally:
		try:
			rmtree(twd)
		except OSError:
			warn("failed to remove temp directory: {!r}".format(twd), 
				RuntimeWarning)

@contextmanager
def tempMove(directory):
	'''Temporarily move from the current working directory to the specified path.'''
	origin = os.getcwd()
	os.chdir( os.path.realpath(directory) )
	try:
		yield
	finally:
		try:
			os.chdir(origin)
		except OSError:
			raise RuntimeError("cannot return to original directory: {!r}".format(origin))

def update(msg):
	'''Print formatted status message to standard output'''
	print("[STATUS] {}".format(msg))

def validateFastqSeqids(read_files, **kwargs):
	'''Validate FASTQ sequence IDs.'''
	
	orphan_reads = kwargs.pop('orphan_reads', None)
		
	if len(read_files) == 0:
		raise ValueError("no FASTQ read files specified") 
	elif len(read_files) > 2:
		raise ValueError("too many FASTQ read files specified")
	
	seqids = dict()
	
	fhs = [ open(read_file, 'r') for read_file in read_files ]
	fastq_iters = [ SeqIO.parse(fh, 'fastq') for fh in fhs ]
  
	for records in izip_longest( *fastq_iters ):

		if any([ record is None for record in records ]):
			raise RuntimeError("record count mismatch in paired input read files")
		
		ids = [ record.id for record in records ]
		
		if not all( i == ids[0] for i in ids[1:] ):
			raise RuntimeError("record ID mismatch in paired input read files")
			
		i = ids[0]
		
		try:
			seqids[i] += 1
		except KeyError:
			seqids[i] = 1
		
	if orphan_reads is not None:
	
		with open(orphan_reads, 'r') as fh:
			
			for record in SeqIO.parse(fh, 'fastq'):
				
				try:
					seqids[ record.id ] += 1
				except KeyError:
					seqids[ record.id ] = 1

	if any( seqids[i] > 1 for i in seqids ):
		raise RuntimeError("duplicate FASTQ headers in input read files")  

################################################################################

class text_tab(csv.Dialect):
	'''Python CSV dialect for a text tab-delimited file.'''
	delimiter = '\t'
	doublequote = True
	lineterminator = '\n'
	quotechar = '"'
	quoting = csv.QUOTE_MINIMAL
	skipinitialspace = False
	strict = False
csv.register_dialect("text-tab", text_tab)

class blast_tab(text_tab):
	'''Python CSV dialect for a BLAST tab-delimited file.'''
	skipinitialspace = True
	strict = True	
csv.register_dialect("blast-tab", blast_tab)	


################################################################################
