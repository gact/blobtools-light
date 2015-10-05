#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Argument processing utilities for the blobtools-light pipeline.'''

import argparse
from collections import OrderedDict
from copy import deepcopy
from ete2 import NCBITaxa
from functools import partial
import os

from NcbiUtil import default_blast_fields
from NcbiUtil import default_blast_task
from NcbiUtil import supported_blast_tasks
from NcbiUtil import validateBlastFields
from Util import default_tax_label_type
from Util import RANKS
from Util import special_taxa
from Util import tax_label_types

# TODO: support for SPADES, VELVET, ABYSS assemblies
# TODO: support for COV, CAS formats

################################################################################

def bounded_integer(argument, argmin=None, argmax=None, argtype=int):
	'''Validate command-line argument as a bounded integer, return as specific type.'''
 
	try:
		# Convert argument to integer.
		value = int(argument)
	except ValueError:
		raise ValueError("argument value is not an integer")

	if argmin is not None:
		
		# Verify argument minimum is of integer type.
		if not isinstance(argmin, int):
			raise TypeError("argument minimum should be of integer type")
			
		# Verify that the integer is not less than minimum bound.
		if value < argmin:
			raise ValueError("value of argument ({!r}) is too low (min={!r})".format(value, argmin))
   
	if argmax is not None:
		
		# Verify argument maximum is of integer type.
		if not isinstance(argmax, int):
			raise TypeError("argument maximum should be of integer type")	
		
		# Verify that the integer is not greater than maximum bound.
		if value > argmax:
			raise ValueError("value of argument ({!r}) is too high (max={!r})".format(value, argmax))

	return argtype(value)

def pathname(argument):
	'''Validate command-line argument as a pathname.'''
	try:
		value = os.path.realpath(argument)
	except AttributeError:
		value = None
	return value

def cline_list(argument, dtype=str):
	'''Validate command-line argument as either a list file or comma-delimited list.'''
	try:
		with open(argument, 'r') as fh:
			values = tuple( dtype( line.rstrip() ) for line in fh )
	except (AttributeError, IOError, OSError, ValueError):
		try:
			values = tuple( map( dtype, argument.split(',') ) )
		except AttributeError:
			raise ValueError("argument is not a list file or comma-delimited list")
		except ValueError:
			raise ValueError("argument is not a list of type {!r}".format(dtype.__name__)) 
	return values 

# Validate argument as a non-negative integer.
non_negative_integer = partial(bounded_integer, argmin=0)

# Validate argument as a positive integer.
positive_integer = partial(bounded_integer, argmin=1)

# Validate argument as a command-line list of integers.
integer_cline_list = partial(cline_list, dtype=int)

################################################################################

# Set parameter groups. All parameters should be in one of these groups; a 
# parameter that does not fit into any group should be added to the general 
# 'keyword options' group. All parameters must be added to 'param_info' below.
param_groups = OrderedDict([
	('input reads', ('forward_reads', 'reverse_reads', 'orphan_reads')), 
	('output reads', ('output_forward_reads', 'output_reverse_reads', 
		'output_orphan_reads')), 
	('assembly options', ('assembly_fasta', 'sga_fasta', 'len_cutoffs')),
	('alignment options', ('sam_file', 'bam_file')),
	('BLAST options', ('blastdb', 'blast_file', 'blast_fields', 'blast_task')),
	('classifier options', ('target_taxa', 'foreign_taxa', 'classifier_ignore')), 
	('filter options', ('include_blobs', 'exclude_blobs', 'exclude_unannotated', 
		'exclude_unmapped')), 
	('Phred options', ('phred33', 'phred64')), 
	('keyword options', ('blob_table', 'blob_stats', 'cov_libs', 'email', 
		'filter_plot', 'filter_table', 'num_threads', 'plot_title', 'tax_label'))
])

# Set mutually exclusive parameter groups. Each mutually exclusive group (MXG) 
# is a group of parameters for which at most one parameter may be specified at 
# a time. Currently, every MXG must be a subset of a parameter group, within a
# hierarchy: parameter group > mutually exclusive parameter group > parameter.
mutually_exclusive_param_groups = OrderedDict([
	('assembly', ('assembly_fasta', 'sga_fasta')),
	('mapping', ('sam_file', 'bam_file')),
	('filter', ('include_blobs', 'exclude_blobs')),
	('Phred offset', ('phred33', 'phred64'))
])

# Creating mappings of: parameter to group, parameter to MXG, MXG to parameter
# group, parameter group to MXG. Also, verify that mappings are unique, and that
# there are no naming conflicts between the parameter groups, MXGs, parameters. 
param2group, param2mxg = OrderedDict(), OrderedDict()
mxg2group, group2mxg = dict(), dict()
for group, params in param_groups.items():
	for p in params:
		assert p not in param_groups, "naming conflict between parameter {!r} and a parameter group".format(p)
		assert p not in mutually_exclusive_param_groups, "naming conflict between parameter {!r} and an MXG".format(p)
		assert p not in param2group, "parameter {!r} matches multiple parameter groups".format(p)
		param2group[p] = group
for mxg in mutually_exclusive_param_groups:
	mxps = mutually_exclusive_param_groups[mxg]
	for p in mxps:
		assert p not in param2mxg, ("parameter {!r} matches multiple MXGs".format(p))
		param2mxg[p] = mxg
	groups = list(set([ param2group[p] for p in mxps ]))
	assert mxg not in param_groups, "naming conflict between MXG {!r} and a parameter group".format(p)
	assert mxg not in mxg2group, ("MXG {!r} matches multiple parameter groups".format(mxg))
	assert groups[0] not in group2mxg, ("parameter group {!r} matches multiple MXGs".format(group))
	mxg2group[mxg] = groups[0]
	group2mxg[ groups[0] ] = mxg

# Set lists of known parameter groups, MXGs, parameters.
known_groups = tuple(param_groups.keys())
known_mxgs = tuple(mutually_exclusive_param_groups.keys())
known_params = tuple(param2group.keys())

# Set info for individual parameters. Each parameter listed here must  
# also be included in 'param_groups', as arguments are added by group.
param_info = {
	'forward_reads': { 
		'flags': ['-1'], 
		'type': pathname,
		'help': 'Input forward read FASTQ file'
	}, 
	'reverse_reads': { 
		'flags': ['-2'], 
		'type': pathname,
		'help': 'Input reverse read FASTQ file'
	}, 
	'orphan_reads': { 
		'flags': ['-0'], 
		'type': pathname,
		'help': 'Input orphan read FASTQ file'
	},	 
	'output_forward_reads': { 
		'flags': ['-f1'], 
		'type': pathname,
		'help': 'Output forward read FASTQ file'
	},
	'output_reverse_reads': { 
		'flags': ['-f2'], 
		'type': pathname,
		'help': 'Output reverse read FASTQ file'
	},
	'output_orphan_reads': { 
		'flags': ['-f0'], 
		'type': pathname,
		'help': 'Output orphan read FASTQ file'
	},
	'assembly_fasta': { 
		'flags': ['-a'], 
		'type': pathname,
		'help': 'Assembly file'
	},
	'sga_fasta': { 
		'flags': ['-sga'], 
		'type': pathname,
		'help': 'SGA assembly file'
	},  
	'len_cutoffs': { 
		'flags': ['-c'], 
		'action': 'append', 
		'default': [],
		'help': 'Minimum blob contig length(s)'
	},
	'sam_file': { 
		'flags': ['-sam'], 
		'type': pathname,
		'help': 'SAM mapping file'
	},
	'bam_file': { 
		'flags': ['-bam'], 
		'type': pathname,
		'help': 'BAM mapping file'
	},
	'blastdb': { 
		'flags': ['-blastdb'], 
		'type': pathname,
		'help': "Path to local copy of NCBI BLAST 'nt' database"
	}, 
	'blast_file': { 
		'flags': ['-blast'], 
		'type': pathname,
		'help': 'BLAST tabular output file'
	},  
	'blast_fields': { 
		'flags': ['-blast_fields'],
		'type': cline_list,
		'default': default_blast_fields,
		'help': 'BLAST fields'
	},  
	'blast_task': { 
		'flags': ['-blast_task'],
		'default': default_blast_task,
		'choices': supported_blast_tasks,
		'help': 'BLAST task (default: {!r})'.format(default_blast_task)	
	},	  
	'target_taxa': { 
		'flags': ['-target_taxa'],
		'type': integer_cline_list,
		'default': (),
		'help': 'NCBI Taxonomy IDs of target taxa'
	},	  
	'foreign_taxa': { 
		'flags': ['-foreign_taxa'], 
		'type': integer_cline_list,
		'default': (),
		'help': 'NCBI Taxonomy IDs of foreign taxa'
	}, 
	'classifier_ignore': { 
		'flags': ['-classifier_ignore'], 
		'type': integer_cline_list,
		'default': (),
		'help': 'List of NCBI Taxonomy IDs to exclude from training set of filter classifier'
	},  
	'include_blobs': { 
		'flags': ['-i', '-include_blobs'], 
		'type': pathname,
		'help': "File listing blob contigs to include, one per line (incompatible with '-e')"
	},  
	'exclude_blobs': { 
		'flags': ['-e', '-exclude_blobs'], 
		'type': pathname,
		'help': "File listing blob contigs to exclude, one per line (incompatible with '-i')"
	},	 
	'exclude_unannotated': { 
		'flags': ['-exclude_unannotated'], 
		'action': 'store_true',
		'default': False,
		'help': 'Exclude blob contigs that do not have taxonomic annotation'
	},	 
	'exclude_unmapped': { 
		'flags': ['-u', '-exclude_unmapped'], 
		'action': 'store_true',
		'default': False,
		'help': 'Exclude unmapped reads'
	},	 
	'phred33': { 
		'flags': ['-phred33'], 
		'action': 'store_true',
		'default': False,
		'help': 'FASTQ files have 33-offset Phred quality scores (default)'
	},  
	'phred64': { 
		'flags': ['-phred64'], 
		'action': 'store_true',
		'default': False,
		'help': 'FASTQ files have 64-offset Phred quality scores'
	},  
	'blob_table': { 
		'flags': ['-blob_table'], 
		'type': pathname,
		'help': 'blobtools-light table file'
	},  
	'blob_stats': { 
		'flags': ['-blob_stats'],
		'type': pathname,
		'help': 'blobtools-light stats file'
	},  
	'cov_libs': { 
		'flags': ['-cov_libs'], 
		'type': cline_list,
		'default': (),
		'help': 'One or more coverage libraries (default: coverage taken from sum of all coverage libraries)'
	},
	'email': { 
		'flags': ['-email'], 
		'help': 'User email (required for any NCBI queries)'
	},
	'filter_table': { 
		'flags': ['-filter_table'], 
		'type': pathname,
		'help': 'Blob contig filter table file'
	},
	'filter_plot': { 
		'flags': ['-filter_plot'], 
		'type': pathname,		
		'help': 'Prefix of blob contig filter plot file'
	},
	'num_threads': { 
		'flags': ['-num_threads'],
		'type': positive_integer,
		'default': 1,
		'help': 'Number of threads for parallel operations'
	},
	'plot_title': {
		'flags': ['-plot_title'],
		'help': 'Title of plot'
	},
	'tax_label': { 
		'flags': ['-tax_label'], 
		'default': default_tax_label_type,
		'choices': tax_label_types,
		'help': "Select taxon label type (default: {!r})".format(default_tax_label_type)
	}	
}

# Verify parameters in 'param_info' match those in 'param_groups'.
assert set(param_info.keys()) == set(known_params), "ArgUtil parameter info mismatch"

################################################################################

def printArgs(args):
	'''Print the specified arguments.'''
	
	kwargs = deepcopy(vars(args))
			
	names = [ k for k in kwargs ]
	values = [ kwargs[k] for k in names ]
	 
	longest_name = max( len(x) for x in names )
	
	template = '{{0: >{}}} : {{1:>}}'.format(longest_name)
	
	print("\n ARGUMENTS \n")
	
	for p in sorted(kwargs):
		
		if hasattr(kwargs[p], '__iter__') and not isinstance(kwargs[p], basestring):
			value = ', '.join( repr(x) for x in kwargs[p] )
		else:
			value = repr(kwargs[p])
			
		print(template.format(p, value))
 
	print('')

def resolveRequiredParameters(param_spec):
	'''Resolve spec of required parameters into parameter groups, MXGs, parameters.'''
	
	groups, mxgs, params = [ set() for _ in range(3) ]

	for x in param_spec:

		if x in known_groups:
			groups.add(x)
			if x in group2mxg:
				mxgs.add(group2mxg[x])
			for p in param_groups[x]:
				params.add(p)
		elif x in known_mxgs:
			mxgs.add(x)
			for p in mutually_exclusive_param_groups[x]:
				params.add(p)
		elif x in known_params:
			if x in param2mxg:
				mxgs.add(param2mxg[x])
			params.add(x)
		else:
			raise ValueError("unknown parameter or group: {!r}".format(x))

	groups = [ g for g in known_groups if g in groups ]
	mxgs = [ mxg for mxg in known_mxgs if mxg in mxgs ]
	params = [ p for p in known_params if p in params ]
	
	return groups, mxgs, params


def resolveSupportedParameters(param_spec):
	'''Resolve spec of supported parameters into parameter groups, MXGs, parameters.'''
	
	groups, mxgs, params = [ set() for _ in range(3) ]

	for x in param_spec:

		if x in known_groups:
			groups.add(x)
			if x in group2mxg:
				mxgs.add(group2mxg[x])
			for p in param_groups[x]:
				params.add(p)
		elif x in known_mxgs:
			if x in mxg2group:
				groups.add(mxg2group[x])
			mxgs.add(x)
			for p in mutually_exclusive_param_groups[x]:
				params.add(p)
		elif x in known_params:
			groups.add(param2group[x])
			if x in param2mxg:
				mxgs.add(param2mxg[x])
			params.add(x)
		else:
			raise ValueError("unknown parameter or group: {!r}".format(x))

	groups = [ g for g in known_groups if g in groups ]
	mxgs = [ mxg for mxg in known_mxgs if mxg in mxgs ]
	params = [ p for p in known_params if p in params ]
	
	return groups, mxgs, params

################################################################################

class InputHandler(object):
	'''Input handler class.'''

	def __init__(self, supported=known_params, required=()):
		'''Init input handler class with supported and required parameters.'''
	
		sup = resolveSupportedParameters(supported)
		req = resolveRequiredParameters(required)
	
		self.supp_groups, self.supp_mxgs, self.supp_params = sup
		self.req_groups, self.req_mxgs, self.req_params = req
	
		for s, r in zip(sup, req):   
			for p in r:
				if p not in s:
					raise ValueError("required parameter {!r} is unknown".format(p))

		self.supp, self.req = dict(), dict()
		for p in known_groups:
			self.supp[p] = p in self.supp_groups
			self.req[p] = p in self.req_groups
		for p in known_mxgs:
			self.supp[p] = p in self.supp_mxgs	
			self.req[p] = p in self.req_mxgs		
		for p in known_params:
			self.supp[p] = p in self.supp_params				
			self.req[p] = p in self.req_params   
		
		if self.supp['reverse_reads'] and not self.supp['forward_reads']:
			raise ValueError("InputHandler cannot support input reverse reads without corresponding forward reads")
		if self.supp['orphan_reads'] and not self.supp['reverse_reads']:
			raise ValueError("InputHandler cannot support input orphan reads without corresponding paired reads")
		if self.supp['output_reverse_reads'] and not self.supp['output_forward_reads']:
			raise ValueError("InputHandler cannot support output reverse reads without corresponding forward reads")
		if self.supp['output_orphan_reads'] and not self.supp['output_reverse_reads']:
			raise ValueError("InputHandler cannot support output orphan reads without corresponding paired reads")			

	def getAssemblyInfo(self, args):
		'''Get assembly type and assembly file from arguments.'''

		if not self.supp['assembly']:
			raise RuntimeError("assembly option not supported")

		assembly_args = [ p for p in mutually_exclusive_param_groups['assembly']
			if hasattr(args, p) and getattr(args, p) is not None ]
		
		# NB: assumes at most one assembly type, since parameter in MXG.
		if len(assembly_args) == 1:
			assembly_arg = assembly_args[0]
		elif not self.req['assembly']:
			assembly_arg = None
		else:
			raise ValueError("no assembly file specified")

		if assembly_arg is not None:
		
			assembly_file = getattr(args, assembly_arg)
			
			if assembly_arg == 'sga_fasta':
				assembly_type = 'SGA'
			elif assembly_arg == 'assembly_fasta':
				assembly_type = 'unknown'
			else:
				raise ValueError("unknown assembly type: {!r}".format(assembly_arg))
			
		else:
			assembly_type = assembly_file = None

		return assembly_type, assembly_file
		
	def getMappingInfo(self, args):
		'''Get mapping type and mapping file from arguments.'''
		
		if not self.supp['mapping']:
			raise RuntimeError("mapping option not supported")
		
		mapping_args = [ p for p in mutually_exclusive_param_groups['mapping']
			if hasattr(args, p) and getattr(args, p) is not None ]

		# NB: assumes at most one mapping type, since parameter in MXG.
		if len(mapping_args) == 1:
			mapping_arg = mapping_args[0]
		elif not self.req['mapping']:
			mapping_arg = None
		else:
			raise ValueError("no mapping file specified")

		if mapping_arg is not None:
		
			mapping_file = getattr(args, mapping_arg)
			
			if mapping_arg == 'sam_file':
				mapping_type = 'SAM'
			elif mapping_arg == 'bam_file':
				mapping_type = 'BAM'
			else:
				raise ValueError("unknown mapping type: {!r}".format(mapping_arg))

		else:
			mapping_type = mapping_file = None

		return mapping_type, mapping_file

	def getFilterInfo(self, args):
		'''Get filter type and filter file from arguments.'''
		
		if not self.supp['filter']:
			raise RuntimeError("filter option not supported")
		
		filter_args = [ p for p in mutually_exclusive_param_groups['filter']
			if hasattr(args, p) and getattr(args, p) is not None ]

		# NB: assumes at most one filter type, since parameter in MXG.
		if len(filter_args) == 1:
			filter_arg = filter_args[0]
		elif not self.req['filter']:
			filter_arg = None
		else:
			raise ValueError("no filter file specified")

		if filter_arg is not None:
		
			filter_file = getattr(args, filter_arg)
			
			if filter_arg == 'include_blobs':
				filter_type = 'include'
			elif filter_arg == 'exclude_blobs':
				filter_type = 'exclude'
			else:
				raise ValueError("unknown filter type: {!r}".format(filter_arg))
			
		else:
			filter_type = default_filter_type
			filter_file = None
			
		return filter_type, filter_file

	def initArgparser(self):
		'''Init argparse ArgumentParser with supported and required parameters.'''
			
		argparser = argparse.ArgumentParser(
			formatter_class=argparse.RawTextHelpFormatter)

		# Parameter group.
		for g in self.supp_groups:
			
			grp = argparser.add_argument_group(g)
			
			# Mutually exclusive group.
			if g in group2mxg:
				mxg = group2mxg[g]
				if self.supp[mxg]:
					gpx = grp.add_mutually_exclusive_group(required=self.req[mxg])
					for p in mutually_exclusive_param_groups[mxg]:
						if self.supp[p]:
							kwargs = deepcopy(param_info[p])
							flags = kwargs.pop('flags')
							kwargs['dest'] = p
							gpx.add_argument(*flags, **kwargs)					
			
			# Add general arguments to main parser, group arguments to argument group.
			handle = argparser if g == 'keyword options' else grp
			
			# Other parameters.
			for p in param_groups[g]:
				if self.supp[p] and p not in param2mxg:
					kwargs = deepcopy(param_info[p])
					flags = kwargs.pop('flags')
					kwargs['dest'] = p
					kwargs['required'] = self.req[p]
					handle.add_argument(*flags, **kwargs)

		# Rename 'optional arguments' to 'keyword options'.
		argparser._optionals.title = "keyword options"

		return argparser		

	def processArgs(self, args):
		'''Process parsed arguments.'''

		a = deepcopy(args)

		if self.supp['input reads']:
			if a.reverse_reads is not None and a.forward_reads is None:
				raise ValueError("input reverse reads specified without corresponding forward reads")
			if a.orphan_reads is not None and a.reverse_reads is None:
				raise ValueError("input orphan reads specified without corresponding paired reads")
			a.main_input_reads = tuple( x for x in 
				(a.forward_reads, a.reverse_reads) 
				if x is not None )
			a.input_reads = tuple( x for x in 
				a.main_input_reads + (a.orphan_reads,)
				if x is not None )
			
		if self.supp['output reads']:
			if a.output_reverse_reads is not None and a.output_forward_reads is None:
				raise ValueError("output reverse reads specified without corresponding forward reads")
			if a.output_orphan_reads is not None and a.output_reverse_reads is None:
				raise ValueError("output orphan reads specified without corresponding paired reads")	  
			a.main_output_reads = tuple( x for x in 
				(a.output_forward_reads, a.output_reverse_reads) 
				if x is not None )
			a.output_reads = tuple( x for x in 
				a.main_output_reads + (a.output_orphan_reads,)
				if x is not None )

		if self.supp['assembly options']:
			if self.supp['assembly']:
				a.assembly_type, a.assembly_file = self.getAssemblyInfo(args)
			if self.supp['len_cutoffs']:
				if len(a.len_cutoffs) > 0:
					a.len_cutoffs = tuple( sorted( set([ non_negative_integer(c) 
						for c in a.len_cutoffs ]) ) )
				else:
					a.len_cutoffs = (100,)
	
		if self.supp['alignment options']: 
			if self.supp['mapping']:	
				a.mapping_type, a.mapping_file = self.getMappingInfo(args)		   
	
		if self.supp['BLAST options']:	
			if self.supp['blast_fields']:
				validateBlastFields(a.blast_fields)	
	
		if self.supp['classifier options']:	 
			
			# Get list of all specified taxa, target and foreign.
			taxa = [ t for k in ('target_taxa', 'foreign_taxa') 
				for t in getattr(a, k) if self.supp[k] ]
			
			if len(taxa) > 0:
				
				taxonomy = NCBITaxa()
				lineages = [ taxonomy.get_lineage(t) for t in taxa ]
				rank_info = taxonomy.get_rank(taxa)
				ranks = [ rank_info[t] for t in taxa ]
			
				for i in range( len(taxa) ):
					
					# Verify that taxon is of a supported rank; make exceptions for special taxa.
					if ranks[i] not in RANKS and taxa[i] not in special_taxa:
						raise ValueError("taxon {!r} is of unsupported rank {!r}".format(taxa[i], ranks[i]))
			
					# Verify that lineages of all taxa are distinct, otherwise this   
					# would introduce conflict or redundancy. Check for duplicates.
					for j in range( i+1, len(taxa) ):
						if taxa[i] in lineages[j]:
							raise ValueError("lineage of taxon {!r} includes taxon {!r}".format(taxa[j], taxa[i]))
						if taxa[j] in lineages[i]:
							raise ValueError("lineage of taxon {!r} includes taxon {!r}".format(taxa[i], taxa[j]))

			# Verify taxa to be ignored are specified as target or foreign taxa.
			if self.supp['classifier_ignore']: 
				for t in a.classifier_ignore:
					if t not in taxa:
						raise RuntimeError("cannot exclude NCBI Taxonomy ID {!r} from "
							"classifier training set; not in training set".format(t))
			
			a.taxa = tuple(sorted(taxa))
			
		if self.supp['filter options']:  
			if self.supp['filter']:
				a.filter_type, a.filter_file = self.getFilterInfo(args)
				
		if self.supp['Phred options']:
			if self.supp['Phred offset']:
				if not a.phred64:
					a.phred33 = True
		
		return a

################################################################################
