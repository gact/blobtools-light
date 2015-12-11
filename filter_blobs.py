#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Filter blob contigs and their corresponding sequencing reads.

Filtering is done in two stages: first generating filter information, then using 
that information to filter read data.

Generating filter information requires a blob contig assembly file, a local copy 
of the NCBI BLAST 'nt' database directory, and a blob table file as output by 
the blobtools-light script 'makeblobs.py'. Target and foreign taxa must also be
specified as NCBI Taxonomy IDs (either as a comma-separated list or as path to a 
file listing one NCBI Taxonomy ID per line). 

NCBI 'nt' database sequences from each taxon are searched with BLAST to assign 
the taxonomy of each contig to one or more of the specified taxa. Each blob 
contig is assigned to a group indicating if it aligns best to target sequences 
('Target'), aligns best to foreign sequences ('Foreign'), or aligns equally well 
to both target and foreign sequences ('Mixed'). Blob contigs without a BLAST hit 
to either group are grouped according to their taxonomic classification in the 
input blob table file (e.g. 'Other', 'no-hit').

Blob contigs confidently assigned as either 'Target' or 'Foreign' are then used 
in the training set of a support vector classifier (SVC), which groups contigs  
into 'Target' and 'Foreign' classes based on their GC-content and coverage depth. 
Blob contigs that were not confidently assigned to 'Target' or 'Foreign' are 
then grouped as 'Target' or 'Foreign' by the classifier. Contigs in the 'Foreign' 
group are marked for removal. Filter information is output as a filter list file 
(listing contigs to include or exclude, one per line); a filter table file; and 
a filter plot file.

The quality of the filter information depends critically on the choice of target 
and foreign taxa, and the results should not be used blindly. Classification of 
blob contigs assumes that the target and foreign taxa can be distinguished by 
their GC-content and coverage. This may not be the case, for example, with 
foreign taxa that have a similar GC-coverage profile to target contigs. In 
such cases, the classifier can be set to ignore such taxa, so that contigs 
from those taxa are taxonomically assigned but not included in the training set.

The filter list file is the key output from filter generation and the key input 
for read filtering. It contains a list of blob contigs, one per line, which 
should be either included or excluded, according to filter type ('include_blobs' 
or 'exclude_blobs'). By default, filter lists indicate the blob contigs to 
exclude. If generating filter information and filtering reads separately, 
filter type must be set consistently.

Reads can be filtered as single-end (SE) or paired-end (PE), with the additional 
option of reading/writing orphan reads. A read mapping file is required to get 
the set of corresponding contigs for each read. A read is removed if most of its 
corresponding contigs are marked for removal. By default, unmapped reads are 
retained, but these can optionally be removed. A paired read with a mate that 
fails the filter will be written to the output orphan read file, if specified.

Note that an email address is required for generating filter information, as it
is required for getting taxon-specific GI lists using the NCBI eutils. 
'''

from __future__ import division
from Bio import SeqIO
from collections import OrderedDict
from copy import deepcopy
from csv import DictWriter
from datetime import datetime
from itertools import izip_longest
from matplotlib import cm
from matplotlib.colors import rgb2hex
import numpy as np
import os
from shutil import copyfile
from sklearn.svm import SVC
import sys
from tempfile import NamedTemporaryFile

from ArgUtil import InputHandler
from ArgUtil import printArgs
from blast_blobs import blastTaxa
from NcbiUtil import getDictOfBestBlastHits
from NcbiUtil import NcbiBlastdb
from NcbiUtil import readBlastTableHeadings
from NcbiUtil import required_blast_fields
import plotblobs
import Util
from Util import default_tax_label_type
from Util import openFastq
from Util import special_plot_colors
from Util import update

# TODO: set maximum number of iterations on classifier

################################################################################

def filterBlobs(args):
	'''Filter contigs by taxonomic affinity and classification.'''

	# Set minimum bitscore difference.
	bitscore_diff = 50.0
	
	# Get arguments.
	assembly_file = args.assembly_file
	blastdb = args.blastdb
	blob_table_file = args.blob_table
	cov_libs = args.cov_libs
	excluding_unannotated = args.exclude_unannotated
	excluding_unmapped = args.exclude_unmapped
	filter_file = args.filter_file
	filter_plot = args.filter_plot	
	filter_table = args.filter_table   
	filter_type = args.filter_type
	foreign_taxa = args.foreign_taxa
	forward_reads = args.forward_reads	
	ignore_taxa = args.classifier_ignore
	main_input_reads = args.main_input_reads
	main_output_reads = args.main_output_reads	
	mapping_file = args.mapping_file
	orphan_reads = args.orphan_reads
	output_forward_reads = args.output_forward_reads
	output_orphan_reads = args.output_orphan_reads	
	output_reverse_reads = args.output_reverse_reads 
	plot_title = args.plot_title
	reverse_reads = args.reverse_reads
	target_taxa = args.target_taxa
	taxa = args.taxa
	tax_label = args.tax_label
	
	# Verify both target and foreign taxa specified.
	if len(target_taxa) == 0:
		raise ValueError("no target taxa specified")
	if len(foreign_taxa) == 0:
		raise ValueError("no foreign taxa specified")
		
	# Assume no action required.
	writing_filter_list = writing_filter_files = filtering_reads = False

	# Set actions required to output filtered reads, if specified.
	if len(main_output_reads) > 0:
		filtering_reads = True
		if filter_file is None or not os.path.exists(filter_file):
			writing_filter_list = True

	# Set action required to write filter list, if specified.
	if filter_file is not None and not os.path.exists(filter_file):
		writing_filter_list = True		
			 
	# Set action required to output filter files, if specified.
	if any([writing_filter_list, filter_table is not None, filter_plot is not None]):
		writing_filter_files = True
	 
	# Verify at least one output specified.
	if not any([writing_filter_files, filtering_reads]):
		raise RuntimeError("no output specified")

	# If filtering reads, verify necessary input files exist.
	if filtering_reads:
		read_filter_input = {
			'input reads': len(main_input_reads) > 0, 
			'read mapping file': mapping_file is not None and os.path.exists(mapping_file)
		} 
		for k, test in read_filter_input.items():
			if test == False:
				raise RuntimeError("cannot filter reads without {}".format(k))			

	# If writing filter files, verify necessary input files exist.
	if writing_filter_files:
		contig_filter_input = {
			'assembly file': assembly_file is not None and os.path.exists(assembly_file), 
			'blob table file': blob_table_file is not None and os.path.exists(blob_table_file), 
			'BLAST database': blastdb is not None
		}
		for k, test in contig_filter_input.items():
			if test == False:
				raise RuntimeError("cannot generate filter without {}".format(k))

	# Verify filter options selected in correct context.
	if excluding_unannotated and not writing_filter_files:
		raise RuntimeError("cannot exclude unannotated contigs - not writing filter files")		
	if excluding_unmapped and not filtering_reads:
		raise RuntimeError("cannot exclude unmapped reads - not filtering reads")

	update("Starting blobtools-light filter pipeline")

	update("Creating temp directory")
	with Util.tempDirectory() as twd:
		update("Created temp directory: {!r}".format(twd))

		# Set filter list filepath to that specified or to a temp file.
		if filter_file is None:
			filter_temp = NamedTemporaryFile(mode='w', dir=twd, delete=False)
			filter_list = filter_temp.name
		else:
			filter_list = filter_file

		if writing_filter_files:
			
			update("Writing filter files")
			
			filters = {'tax_filter': dict(), 'class_filter': dict()}
			tax_group_names = ('target', 'foreign')
			target, foreign = (0, 1)
					
			update("Checking NCBI BLAST 'nt' database in directory: {!r}".format(blastdb))
			if not NcbiBlastdb.present(blastdb, 'nt'):
				NcbiBlastdb('nt').download(blastdb)

			update("Reading blob table file {!r}".format(blob_table_file))
			blob_data = Util.readBlobTable(blob_table_file, cov_libs=cov_libs)
			
			update("Fetching contig IDs from {!r}".format(assembly_file))
			contig_ids = Util.getFastaSeqids(assembly_file)
						
			for contig_id in blob_data.keys():
				if contig_id not in contig_ids:
					raise ValueError("blob contig {!r} not found in assembly".format(contig_id))
						
			update("Starting check of taxonomic affinity for target and foreign taxa")
			anno = dict()
			for g, tax_group in enumerate([target_taxa, foreign_taxa]):
				
				for taxid in tax_group:
					
					update("Checking taxonomic affinity for {} taxon {!r}".format(
						tax_group_names[g], taxid))
					
					# Set BLAST arguments.
					a = deepcopy(args)
					a.blast_fields = required_blast_fields
					a.blast_file = os.path.join(twd, '{}.blast'.format(taxid))
					a.blast_task = 'blastn'
					a.target_taxa = [taxid]
				
					# Run BLAST on taxon.
					blastTaxa(a)
					
					# Get best BLAST hits for this taxon.
					hits = getDictOfBestBlastHits(a.blast_file)
					
					# Add taxonomic affinity annotation for this taxon.
					for c in hits:
					
						# Ensure annotation entry exists for contig,
						# both for target and foreign annotations.
						anno.setdefault(c, [None, None])
					
						# Get best E-value and bitscore for contig in this taxon.
						evalue, bitscore = [ hits[c][k] 
							for k in ('evalue', 'bitscore') ]
						
						# If previous BLAST hit for this contig, compare E-values 
						# and bitscores. Skip new hit if its scores are not at 
						# least equal to the best hit for this contig in this 
						# taxon group. If scores match best hit, add Taxonomy ID.
						if anno[c][g] is not None:
							e, s = [ anno[c][g][k] for k in ('evalue', 'bitscore') ]
							if evalue > e:
								continue					
							elif evalue == e:
								if bitscore < s:
									continue
								elif bitscore == s:
									anno[c][g]['taxids'].add(taxid)
									continue
					
						# Set taxonomic affinity annotation for this contig.
						anno[c][g] = { 'taxids': set([taxid]), 'evalue': evalue,
							'bitscore': bitscore }
			
			update("Completed check of taxonomic affinity for target and foreign taxa")
					
			update("Grouping contigs by taxonomic affinity")
			training_contigs = list()
			for c in contig_ids:
				
				# If contig has affinity to any target or foreign taxa, set
				# as target, foreign, or mixed annotation as appropriate..
				if c in anno:
										
					if all( a is not None for a in anno[c] ):
						
						e = [ a['evalue'] for a in anno[c] ]
						s = [ a['bitscore'] for a in anno[c] ]
						
						if e[target] < e[foreign] or (e[target] == e[foreign] 
							and s[target] - s[foreign] >= bitscore_diff):
							status = 'Target'
						elif e[foreign]< e[target] or (e[foreign] == e[target] 
							and s[foreign] - s[target] >= bitscore_diff):
							status = 'Foreign'
						else:
							status = 'Mixed'
							
					elif anno[c][target] is not None:
					
						status = 'Target'
						
					elif anno[c][foreign] is not None:
					
						status = 'Foreign'
				
				# ..otherwise distinguish contigs annotated as another taxon, 
				# contigs with ambiguous annotation, and contigs with none.
				elif blob_data[c]['tax'] in ('no-hit', 'ambig-hit'):
					status = blob_data[c]['tax']
				else:
					status = 'Other'
				
				# If contig has clear annotation and classifier is 
				# not ignoring annotated taxon, add to training set.
				if ( (status == 'Target' and not all(t in ignore_taxa for t in anno[c][target]['taxids'])) or 
					(status == 'Foreign' and not all(t in ignore_taxa for t in anno[c][foreign]['taxids'])) ):
					training_contigs.append(c)
				
				# Set annotation filter status for contig.
				filters['tax_filter'][c] = status
	
			# Set training data.
			training_data = np.zeros([len(training_contigs), 2], dtype=np.float)
			training_taxa = np.zeros([len(training_contigs), ], dtype=np.int)
			for i, c in enumerate(training_contigs):
				training_data[i] = [ blob_data[c]['gc'], blob_data[c]['cov'] ]
				training_taxa[i] = 1 if filters['tax_filter'][c] == 'Target' else -1
		
			# Verify that training data has both target and foreign taxa.
			if 1 not in training_taxa:
				raise RuntimeError("classifier training set has no data from target taxa")
			if -1 not in training_taxa:
				raise RuntimeError("classifier training set has no data from foreign taxa")
		
			# Set test data.
			contig_gc = np.array([blob_data[c]['gc'] for c in contig_ids ], dtype=np.float)
			contig_cov = np.array([blob_data[c]['cov'] for c in contig_ids ], dtype=np.float)
			test_data = np.array( zip(contig_gc, contig_cov), dtype=np.float)
		
			# Init classifier; set class_weight to 'auto' (weight by class size).
			classifier = SVC(class_weight='auto')
			
			update("Training classifier")
			classifier.fit(training_data, training_taxa)

			update("Predicting taxonomic affinity of contigs based on coverage and GC-content")
			predictions = classifier.predict(test_data)

			# Set classifier filter value.
			for c, prediction in zip(contig_ids, predictions):
				if prediction == 1:
					filters['class_filter'][c] = 'include'
				elif prediction == -1:
					filters['class_filter'][c] = 'exclude'
				else:
					raise RuntimeError("unknown classifier label: {!r}".format(prediction))				
			   
			update("Combining filter classifications")
			filt = dict()
			for c in contig_ids:
				if filters['tax_filter'][c] == 'Target':
					filt[c] = 'include'
				elif filters['tax_filter'][c] == 'Foreign':
					filt[c] = 'exclude'
				elif excluding_unannotated and filters['tax_filter'][c] == 'no-hit':
					filt[c] = 'exclude'
				else:   
					filt[c] = filters['class_filter'][c]
			
			if writing_filter_list:
				if filter_file is not None:
					update("Writing list of contigs to {}: {!r}".format(
						filter_type, filter_file))
				with open(filter_list, 'w') as fh:
					for c in contig_ids:
						if filt[c] == filter_type:
							fh.write("{}\n".format(c))
									
			if filter_table is not None:
				
				update("Writing filter table to {!r}".format(filter_table))
				
				filter_headings = ('contig_id', 'gc', 'cov', 'tax_filter', 
					'class_filter', 'filter')
								
				with NamedTemporaryFile(mode='w', dir=twd, suffix='.tab', delete=False) as ftab:
				
					filter_table_file = ftab.name
				
					ftab.write("# {}\n# {}\n# {}\n# {}\n".format(sys.argv[0],
						datetime.now().strftime("%m/%d/%Y %H:%M:%S"), 
						' '.join(sys.argv), '\t'.join(filter_headings)))
				
					writer = DictWriter(ftab, fieldnames=filter_headings, 
						dialect='text-tab')
								
					for c in contig_ids:
						row = { k: blob_data[c][k] for k in ('contig_id', 'gc', 'cov') }
						row.update({ k: filters[k][c] for k in ('tax_filter', 'class_filter') })
						row.update({ 'filter': filt[c] })
						writer.writerow(row)
						
				copyfile(filter_table_file, filter_table)
							   
			if filter_plot is not None:
				
				update("Writing filter plot to {!r}".format(filter_plot))
				
				tax_groups = getFilterPlotTaxGroups(target_taxa, foreign_taxa)
				color_dict = getFilterPlotColorDict()
				label_dict = getFilterPlotLabelDict(target_taxa, foreign_taxa, 
					tax_label=tax_label)
								
				# Set plot data.
				plot_data = np.array([ [ d['contig_id'], d['length'], d['gc'], 
					filters['tax_filter'][b] ] for b, d in blob_data.items() ])
				plot_cov = np.array([ d['cov'] for d in blob_data.values() ])
				
				plotblobs.plot(plot_data, plot_cov, filter_plot, plot_title, 
					tax_groups=tax_groups, color_dict=color_dict, 
					label_dict=label_dict, classifier=classifier)
			
			update("Filter files written")
				
		if filtering_reads:
			
			update("Starting read filtering")
			
			# Init dictionary of FASTQ read IDs (for finding duplicates).
			seqids = dict()
			
			reads_per_fragment = len(main_input_reads)
			read_type = 'PE' if reads_per_fragment == 2 else 'SE'
			
			with open(filter_list, 'r') as fh:
				update("Reading list of contigs to {}: {!r}".format(filter_type, 
					filter_list))
				filter_set = set([ line.rstrip() for line in fh ])
			
			reads, read2ctgs, contigs = Util.readMappingInfo(mapping_file)
			
			if filter_type == 'exclude':
				contig_passes = { c: c not in filter_set for c in contigs }
			elif filter_type == 'include':
				contig_passes = { c: c in filter_set for c in contigs }
			
			temp_output_reads = [ os.path.join( twd, os.path.basename(x) ) 
				for x in main_output_reads ]
			
			try:
				fin = [ openFastq(i, 'r') for i in main_input_reads ]
				fout = [ openFastq(o, 'w') for o in temp_output_reads ]
				
				if output_orphan_reads is not None:
			
					temp_orphan_reads = os.path.join(twd, 
						os.path.basename(output_orphan_reads))
				
					foo = openFastq(temp_orphan_reads, 'w')
				
					if orphan_reads is not None:
						
						update("Filtering existing orphan reads")
						
						with openFastq(orphan_reads, 'r') as fio:
						
							for record in SeqIO.parse(fio, 'fastq'):
						
								r = record.id
								
								# Increment read ID count.
								seqids.setdefault(r, 0)
								seqids[r] += 1
								
								# Get mapping of read to corresponding contigs.
								try:
									read_contig_map = read2ctgs[r]
								except KeyError as e:
									if not excluding_unmapped:
										read_contig_map = None
									else:
										continue
								
								# If read maps to a contig, filter by contig status.
								if read_contig_map is not None:
									
									# Verify correct read type in read-contig mapping.
									if len(read_contig_map) != 1:
										raise ValueError("read type mismatch for {!r}".format(r))
									
									# If most corresponding contigs pass, read passes.
									ctgs = read_contig_map[0].keys()
									if len(ctgs) == 1:
										if not contig_passes[ ctgs[0] ]:
											continue
									else:
										passes = sum( 1 if contig_passes[c] else 0 for c in ctgs )
										if passes / len(ctgs) < 0.5:
											continue
								
								SeqIO.write(record, foo, 'fastq')
				
				
				update("Filtering {} reads".format(read_type))
				
				fastq_iters = [ SeqIO.parse(fh, 'fastq') for fh in fin ]
				
				for records in izip_longest( *fastq_iters ):

					if any([ record is None for record in records ]):
						raise RuntimeError("record count mismatch in paired input read files")
	
					rids = [ record.id for record in records ]
	
					if not all( r == rids[0] for r in rids[1:] ):
						raise RuntimeError("record ID mismatch in paired input read files")
		
					r = rids[0]
			
					# Increment read ID count.
					seqids.setdefault(r, 0)
					seqids[r] += 1
					
					# Get mapping(s) of read(s) to corresponding contigs.
					try:
						read_contig_map = read2ctgs[r]
					except KeyError as e:
						if not excluding_unmapped:
							read_contig_map = None
						else:
							continue
					
					# If read maps to a contig, filter by contig status.
					if read_contig_map is not None:
						
						# Verify correct read type in read-contig mapping.
						if len(read_contig_map) != reads_per_fragment:
							raise ValueError("read type mismatch for {!r}".format(r))
						
						# For each read in fragment, if most of the
						# corresponding contigs pass, read passes.
						results = [False] * reads_per_fragment
						for i in range(reads_per_fragment):
							try:
								ctgs = read_contig_map[i].keys()
							except AttributeError:
								results[i] = False
								continue
							if len(ctgs) == 1:
								results[i] = contig_passes[ ctgs[0] ]
							else:
								passes = sum( 1 if contig_passes[c] else 0 for c in ctgs )
								results[i] = passes / len(ctgs) >= 0.5
					
					# ..otherwise pass all reads.
					else:
						results = [True] * reads_per_fragment
					
					# If SE read passes or PE read pair pass, write to output..
					if all( result == True for result in results ):
						for record, fh in zip(records, fout): 
							SeqIO.write(record, fh, 'fastq')
					# ..otherwise if none pass, continue to next read(s)..
					elif all( result == False for result in results ):
						continue
					# ..otherwise if writing orphan reads,  
					# output each orphan read that passes.
					elif output_orphan_reads:
						for result, record in zip(results, records):
							if result == True:
								SeqIO.write(record, foo, 'fastq')
			finally:
				if output_orphan_reads is not None: foo.close()
				for fh in fin + fout: fh.close()
			
			# Verify no duplicate read IDs.
			if any( seqids[r] > 1 for r in seqids ):
				raise RuntimeError("duplicate FASTQ headers in input read files")  
			
			update("Writing main output reads to {!r}".format(', '.join(main_output_reads)))
			for temp_file, output_file in zip(temp_output_reads, main_output_reads):
				copyfile(temp_file, output_file)
			
			if output_orphan_reads is not None:
				update("Writing filtered orphan reads to {!r}".format(output_orphan_reads))
				copyfile(temp_orphan_reads, output_orphan_reads)
	
			update("Completed read filtering")
			
	update("Completed blobtools-light filter pipeline")

def getFilterPlotColorDict():
	'''Get colour dict for filter plot.'''

	# Set colours for filter categories.
	colormap = cm.get_cmap('coolwarm_r')
	target_rgba = colormap(1.0)
	foreign_rgba = colormap(0.0)
	both_rgba = np.array( (target_rgba, foreign_rgba) )
	mix_rgba = np.mean(both_rgba, axis=0).tolist()
	
	# Set filter classification colour dict.
	color_dict = {
		'Target':    rgb2hex(target_rgba), 
		'Foreign':   rgb2hex(foreign_rgba),
		'Mixed':     rgb2hex(mix_rgba),
		'other-taxa':special_plot_colors['other-taxa'],
		'ambig-hit': special_plot_colors['ambig-hit'],
		'no-hit':    special_plot_colors['no-hit']		
	}	

	return color_dict

def getFilterPlotLabelDict(target_taxa, foreign_taxa, 
	tax_label=default_tax_label_type):
	'''Get label dict for filter plot.'''
	
	# Get list of all taxa.
	taxa = [ t for t in target_taxa + foreign_taxa ]
	
	# Map each taxon to a taxon label.
	tax2label = { t: Util.formatTaxon(t, tax_label=tax_label) for t in taxa }
				
	# Set filter classification label dict.
	label_dict = {
		'Target': "Target taxa [{}]".format(', '.join( tax2label[t] for t in target_taxa )), 
		'Foreign': "Foreign taxa [{}]".format(', '.join( tax2label[t] for t in foreign_taxa )),
		'Mixed': "Mixed annotation",
		'other-taxa': "Other annotation",
		'ambig-hit': "Ambiguous annotation ['ambig-hit']",
		'no-hit': "Unannotated ['no-hit']"
	}
	
	return label_dict

def getFilterPlotTaxGroups(target_taxa, foreign_taxa):
	"""Return OrderedDict of tax groups for filter plot."""
	return OrderedDict([ (g, [g]) for g in ('Target', 'Foreign', 'Mixed', 
		'other-taxa', 'ambig-hit', 'no-hit') ])

################################################################################

supported_params = ('assembly', 'blastdb', 'blob_table', 'classifier options',  
	'cov_libs', 'email', 'filter options', 'filter_plot', 'filter_table', 
	'input reads', 'mapping', 'num_threads', 'output reads', 'plot_title', 
	'tax_label')
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
		
	filterBlobs(args)

	update("Done.")

################################################################################
	
if __name__ == '__main__':
	main(sys.argv)

################################################################################
