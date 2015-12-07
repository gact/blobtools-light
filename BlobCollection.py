#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   		: BlobCollection.py
Version 	: 0.1
Author 		: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs 		: ?
To do 		: Add to CAS parser that it checks that the assembly is the correct one
"""

from __future__ import division
import commands
from ete2 import NCBITaxa
import numpy
import re
import subprocess
import sys
import time

from MiscFunctions import keyWithMaxVal
from MiscFunctions import n50
from NcbiUtil import getDictOfBestBlastHits
from Util import formatTaxon
from Util import getConsensusTaxonomyID
from Util import getCorrectedNucSeqLength
from Util import imprecise_sciname_prefixes
from Util import pseudo_domains
from Util import RANKS
from Util import radix_vitae
from Util import special_taxa

class Blob():
	'''
	Blobs contain the information of a contig in an assembly.
	Each line that eventually gets printed into the blobplot file is a blob. 
	'''

	def __init__(self, header, seq):
		self.name = header
		self.seq = seq.upper()
		self.length = len(seq)
		self.corrected_length = getCorrectedNucSeqLength(self.seq)
		self.gc = float((self.seq.count('G') + self.seq.count('C') ) / self.corrected_length ) if self.corrected_length > 0 else 0.0
		self.covs = dict()
		self.tax = dict() 

class BlobCollection():
	'''
	The BlobCollection class contains the blobs and all methods needed 
	for parsing the input files and doing the necessary computations 
	for arriving at the final results. 
	'''
	def __init__(self, rank_level=None, tax_label='scientific name'):
		self.contigs = dict()
		self.index = list()
		self.outfile = str()
		self.stats = dict()
		self.cov_libs = list()
		self.blast_libs = list()
		self.rank_level = rank_level
		self.tax_label = tax_label
		self.blast_order = list()

	def addBlob(self, blob):
		'''
		Adds a blob to the BlobCollection
			- puts it in the contigs dict with header as the key
			- puts it header in the index list so that it can be found by index (this is for CAS files and printing blobs in order)
		'''
		if not blob.name in self.contigs: 
			self.contigs[blob.name] = blob
		else: 
			sys.exit("[ERROR] - Sequence header {} occurs more than once".format(blob.name))
		self.index.append(blob.name)

	def getBlobsFromAssembly(self, assembly_file, assembly_type, exclude_assembly_cov):
		print "[STATUS] - Parsing assembly %s" % (assembly_file)
		header, seq = '', ''
		with open(assembly_file) as fh:
			for line in fh:
				line_data = line.rstrip("\n")
				if line.startswith('>'):
					if (seq): 
						blob = Blob(header, seq) 
						self.addBlob(blob)
						seq = ''
						if assembly_type == 'unknown' or exclude_assembly_cov:
							pass
						else:
							cov = float(self.parseCovFromHeader(header, assembly_type))
							self.addBlobCov(header, assembly_type, cov)
					header = line.rstrip("\n").lstrip(">")
				else:
					seq += line.rstrip("\n").upper() 
			blob = Blob(header, seq) 
			self.addBlob(blob)
			if assembly_type == 'unknown' or exclude_assembly_cov:
				pass
			else:
				cov = float(self.parseCovFromHeader(header, assembly_type))
				self.addBlobCov(header, assembly_type, cov)
		if not assembly_type == 'unknown' and not exclude_assembly_cov:
			self.cov_libs.append(assembly_type)

	def addBlobCov(self, header, mapping_lib, cov):
		'''
		Adds coverage to the covs dict: key = cov_lib name, value = coverage
		'''
		#if mapping_lib in self.contigs[header].covs:
		#	sys.exit("[ERROR] - Contig {} received more than one coverage from the {} assembly file".format( header, mapping_lib))
		#else:
		self.contigs[header].covs[mapping_lib] = float(cov)

	def parseCovFromHeader(self, header, assembly_type):
		''' 
		Returns the coverage from the header of a FASTA 
		sequence depending on the assembly type
		'''
		if assembly_type == 'spades':
			return header.split("_")[-3]
		elif assembly_type == 'velvet':
			return header.split("_")[-1]
		elif assembly_type == 'abyss':
			temp = header.split(" ")
			return temp[2]/(temp[1]+1-75)
		else:
			sys.exit("[ERROR] - Coverage could not be parsed from header {} of type {}".format( header, assembly_type ))

	def getCovForBlobs(self, mapping_files):
		'''
		- Coordinates parsing of coverages from different mapping files
		- If more than one cov_lib is specified then a new cov_lib is created for each blob which contain the sum of the other cov_libs 
		'''
		for cov_lib, mapping_file in mapping_files.items():
			self.cov_libs.append(cov_lib)
			print "[STATUS] - Parsing coverage from cov_lib \"{}\" from file {}".format(cov_lib, mapping_file)
			if cov_lib.startswith("CAS"):
				self.parseCovFromCasFile(cov_lib, mapping_file)
			elif cov_lib.startswith("BAM"):
				self.parseCovFromBAMFile(cov_lib, mapping_file)
			elif cov_lib.startswith("SAM"):
				self.parseCovFromSAMFile(cov_lib, mapping_file)
			elif cov_lib.startswith("COV"):
				self.parseCovFromCovFile(cov_lib, mapping_file)
			else:
				sys.exit("[ERROR] - Unknown cov_lib type {} in {}".format(cov_lib, mapping_file))	
	 
		for contig_name in self.contigs:
			cov_sum = 0.0
			for cov_lib in sorted(self.cov_libs):
				if cov_lib in self.contigs[contig_name].covs:	
					cov_sum += self.contigs[contig_name].covs[cov_lib]	
				else:
					self.contigs[contig_name].covs[cov_lib] = 0.0		
			if len(self.cov_libs) > 1:
				self.contigs[contig_name].covs['SUM'] = cov_sum
		if len(self.cov_libs) > 1:
			self.cov_libs.append("SUM")

	def parseCovFromCasFile(self, lib_name, cas_file):
		'''
		Parse coverage from CAS file
		'''
		error, message = commands.getstatusoutput("clc_mapping_info -s -f " + cas_file)
		if (error):
			sys.exit("[ERROR] - Please add clc_mapping_info to you PATH variable.") 
		p = subprocess.Popen("clc_mapping_info -n " + cas_file , stdout=subprocess.PIPE, bufsize=1, shell=True)
		read_counter = 1
		cas_line_re = re.compile(r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d{2})\s+(\d+)\s+(\d+.\d{2})")
		for line in iter(p.stdout.readline, b''):
			match = cas_line_re.search(line)
			if match:
				contig_index = int(match.group(1)) - 1 # -1 because index of contig list starts with zero 
				contig_cov = float(match.group(4))
				try:
					contig_id = self.index[contig_index]
					self.addBlobCov(contig_id, lib_name, contig_cov)
				except IndexError, e:
					# throw IndexError if there is a different amount of contigs in CAS than there is in the assembly
					sys.exit("[ERROR] : There are more contigs in the CAS file than in the assembly. It seems that the mapping has not been performed against this assembly.")	

	def parseCovFromBAMFile(self, lib_name, bam_file):
		'''
		Parse coverage from BAM file
		'''
		contig_base_cov = {}
		bam_line_re = re.compile(r"\S+\s+\d+\s+(\S+)\s+\d+\s+\d+\s+(\S+)")
		cigar_match_re = re.compile(r"(\d+M)") # only counts M's
		error, total_read_count = commands.getstatusoutput("samtools view -c -F 4 " + bam_file)
		if not (total_read_count):
			sys.exit("[ERROR] - Please add samtools to you PATH variable.") 
		p = subprocess.Popen("samtools view -F 4 " + bam_file , stdout=subprocess.PIPE, bufsize=1, shell=True)
		read_counter = 0
		for line in iter(p.stdout.readline, b''):
			match = bam_line_re.search(line)
			if match:
				read_counter += 1
				contig_name = match.group(1)
				contig_cigar_string = match.group(2)
				matchings = cigar_match_re.findall(contig_cigar_string)
				sum_of_matchin_bases = 0	
				for matching in matchings:
					sum_of_matchin_bases += int(matching.rstrip("M"))
				contig_base_cov[contig_name] = contig_base_cov.get(contig_name, 0) + sum_of_matchin_bases
				if read_counter % 5000 == 0:		
					sys.stdout.write('\r')
					progress = int(read_counter)/int(total_read_count)
					print "\t[PROGRESS] - " + format(float(progress),'.2%'),
					sys.stdout.flush()
		sys.stdout.write('\r')
		print "\t[PROGRESS] - 100.00%"
		for contig_id, base_cov in contig_base_cov.items():
			if contig_id not in self.contigs:
				sys.exit("[ERROR] - Sequence header %s in %s does not seem to be part of the assembly. FASTA header of sequence in assembly MUST equal reference sequence name in BAM file. Please check your input files." % (contig_id, bam_file)) 
			else:
				contig_cov = base_cov / self.contigs[contig_id].corrected_length
				self.addBlobCov(contig_id, lib_name, contig_cov)

	def parseCovFromSAMFile(self, lib_name, sam_file):
		'''
		Parse coverage from SAM file
		'''
		contig_base_cov = {}
		sam_line_re = re.compile(r"\S+\s+(\d)+\s+(\S+)\s+\d+\s+\d+\s+(\S+)")
		cigar_match_re = re.compile(r"(\d+M)") # only counts M's
		with open(sam_file) as fh:
			for line in fh:
				match = sam_line_re.search(line)
				if match:
					sam_flag = match.group(1)
					contig_name = match.group(2)
					if contig_name == '*' or sam_flag == 4:
						pass
					else:
						contig_cigar_string = match.group(3)
						matchings = cigar_match_re.findall(contig_cigar_string)
						sum_of_matchin_bases = 0	
						for matching in matchings:
							sum_of_matchin_bases += int(matching.rstrip("M"))
						contig_base_cov[contig_name] = contig_base_cov.get(contig_name, 0) + sum_of_matchin_bases
		for contig_id, base_cov in contig_base_cov.items():
			if contig_id not in self.contigs:
				sys.exit("[ERROR] - Sequence header %s in %s does not seem to be part of the assembly. FASTA header of sequence in assembly MUST equal reference sequence name in BAM file. Please check your input files." % (contig_id, sam_file)) 
			else:
				contig_cov = base_cov / self.contigs[contig_id].corrected_length
				self.addBlobCov(contig_id, lib_name, contig_cov)

	def parseCovFromCovFile(self, lib_name, cov_file):
		'''
		Parse coverage from COV file
		'''
		cov_line_re = re.compile(r"^(\S+)\t(\d+\.\d+)")
		with open(cov_file) as fh:
			for line in fh:
				match = cov_line_re.search(line)
				if match:
					contig_id, contig_cov = match.group(1), float(match.group(2))
					self.addBlobCov(contig_id, lib_name, contig_cov)

	def getTaxForBlobs(self, blast_files, target_taxa=()):
		'''Calculate consensus taxonomy for blob contigs within each BLAST library.
		
		For each contig within each BLAST file, take the Taxonomy IDs for the
		set of hits with the best BLAST E-value and score, get a consensus 
		Taxonomy ID of the best BLAST hit sequences (the consensus taxid, 
		E-value, and bitscore are stored in blob.tax[blast_lib]). 
		
		The consensus taxid is set to None for blob contigs without a BLAST hit,
		and later converted to 'no-hit'. If the consensus taxid is '1' (root of 
		all life), this indicates that an unambiguous consensus taxid could not
		be found, so this is later converted to 'ambig-hit'.
		'''
		
		print("[STATUS] - Getting taxonomy of blob contigs")
		
		# Cache reciprocals to avoid needless division.
		reciprocal = dict()
		
		# Load NCBI Taxonomy data.
		taxonomy = NCBITaxa()
		
		# Get all taxids in the clades defined by target taxa.
		if len(target_taxa) > 0:
			target_taxids = set()
			for target_taxon in target_taxa:
				clade_leaves = taxonomy.get_descendant_taxa(target_taxon)
				for clade_leaf in clade_leaves:
					l = taxonomy.get_lineage(clade_leaf)
					for t in l[l.index(target_taxon):]:
						target_taxids.add(t)

		# Process each BLAST result file.
		for blast_lib, blast_file in blast_files.items():
		
			self.blast_libs.append(blast_lib)
						
			hits = getDictOfBestBlastHits(blast_file)
			
			# If target taxa specified, give preference to hits within target taxa.
			if len(target_taxa) > 0:
				for q in hits:
					if any( t in target_taxids for t in hits[q]['taxa'] ):
						rejected_taxids = [ t for t in hits[q]['taxa'] 
							if t not in target_taxids ] 
						for t in rejected_taxids:
							del hits[q]['taxa'][t]

			# Get Taxonomy ID for each contig.
			for contig_name in self.contigs:
			
				taxid = evalue = bitscore = None
			
				if contig_name in hits:
				
					# Get taxids, E-value, bitscore.
					taxids = hits[contig_name]['taxa'].keys()
					evalue = hits[contig_name]['evalue']
					bitscore = hits[contig_name]['bitscore']
													
					# If there are multiple Taxonomy IDs, get a consensus taxid..
					if len(taxids) > 1:
						
						# Map BLAST result line indices to taxids.
						index2taxids = dict()
						for t in taxids:
							for i in hits[contig_name]['taxa'][t]:
								index2taxids.setdefault(i, set())
								index2taxids[i].add(t)
						
						# Weight taxa by number of BLAST result lines. Each line 
						# adds weight 1 to taxon. If many taxa on a single BLAST 
						# result line, the weight added to each taxon is the 
						# reciprocal of the number of taxids in that line.
						taxid_weights = { t: 0 for t in taxids }
						for i in index2taxids:
							for t in index2taxids[i]:
								n = len(index2taxids[i])
								try:
									taxid_weights[t] += reciprocal[n] 
								except KeyError:
									reciprocal[n] = 1 / n
									taxid_weights[t] += reciprocal[n]
						
						# Get the consensus Taxonomy ID.
						taxid = getConsensusTaxonomyID(taxids, 
							taxid_weights=taxid_weights)
						
						if taxid == radix_vitae: 
							
							# Get minimal topology of the specified taxids;
							# its root will be the root of life.
							tree = taxonomy.get_topology(taxids)
							
							# Get children of the root of life.
							kids = [ kid for kid in tree.get_children() ]
							
							# Filter children: keep only real domains of life.
							domains = [ node for node in kids 
								if node.taxid not in pseudo_domains ]
							
							# If one domain remains, get the consensus 
							# taxid for the taxids within that domain.
							if len(domains) == 1:
								ancestor = domains.pop()
								taxids = [ node.taxid 
									for node in ancestor.get_descendants() ]
								for t in taxid_weights:
									if t not in taxids:
										del taxid_weights[t]
								taxid = getConsensusTaxonomyID(taxids, 
									taxid_weights=taxid_weights)
								
					# ..otherwise take the single Taxonomy ID.														
					else:
						taxid = taxids.pop()
			
				# If unambiguous consensus taxid found, adjust as appropriate.
				if taxid not in (None, radix_vitae):
					
					# Get lineage of taxid: from taxid to root of life, inclusive.
					lineage = [ t for t in reversed(taxonomy.get_lineage(taxid)) ]
				
					# Check if taxid is in a 'special' taxon (e.g. vectors).
					special_taxids = [ t for t in lineage if t in special_taxa ]
				
					# If special taxids found, set taxid to closest in lineage..
					if len(special_taxids) > 0:
					
						taxid = special_taxids[0]
					
					# ..otherwise adjust taxid as normal.
					else:
						tax2sciname = taxonomy.get_taxid_translator(lineage)
						tax2rank = taxonomy.get_rank(lineage)
						
						# Check for 'imprecise' scientific names (e.g. 'environmental samples'). 
						# Set taxid to that of first ancestral taxon in a supported rank.
						for i, t in enumerate(lineage):
							if any( tax2sciname[t].startswith(p) for p in imprecise_sciname_prefixes ):
								for j, u in enumerate(lineage[i+1:]):
									if tax2rank[u] in RANKS:
										lineage = lineage[j:]
										taxid = u
										break
								break
								
						# Set taxid to first lineage taxon in a supported rank.
						for i, t in enumerate(lineage):
							if tax2rank[t] in RANKS:
								lineage = lineage[i:]
								taxid = t
								break
						
						# If preferred rank specified, seek taxid for that rank.
						if self.rank_level is not None:
							for t in lineage:
								if tax2rank[t] == self.rank_level:
									taxid = t
									break	

				# Set consensus taxid for this contig in this BLAST library.
				self.contigs[contig_name].tax[blast_lib] = { 'taxid': taxid,
					'evalue': evalue, 'bitscore': bitscore }

	def getConsensusTaxForBlobs(self, taxrule, blast_order):
		'''Get consensus taxonomy for blob contigs across BLAST libraries.
		
		Based on taxrule ("A" or "B") and the blast_order (list in order in  
		which blast files where specified) it calculates the consensus taxonomy  
		for each blob contig. 
		- if taxrule == A:
			- if a blob contig has hits in more than one BLAST file, the BLAST 
			  hit with the best E-value/bitscore is taken as the best hit, and 
			  the taxid of this BLAST hit is taken as the consensus taxid
		- if taxrule == B:
		    - if a blob contig has hits in more than one BLAST file, the BLAST 
		      hit is taken from the first BLAST result file in which it is found, 
		      and the taxid of this hit is taken as the consensus taxid
		'''
		
		print("[STATUS] - Getting consensus taxonomy of blob contigs")
		
		for contig_name in self.contigs:
		
			consensus = { 'taxid': None, 'evalue': None, 'bitscore': None }
			
			for blast_lib in blast_order:
			
				# Get best BLAST hit for this contig in this BLAST library.
				taxid, evalue, bitscore = [ 
					self.contigs[contig_name].tax[blast_lib][k] 
					for k in ('taxid', 'evalue', 'bitscore') ]
				
				if taxrule == 'A':
					
					# If BLAST hit found for this contig in this BLAST library 
					# and E-value and bitscore are better than previous best, 
					# set consensus taxid, E-value, bitscore from this hit.
					if taxid is not None and consensus['taxid'] is not None:
						e, s = [ consensus[k] for k in ('evalue', 'bitscore') ]
						if evalue >= e or (evalue == s and bitscore <= s):
							continue
					consensus = { 'taxid': taxid, 'evalue': evalue, 
						'bitscore': bitscore }
						
				elif taxrule == 'B':
					
					# If BLAST hit found for this contig in this BLAST library, 
					# set consensus taxid, E-value, bitscore; then continue to 
					# next contig without looking in subsequent BLAST libraries.
					if taxid is not None:
						consensus = { 'taxid': taxid, 'evalue': evalue, 
							'bitscore': bitscore }
						break
			
			# Set consensus taxid for this contig.
			self.contigs[contig_name].tax['tax'] = consensus
		
		# Add consensus BLAST library 'tax' to BLAST libraries.
		self.blast_libs.append('tax')
	
	def printCOVToFiles(self, mapping_files):
		for cov_lib, cov_file in mapping_files.items():
			if cov_lib.startswith('BAM') or cov_lib.startswith('SAM') or cov_lib.startswith('CAS'):
				cov_fh = open(cov_file + ".cov", 'w') 
				for contig_name in self.index:
					blob = self.contigs[contig_name]
					cov_fh.write(blob.name + "\t" + str(blob.covs[cov_lib]) + "\n")
				cov_fh.close() 
			else:
				pass

	def getStats(self):
		'''
		Calculates all different kind of stats for each taxonomic group based on each BLAST file and the consensus... 
		'''
		self.stats['count'] = {}
		self.stats['span']= {}
		self.stats['n50']= {}
		self.stats['lengths']= {}
		self.stats['gc']= {}
		self.stats['cov'] = {}
		self.stats['total_count'] = 0
		self.stats['total_span'] = 0
		self.stats['total_n50'] = 0
		self.stats['total_lengths'] = []
		self.stats['total_cov'] = {}
		self.stats['total_gc'] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}
		self.stats['cov_libs'] = []

		for contig_name in self.contigs:
			blob = self.contigs[contig_name]
			self.stats['total_count'] += 1
			self.stats['total_span'] += blob.length
			self.stats['total_lengths'].append(blob.length)
			self.stats['total_gc']['raw'].append(blob.gc)

			for blast_lib in self.blast_libs:
				
				taxid = blob.tax[blast_lib]['taxid']
				if not blast_lib in self.stats['count']:
					self.stats['count'][blast_lib] = {}
					self.stats['span'][blast_lib] = {}
					self.stats['lengths'][blast_lib] = {}
					self.stats['gc'][blast_lib] = {}
					self.stats['cov'][blast_lib] = {}
				self.stats['count'][blast_lib][taxid] = self.stats['count'][blast_lib].get(taxid, 0) + 1	
				self.stats['span'][blast_lib][taxid] = self.stats['span'][blast_lib].get(taxid, 0) + blob.length

				if not taxid in self.stats['gc'][blast_lib]:
					self.stats['gc'][blast_lib][taxid] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}
					self.stats['lengths'][blast_lib][taxid] = []	
				self.stats['gc'][blast_lib][taxid]['raw'].append(blob.gc) 
				self.stats['lengths'][blast_lib][taxid].append(blob.length)
				
				for cov_lib, cov in blob.covs.items():
					if not cov_lib in self.stats['cov'][blast_lib]:
						self.stats['total_cov'][cov_lib] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}
						self.stats['cov'][blast_lib][cov_lib]={}
					if not taxid in self.stats['cov'][blast_lib][cov_lib]:
						self.stats['cov'][blast_lib][cov_lib][taxid] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0} 
					self.stats['cov'][blast_lib][cov_lib][taxid]['raw'].append(cov)
			
			for cov_lib, cov in blob.covs.items():
				self.stats['total_cov'][cov_lib]['raw'].append(cov)

		for blast_lib in self.blast_libs:
			# calculate N50
			for tax, list_of_lengths in self.stats['lengths'][blast_lib].items():
				if not blast_lib in self.stats['n50']:
					self.stats['n50'][blast_lib] = {}
				self.stats['n50'][blast_lib][tax] = n50(list_of_lengths)
			self.stats['total_n50'] = n50(self.stats['total_lengths'])

			# calculate total gc mean/stdev
			for tax in self.stats['gc'][blast_lib]:
				self.stats['gc'][blast_lib][tax]['mean'] = "{0:.2f}".format(numpy.mean(self.stats['gc'][blast_lib][tax]['raw']))
				self.stats['gc'][blast_lib][tax]['stdev'] = "{0:.2f}".format(numpy.std(self.stats['gc'][blast_lib][tax]['raw']))

			# calculate total cov mean/stdev
			for cov_lib in self.stats['cov'][blast_lib]:
				self.stats['total_cov'][cov_lib]['mean'] = "{0:.2f}".format(numpy.mean(self.stats['total_cov'][cov_lib]['raw']))
				self.stats['total_cov'][cov_lib]['stdev'] = "{0:.2f}".format(numpy.std(self.stats['total_cov'][cov_lib]['raw']))
				
				# calculate tax-specific cov mean/stdev
				for tax in self.stats['cov'][blast_lib][cov_lib]:
					self.stats['cov'][blast_lib][cov_lib][tax]['mean'] = "{0:.2f}".format(numpy.mean(self.stats['cov'][blast_lib][cov_lib][tax]['raw']))
					self.stats['cov'][blast_lib][cov_lib][tax]['stdev'] = "{0:.2f}".format(numpy.std(self.stats['cov'][blast_lib][cov_lib][tax]['raw']))
		
		self.stats['total_gc']['mean'] = "{0:.2f}".format(numpy.mean(self.stats['total_gc']['raw']))
		self.stats['total_gc']['stdev'] = "{0:.2f}".format(numpy.std(self.stats['total_gc']['raw']))

	def writeOutput(self, version):
		
		'''
		Writes outputfiles:
		- stats.txt which contains tables with the stats calculated through getStats()
		- blobs.txt which contains the blobs
		'''

		header = '# makeblobs.py v{}\n# {} {}\n# {}\n'.format(version, time.strftime("%Y-%m-%d"), time.strftime("%H:%M:%S"), " ".join(sys.argv))
		
		# Init NCBI Taxonomy if necessary.
		taxonomy = NCBITaxa() if self.tax_label == 'scientific name' else None

		# Writing of stats.txt
		stats_fh = open(self.outfiles['stats'], 'w')
		stats_fh.write(header)
		stats_string = ''
		for blast_lib in self.blast_libs:
			stats_string += "\tTAX:{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:^10}".format(blast_lib, "contigs", "span", "N50", "GC")
			for cov_lib in self.cov_libs:
				stats_string += "\t{:<10}".format(cov_lib)
			stats_string += "\n\t" + (("-" * 10) + "\t") * (5 + len(self.cov_libs)) + "\n"
			for tax, score in sorted(self.stats['span'][blast_lib].items(), key=lambda x: x[1], reverse = True):
				taxon = formatTaxon(tax, tax_label=self.tax_label, taxonomy=taxonomy)
				stats_string += "\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}".format(taxon, self.stats['count'][blast_lib][tax], self.stats['span'][blast_lib][tax], self.stats['n50'][blast_lib][tax],  str(self.stats['gc'][blast_lib][tax]['mean']) + " SD:" + str(self.stats['gc'][blast_lib][tax]['stdev']))
				for cov_lib in self.cov_libs:
					stats_string += "\t{:>10}".format( str(self.stats['cov'][blast_lib][cov_lib][tax]['mean']) + " SD:" + str(self.stats['cov'][blast_lib][cov_lib][tax]['stdev']))
				stats_string += "\n"
			stats_string += "\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}".format("Total", self.stats['total_count'], self.stats['total_span'], self.stats['total_n50'],  str(self.stats['total_gc']['mean']) + " SD:" + str(self.stats['total_gc']['stdev']))
			for cov_lib in self.cov_libs:
				stats_string += "\t{:>10}".format( str(self.stats['total_cov'][cov_lib]['mean']) + " SD:" + str(self.stats['total_cov'][cov_lib]['stdev']))
			stats_string += "\n\n"
		stats_fh.write(stats_string)
		stats_fh.close()

		# Writing of blobs.txt
		blobs_fh = open(self.outfiles['blobs'], 'w')
		blobs_fh.write(header)
		blobs_fh.write("# contig_id\tlength\tgc\tcov\ttaxonomy\n")
		blobs_string = ''
		for contig_name in self.index:
			blob = self.contigs[contig_name]
			blobs_string += "{}\t{}\t{:.3f}".format(blob.name, blob.length, blob.gc)
			cov_string, tax_string = '\t','\t'
			for cov_lib, cov in blob.covs.items():
				cov_string += "{}={};".format(cov_lib, cov) 
			blobs_string += cov_string[0:-1]
			for blast_lib in blob.tax:
				taxon = formatTaxon(blob.tax[blast_lib]['taxid'], 
					tax_label=self.tax_label, taxonomy=taxonomy)
				tax_string += "{}={};".format(blast_lib, taxon)
			blobs_string += tax_string[0:-1]
			blobs_string += "\n"
		blobs_fh.write(blobs_string)
		blobs_fh.close()
