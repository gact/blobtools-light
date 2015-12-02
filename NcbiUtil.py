#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Utilities for NCBI databases and tools.'''

from Bio import Entrez
from collections import MutableSequence
from csv import DictReader
from csv import DictWriter
from datetime import datetime
import os
import re
from socket import error as SocketError
from subprocess import check_call
import tarfile
from tempfile import NamedTemporaryFile
from time import sleep
from urllib import urlretrieve
from urllib2 import HTTPError
from urllib2 import URLError
from urlparse import urljoin

import Util

################################################################################

blast_dbtypes = ('nucl', 'prot')

supported_blast_fields = ('qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 
	'sallseqid', 'sgi', 'sallgi', 'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 
	'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 
	'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 
	'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 
	'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 
	'sstrand', 'qcovs', 'qcovhsp')

default_blast_fields = ('qseqid', 'staxids', 'bitscore', 'evalue', 'pident',
	'sseqid', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send')

required_blast_fields = ('qseqid', 'staxids', 'bitscore', 'evalue')

supported_blast_tasks = ('megablast', 'dc-megablast', 'blastn')

default_blast_task = 'megablast'

################################################################################

def entrezEsearch(db, term, email, **kwargs):
	'''Wrapper function for Entrez esearch.'''

	faults = ('WarningList', 'ErrorList')
	max_attempts = 10
	result = None
	attempts = 0
	
	# Set Entrez email attribute.
	Entrez.email = email
	
	while attempts < max_attempts:
		try:
			request = Entrez.esearch(db, term, **kwargs)
			result = Entrez.read(request)
		except (HTTPError, RuntimeError, SocketError, URLError):
			sleep(60.0)
		else:
			break
	
	if result is None:
		raise RuntimeError("Entrez esearch failed after {!r} attempts".format(attempts))
	
	if any( fault in result for fault in faults ):
		msg = "check parameters"
		for fault in faults:
			if 'OutputMessage' in result[fault]:
				msg = '\n'.join( result[fault]['OutputMessage'] )
				break
		raise RuntimeError("Entrez esearch failed: {}".format(msg))
	
	return result

def filterBestBlastHits(full_blast_file, filt_blast_file, blast_fields=None):
	'''Filter BLAST tabular result file to retain only best BLAST hits.
	
	Where ties occur, all tied hits are retained.
	'''

	headings = readBlastTableHeadings(full_blast_file)

	if blast_fields is not None:
		fields = validateBlastFields(blast_fields)
		if headings is not None and ( len(headings) != len(fields) or 
			any( x != y for x, y in zip(headings, fields) ) ):
			raise ValueError("BLAST field mismatch")	
	elif headings is not None:
		fields = headings
	else:
		fields = None
	
	hits = getDictOfBestBlastHits(full_blast_file)
	
	# Get set of indices for lines containing best BLAST hits.
	indices = set()
	for q in hits:
		for t in hits[q]['taxa']:
			indices.update(hits[q]['taxa'][t])
	
	if fields is not None:
		writeBlastTableHeadings(fields, filt_blast_file)
		
	with open(filt_blast_file, 'a') as fout:
		
		writer = DictWriter(fout, fieldnames=fields, dialect='blast-tab')
				
		with open(full_blast_file, 'r') as fin:
			
			reader = DictReader(fin, fieldnames=fields, dialect='blast-tab')
			
			if headings is not None:
				next(reader)
			
			for i, row in enumerate(reader):
				if i in indices:
					writer.writerow(row)

def getDictOfBestBlastHits(blast_file):
	'''Get dictionary of best BLAST hits from a BLAST tabular result file.
	
	Where ties occur, all tied hits are retained. The keys of the returned dict
	are query sequence IDs, and each value is a dictionary with three keys: 
	'evalue', 'bitscore', 'taxa'. The E-value and bitscore are the best scores 
	found for the given query sequence, while each element of the 'taxa' dict
	contains the taxids of database sequences with the given E-value and bitscore, 
	each a key to a set of indices indicating the lines in the BLAST result file 
	containing matches between the given query and that taxon.
	'''

	hits = dict()
	
	# Read BLAST headings from BLAST result file.
	headings = readBlastTableHeadings(blast_file)
	
	# If BLAST headings found, set BLAST fields from these..
	if headings is not None:
		blast_fields = headings
	# ..otherwise assume required BLAST fields are present.
	else:
		blast_fields = required_blast_fields

	with open(blast_file, 'r') as fh:
				
		reader = DictReader(fh, fieldnames=blast_fields, dialect='blast-tab')
						
		if headings is not None:
			next(reader)
		
		for i, row in enumerate(reader):
			
			# Get BLAST field values.
			q, staxids, bitscore, evalue = [ row[k] for k in required_blast_fields ]
			
			# Validate BLAST field values.
			try:
				assert q is not None
				taxids = staxids.split(";")
				bitscore = float(bitscore)
				evalue = float(evalue)
			except (AssertionError, AttributeError, TypeError, ValueError):
				raise RuntimeError("BLAST results do not seem to be in the right"
					" format ('6 qseqid staxids bitscore evalue ... ')")
			
			# Validate taxids field.
			try:
				taxids = [ int(x) for x in taxids ]
			except ValueError:
				if len(taxids) != 1 or taxids[0] != 'N/A':
					raise ValueError("invalid NCBI Taxonomy ID: {!r}".format(staxids))
				continue
				
			# Compare this hit to the best previous hit(s), if any.
			if q in hits:
				
				# Get info for previous best hit(s).
				e, s = [ hits[q][k] for k in ('evalue', 'bitscore') ]
			
				# Continue to next hit unless this hit has an E-value or
				# bitscore equal to or better than the previous best for 
				# this contig. If this hit has scores equal to the best 
				# previous hit(s), add taxids and line index to hit info.
				if evalue > e:
					continue					
				elif evalue == e:
					if bitscore < s:
						continue
					elif bitscore == s:
						for taxid in taxids:
							hits[q]['taxa'].setdefault(taxid, set())
							hits[q]['taxa'][taxid].add(i)
						continue
							
			# This hit is first or best, so set new hit info.
			hits[q] = { 'evalue': evalue, 'bitscore': bitscore }
			hits[q]['taxa'] = { taxid: set([i]) for taxid in taxids }

	return hits

def readBlastTableHeadings(blast_file):
	'''Read header row of NCBI BLAST tabular result file.
	
	Returns None if no header row found.
	'''

	blast_headings = None
	
	with open(blast_file, 'rU') as fh:
		
		line = fh.readline().strip()
		
		if line.startswith('#'):
			
			fields = tuple( line[1:].strip().split('\t') )

			try:
				blast_headings = validateBlastFields(fields)
			except ValueError as e:
				pass
				
	return blast_headings

def validateBlastFields(blast_fields):
	'''Validate BLAST field names.'''
 
	for field in blast_fields:
		if field not in supported_blast_fields:
			raise ValueError("unknown BLAST field: {!r}".format(field)) 
	
	for i, required_field in enumerate(required_blast_fields):
		if blast_fields[i] != required_field:
			raise ValueError("BLAST field in column {!r} should be {!r}, not {!r}".format(i+1, required_field, blast_fields[i]))
	
	return tuple( x for x in blast_fields )

def writeBlastTableHeadings(blast_headings, blast_file):
	'''Write header row of NCBI BLAST tabular result file.'''
	with open(blast_file, 'w') as fh:
		fh.write("# {}\n".format('\t'.join(blast_headings)))

################################################################################

class NcbiBlastdb(object):
	'''NCBI BLAST database.'''

	@classmethod
	def check(cls, database_path, database_name):
		'''Check for presence of specified NCBI BLAST database in given location.
		
		Raises RuntimeError if not present.
		'''
		alias = NcbiBlastdbAlias.fromBlastdb(database_path, database_name)
		for database_part in alias.dblist:
			NcbiBlastdbPart.check(database_path, database_part)

	@classmethod
	def present(cls, database_path, database_name):
		'''Check if specified NCBI BLAST database is present in given location.
		
		Returns True if present, False otherwise.
		'''
		try:
			NcbiBlastdb.check(database_path, database_name)
			return True
		except RuntimeError:
			return False
			
	def __init__(self, database_name):
		'''Init NCBI BLAST database.'''
		self.name = database_name
		
	def download(self, database_path):
		'''Download NCBI BLAST database to the given location.'''
		
		# Assume that the first database part name is either the database name
		# or the database name appended with a two-digit numerical suffix. 
		zero_parts = ( NcbiBlastdbPart(self.name + '.00'), 
			NcbiBlastdbPart(self.name) )
			
		zero_part_found = False
		
		for zero_part in zero_parts:
		
			try:
				zero_part.download(database_path)
			except IOError:
				sleep(5.0)
			else:
				zero_part_found = True
				break 
	
		if not zero_part_found:
			raise RuntimeError("initial database part not found for NCBI BLAST database {!r}".format(self.name)) 
					
		alias_ext = NcbiBlastdbAlias.aliasFileExtensionFromDbtype(zero_part.dbtype)		
		alias_path = os.path.join(database_path, self.name + alias_ext)
		alias = NcbiBlastdbAlias(alias_path)
		
		for database_part in alias.dblist[1:]:
			dbpart = NcbiBlastdbPart(database_part)
			dbpart.download(database_path)
			
		NcbiBlastdb.check(database_path, self.name)


class NcbiBlastdbAlias(object):
	'''NCBI BLAST database alias.'''
	
	supported_attributes = ('TITLE', 'DBLIST', 'GILIST', 'LENGTH', 'NSEQ', 'MEMB_BIT')
	required_attributes = ('DBLIST',)
	
	@classmethod
	def aliasFileExtensionFromDbtype(cls, dbtype):
		'''Get alias file extension corresponding to the given BLAST database type.'''
		
		if dbtype == 'nucl':
			ext = '.nal'
		elif dbtype == 'prot':
			ext = '.pal'
		else:
			raise ValueError("NCBI BLAST database database type must be 'nucl' or 'prot'")
		
		return ext   
	
	@classmethod
	def dbtypeFromAliasFile(cls, alias_file):
		'''Get BLAST database type corresponding to the given alias file name.'''
		
		if alias_file.endswith('.nal'):
			dbtype = 'nucl'
		elif alias_file.endswith('.pal'):
			dbtype = 'prot'
		else:
			raise ValueError("invalid extension on NCBI BLAST database alias filename {!r}".format(alias_file))
		
		return dbtype
		
	@classmethod
	def dbtypeFromBlastdb(cls, database_path, database_name):
		'''Get database type of the given NCBI BLAST database.'''

		alias_paths = {
			'nucl': os.path.join(database_path, '{}.nal'.format(database_name) ),
			'prot': os.path.join(database_path, '{}.pal'.format(database_name) )
		}
		
		dbtypes = [ k for k in alias_paths if os.path.exists(alias_paths[k]) ]
		
		if len(dbtypes) == 1:
			dbtype = dbtypes[0]
		elif len(dbtypes) > 1:
			raise RuntimeError("multiple alias files found for NCBI BLAST database {!r}".format(database_name))
		else:
			raise RuntimeError("no alias file found for NCBI BLAST database {!r}".format(database_name))
			
		return dbtype
		
	@classmethod
	def fromBlastdb(cls, database_path, database_name):
		'''Load alias from given NCBI BLAST database.'''
		
		alias_paths = {
			'nucl': os.path.join(database_path, '{}.nal'.format(database_name) ),
			'prot': os.path.join(database_path, '{}.pal'.format(database_name) )
		}
		
		dbtypes = [ k for k in alias_paths if os.path.exists(alias_paths[k]) ]
				
		if len(dbtypes) == 1:
			alias = NcbiBlastdbAlias( alias_paths[ dbtypes[0] ] )
		elif len(dbtypes) > 1:
			raise RuntimeError("multiple alias files found for NCBI BLAST database {!r}".format(database_name))
		else:
			raise RuntimeError("no alias file found for NCBI BLAST database {!r}".format(database_name))	   

		return alias		

	@classmethod
	def fromAliasFile(cls, alias_file):
		'''Load alias from given alias file.'''
		return NcbiBlastdbAlias(alias_file)
		
	def __init__(self, *args, **kwargs):
		'''Init NCBI BLAST database alias.'''
		if len(args) == 1:
			self.read(args[0])
		elif len(args) > 1:
			raise RuntimeError("too many arguments passed to NcbiBlastdbAlias()")
		for key, value in kwargs.items():
			setattr(self, key, value)
			
	def __setattr__(self, name, value):
		'''Set attribute of NCBI BLAST database alias.'''
		
		NOM, nom = name.upper(), name.lower()
	   
		if NOM in NcbiBlastdbAlias.supported_attributes:
		
			if NOM in ('TITLE', 'GILIST'):
				value = str(value)
			elif NOM in ('LENGTH', 'NSEQ'):
				value = int(value)
			elif NOM == 'DBLIST':
				value = NcbiBlastdbList(value)
				value.validate()
			elif value == 'MEMB_BIT':
				if value not in ('0', '1'):
					raise ValueError("NCBI BLAST database alias attribute 'MEMB_BIT' must be 0 or 1")
				value = int(value)
				
		elif nom == 'dbtype':
		
			if value.lower() not in blast_dbtypes:
				raise ValueError("NCBI BLAST database alias dbtype must be 'nucl' or 'prot'")
				
		else:
			raise ValueError("unsupported NCBI BLAST database alias attribute: {!r}".format(name))			
						   
		self.__dict__[ name.lower() ] = value
		
	def read(self, alias_file):
		'''Read alias attributes from given alias file.'''
		
		for attr in NcbiBlastdbAlias.supported_attributes:
			try:
				delattr(self, attr)
			except AttributeError:
				pass				
		
		self.dbtype = NcbiBlastdbAlias.dbtypeFromAliasFile(alias_file)
		
		with open(alias_file, 'r') as f:
		
			attr_regex = re.compile('^([^\s]+)\s+(.+)$')

			for line in f:
				
				if line.startswith('#'):
					continue
				
				line = line.rstrip()
				
				if line == '':
					continue
				
				valid_alias_attr = False

				m = attr_regex.match(line)
				
				if m is None:
					raise ValueError("failed to parse NCBI BLAST database alias file {!r}".format(alias_file))
				
				name = m.group(1)
				value = m.group(2)
			 
				NOM, nom = name.upper(), name.lower()
					
				if hasattr(self, nom):
					raise ValueError("duplicate key {!r} found in NCBI BLAST database alias file {!r}".format(NOM, alias_file))
				
				setattr(self, nom, value)
			
		for required in NcbiBlastdbAlias.required_attributes + ('DBTYPE',):
			if not hasattr(self, required.lower()):
				raise ValueError("required attributes not found in NCBI BLAST database alias file {!r}".format(alias_file))
					
	def write(self, alias_file):
		'''Write alias attributes to given alias file.'''
		
		dbtype = NcbiBlastdbAlias.dbtypeFromAliasFile(alias_file)
					
		if not hasattr(self, 'dbtype'):
			self.dbtype = dbtype
		elif dbtype != self.dbtype:
			raise ValueError("NCBI BLAST database type mismatch")
		
		for required in NcbiBlastdbAlias.required_attributes + ('DBTYPE',):
			if not hasattr(self, required.lower()):
				raise ValueError("cannot write NCBI BLAST database alias without attribute {!r}".format(required))

		with open(alias_file, 'w') as f:

			f.write("#\n")
			f.write("# Alias file created {}\n".format( 
				datetime.now().strftime("%m/%d/%Y %H:%M:%S")))
			f.write("#\n")

			for name in NcbiBlastdbAlias.supported_attributes:
				
				NOM, nom = name.upper(), name.lower()
				
				if hasattr(self, nom):
					
					value = getattr(self, nom)
					
					if NOM == 'DBLIST':
						value.validate()
					elif NOM == 'GILIST':
						value = NcbiGilist.resolveGilistFile(value, 
							alias_dir=os.path.dirname(alias_file), dbtype=dbtype)
						
					f.write("{} {}\n".format(NOM, str(value)))
		
				
class NcbiBlastdbList(MutableSequence):
	'''NCBI BLAST database list.'''

	@classmethod
	def fromString(cls, string):
		'''Get NCBI BLAST database list from the given string.'''
		
		if not isinstance(string, basestring):
			raise ValueError("NcbiBlastdbList::fromString argument must be of string type")
		
		dbpart = re.compile("""(?:["'][^"']+["'])|(?:[^\s]+)""")
		
		iterable = dbpart.findall(string)
		
		return NcbiBlastdbList.fromIterable(iterable)

	@classmethod
	def fromIterable(cls, iterable):
		'''Get NCBI BLAST database list from the given iterable.'''
		
		if isinstance(iterable, basestring):
			raise ValueError("NcbiBlastdbList::fromIterable argument cannot be of string type")
		
		quotes = ('"', "'")
				
		obj = NcbiBlastdbList()
		
		for element in iterable:
			
			if not isinstance(element, basestring):
				raise ValueError("NcbiBlastdbList elements must be of string type")
			
			for quote in quotes:
				if element.startswith(quote) and element.endswith(quote):
					element = element[1:-1]
					break
						
			if any( quote in element for quote in quotes ):
				raise ValueError("NcbiBlastdbList elements must not contain internal quotes")
		  
			obj.append(element)
				 
		return obj
		
	def __init__(self, dblist=None):
		'''Init NCBI BLAST database list.'''
		
		if dblist is None:
			dblist = list()
		elif isinstance(dblist, basestring):
			dblist = NcbiBlastdbList.fromString(dblist)
		elif hasattr(dblist, '__iter__'):
			dblist = NcbiBlastdbList.fromIterable(dblist)
		else:
			raise ValueError("NcbiBlastdbList must be a string or iterable of strings")
		
		self.dblist = list(dblist)
		
	def __add__(self, other):
		'''Concatenate this list of NCBI BLAST database parts with another.'''
		other = NcbiBlastdbList.fromIterable(other)
		result = NcbiBlastdbList.fromIterable(self.dblist + other.dblist)
		return result

	def __delitem__(self, key):
		'''Delete NCBI BLAST database parts from this list.'''
		del self.dblist[key]

	def __getitem__(self, key):
		'''Get NCBI BLAST database parts from this list.'''
		return self.dblist[key]

	def __iadd__(self, other):
		'''Append NCBI BLAST database parts to this list.'''
		other = NcbiBlastdbList.fromIterable(other)
		self = NcbiBlastdbList.fromIterable(self.dblist + other.dblist)
		return self
			
	def __len__(self):
		'''Get length of this NCBI BLAST database list.'''
		return len(self.dblist)

	def __setitem__(self, key, value):
		'''Set NCBI BLAST database parts in this list.'''
		if isinstance(key, int):
			self.dblist[key] = NcbiBlastdbList.fromString(value)
		elif isinstance(key, slice):
			self.dblist[key] = NcbiBlastdbList(value)
		else:
			raise ValueError("cannot set item of NcbiBlastdbList with key of type {!r}".format(key.__class__.__name__))

	def __str__(self):
		'''Convert NCBI BLAST database list to string.'''
		spaces = (' ', '\t')
		dblist = self.dblist
		if any( space in dbpart for dbpart in dblist for space in spaces ):
			dblist = [ '''"{}"'''.format(dbpart) for dbpart in dblist ]
		return ' '.join(dblist)
		
	def insert(self, index, element):
		'''Insert NCBI BLAST database part in this list.'''
		if not isinstance(element, basestring):
			raise ValueError("NcbiBlastdbList elements must be of string type")
		self.dblist.insert(index, element)
	
	def validate(self):
		'''Validate NCBI BLAST database list.'''
		
		if not hasattr(self, 'dblist') or len(self.dblist) == 0:
			raise RuntimeError("NcbiBlastdbList is empty")
			
		if len(self.dblist) > 1:
			
			part_names = [ NcbiBlastdbPart.parseName(x) for x in self.dblist ]
			part_nums = [ int(part_num) for database_name, part_num in part_names ]
			
			if any( i != num for i, num in enumerate(part_nums) ):
				raise RuntimeError("NcbiBlastdbList elements are in an unexpected order")
		
		
class NcbiBlastdbPart(object):
	'''NCBI BLAST database part.'''
 
	@classmethod
	def checkName(cls, database_part_name):
		'''Check validity of NCBI BLAST database part name.
		
		Raises ValueError if database part name appears to be invalid.
		'''
		cls.parseName(database_part_name)
	
	@classmethod
	def parseName(cls, database_part_name):
		'''Parse name of NCBI BLAST database part.
		
		Raises ValueError if database part name appears to be invalid.
		'''
		
		database_part_pattern = re.compile('^([^\s.]+)(?:[.](\d+))?$')
		
		m = database_part_pattern.match(database_part_name)
		
		try:
			database_name = m.group(1)
			part_num = m.group(2)
		except AttributeError:
			raise ValueError("invalid database part name: {!r}".format(database_part_name))
		
		if part_num is not None and ( len(part_num) < 2 or not part_num.isdigit() ):
			raise ValueError("invalid number on database part: {!r}".format(database_part_name))
			
		return database_name, part_num
		
	@classmethod
	def check(cls, database_path, database_part):
		'''Check for presence of specified NCBI BLAST database part in given location.
		
		Raises RuntimeError if not present.
		'''

		database_name, part_num = NcbiBlastdbPart.parseName(database_part)

		if not os.path.exists(database_path):
			raise RuntimeError("NCBI BLAST database directory not found: {!r}".format(database_path))
	
		dbtype = NcbiBlastdbAlias.dbtypeFromBlastdb(database_path, database_name)
		alias_ext = NcbiBlastdbAlias.aliasFileExtensionFromDbtype(dbtype)
		alias_path = os.path.join(database_path, database_name + alias_ext)
	  
		dbext = {
			'nucl': ('.nhd', '.nhi', '.nhr', '.nin', '.nnd', '.nni', '.nog', '.nsd', '.nsi', '.nsq'),
			'prot': ('.phd', '.phi', '.phr', '.pin', '.pnd', '.pni', '.pog', '.psd', '.psi', '.psq')
		}
	
		dbfiles = [ os.path.join(database_path, database_part + ext) 
			for ext in dbext[dbtype] ]
			
		if any( not os.path.exists(dbfile) for dbfile in dbfiles ):
			raise RuntimeError("NCBI BLAST database part files not found: {!r}".format(database_part))
				
		alias_mtime = os.path.getmtime(alias_path)
		
		# If any database files have been modified more recently
		# than the database alias file, they need to be updated.
		if any( os.path.getmtime(dbfile) > alias_mtime for dbfile in dbfiles ):
			raise RuntimeError("NCBI BLAST database files changed since last modification of alias file ({!r})".format(alias_file))
				
	@classmethod
	def present(cls, database_path, database_part):
		'''Check if specified NCBI BLAST database part is present in given location.
		
		Returns True if present, False otherwise.
		'''
		try:
			NcbiBlastdbPart.check(database_path, database_part)
			return True
		except RuntimeError:
			return False

	@property
	def name(self):
		'''Name of NCBI BLAST database part.'''
		name = self.dbname
		if self.dbpart is not None: name = '{}.{}'.format(name, self.dbpart)
		return name
		
	def __init__(self, database_part):
		'''Init NCBI BLAST database part.'''		
		database_name, part_num = NcbiBlastdbPart.parseName(database_part)
		self.dbname = database_name
		self.dbpart = part_num
	
	def download(self, database_path):
		'''Download NCBI BLAST database part to the given location.'''
				
		if not os.path.exists(database_path):
			os.makedirs(database_path)
		elif not os.path.isdir(database_path):
			raise ValueError("cannot download NCBI BLAST database part - destination is not a directory")
							
		blastdb_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/'
		
		# Download NCBI BLAST database part.
		blastdb_part_file = '{}.tar.gz'.format(self.name)
		blastdb_part_url = urljoin(blastdb_url, blastdb_part_file)
		blastdb_part_path = os.path.join(database_path, blastdb_part_file)
		urlretrieve(blastdb_part_url, filename=blastdb_part_path)
		
		# Download checksum file for NCBI BLAST database part.
		blastdb_sum_file =  '{}.md5'.format(blastdb_part_file)
		blastdb_sum_url = urljoin(blastdb_url, blastdb_sum_file)
		blastdb_sum_path = os.path.join(database_path, blastdb_sum_file)
		urlretrieve(blastdb_sum_url, filename=blastdb_sum_path)
		
		md5file_regex = re.compile("^([^\s]+)\s+(.+)$")
		
		with open(blastdb_sum_path, 'r') as f:
			md5file_content = f.read()
			m = md5file_regex.match(md5file_content)
		
			try:
				expected_checksum = m.group(1)
				expected_filename = m.group(2)
			except AttributeError:
				raise RuntimeError("failed to parse NCBI BLAST database checksum file: {!r}".format(blastdb_sum_path))  

		blastdb_part_checksum = Util.fileMD5(blastdb_part_path)
	
		if blastdb_part_file != expected_filename:
			raise RuntimeError("NCBI BLAST database archive filename mismatch: {!r}".format(blastdb_part_file))
			
		if blastdb_part_checksum != expected_checksum:
			raise RuntimeError("NCBI BLAST database archive checksum mismatch: {!r}".format(blastdb_part_checksum))
		
		os.remove(blastdb_sum_path)
		
		alias_file = None
		
		with tarfile.open(blastdb_part_path) as tar:
		
			blastdb_part_members = [ x.name for x in tar.getmembers() ]
			
			for blastdb_part_member in blastdb_part_members:
			
				if blastdb_part_member.endswith('.nal'):
					alias_file = blastdb_part_member
					self.dbtype = 'nucl'
				elif blastdb_part_member.endswith('.pal'):
					alias_file = blastdb_part_member
					self.dbtype = 'prot'
			
			if alias_file is None:
				raise RuntimeError("no alias file found for NCBI BLAST database part: {!r}".format(self.name))

		with tarfile.open(blastdb_part_path) as tar:
			tar.extractall(database_path)

		# Ensure alias file has later modified time than all other NCBI BLAST database
		# files. This will be used later to check if database needs to be updated.
		alias_path = os.path.join(database_path, alias_file)
		with open(alias_path, 'a'):
			os.utime(alias_path, None)
		
		NcbiBlastdbPart.check(database_path, self.name)
		
		os.remove(blastdb_part_path)


class NcbiGilist(MutableSequence):
	'''NCBI GI list.'''

	@classmethod
	def dbtypeFromGilistFilename(cls, gilist_file):
		'''Get BLAST database type corresponding to the given NCBI GI list file name.'''
		
		if gilist_file.endswith('.n.gil'):
			dbtype = 'nucl'
		elif gilist_file.endswith('.p.gil'):
			dbtype = 'prot'
		else:
			dbtype = None
		
		return dbtype 
		
	@classmethod
	def gilistFileExtensionFromDbtype(cls, dbtype=None):
		'''Get NCBI GI list file extension corresponding to the given BLAST database type.'''
		
		if dbtype == 'nucl':
			ext = '.n.gil'
		elif dbtype == 'prot':
			ext = '.p.gil'
		elif dbtype is None:
			ext = '.gil'
		else:
			raise ValueError("NCBI BLAST database database type must be 'nucl' or 'prot'")
		
		return ext		 

	@classmethod
	def fromTaxid(cls, taxid, dbtype, email):
		'''Get NCBI GI list from given NCBI Taxonomy ID.'''
		
		taxid = str(taxid)
		
		if not taxid.isdigit():
			raise ValueError("invalid NCBI Taxonomy ID: {!r}".format(taxid))		
		
		if dbtype not in blast_dbtypes:
			raise ValueError("NCBI BLAST database database type must be 'nucl' or 'prot'")
		
		db = 'nuccore' if dbtype == 'nucl' else 'protein'
		
		obj = NcbiGilist()

		term = "txid{}[Organism:exp]".format(taxid)
		
		result = entrezEsearch(db, term, email, retstart=0, retmax=1)
		
		count = int(result['Count'])
		
		batch_size = 1024
		batch_offset = 0
		
		for batch_offset in range(0, count, batch_size):
			result = entrezEsearch(db, term, email, rettype='gi', 
			  retmode='text', retstart=batch_offset, retmax=batch_size)
			obj += result['IdList']
			sleep(1.0)

		return obj

	@classmethod
	def fromIterable(cls, iterable):
		'''Get NCBI GI list from the given iterable.'''
		
		obj = NcbiGilist()
				
		if isinstance(iterable, basestring):
			raise ValueError("cannot set NcbiGilist from string") 
				
		for element in iterable:
			
			try:
				gi = int(element)
			except ValueError:
				raise ValueError("invalid GI number: {!r}".format(element))
			
			obj.append(gi)
				 
		return obj
 
	@classmethod
	def resolveGilistFile(cls, gilist_file, alias_dir=os.getcwd(), dbtype=None):
		'''Get path to the specified GI list file, relative to the given alias file.
		
		Get the resolved path of the GI list file, relative to the directory 
		containing the given alias file. Additionally, ensure that the GI list
		file is in the binary format taken as input by NCBI BLAST.
		'''
				
		gilist_path = os.path.normpath( os.path.join(alias_dir, gilist_file) )
						
		if dbtype is not None and dbtype not in blast_dbtypes:
			raise ValueError("NCBI BLAST database database type must be 'nucl' or 'prot'")		
				
		if gilist_path.endswith('.gil'):
			
			gilist_dbtype = NcbiGilist.dbtypeFromGilistFilename(gilist_path)
			
			if dbtype is not None and gilist_dbtype is not None:
				if gilist_dbtype != dbtype:
					raise ValueError("NCBI BLAST database type mismatch")
		else:
		
			gilist = NcbiGilist()
		
			gilist.read(gilist_path)
				
			gilist_prefix, old_ext = os.path.splitext(gilist_path)
			
			new_ext = NcbiGilist.gilistFileExtensionFromDbtype(dbtype)
			
			gilist_path = '{}{}'.format(gilist_prefix, new_ext)
			
			gilist.write(gilist_path)
				
		gilist_path = os.path.relpath(gilist_path, alias_dir)
				
		return gilist_path
		
	def __init__(self, gilist=None):
		'''Init NCBI GI list.'''
	
		if gilist is None:
			self.gilist = list()
		else:
			NcbiGilist.fromIterable(gilist)
		
	def __add__(self, other):
		'''Concatenate this list of GI numbers with another.'''
		other = NcbiGilist.fromIterable(other)
		result = NcbiGilist.fromIterable(self.gilist + other.gilist)
		return result

	def __delitem__(self, key):
		'''Delete GI numbers from this list.'''
		del self.gilist[key]

	def __getitem__(self, key):
		'''Get GI numbers from this list.'''
		return self.gilist[key]

	def __iadd__(self, other):
		'''Append GI numbers to this list.'''
		other = NcbiGilist.fromIterable(other)
		self = NcbiGilist.fromIterable(self.gilist + other.gilist)
		return self
			
	def __len__(self):
		'''Get length of this GI list.'''
		return len(self.gilist)

	def __setitem__(self, key, value):
		'''Set GI numbers in this list.'''
		if isinstance(key, int):
			try:
				self.gilist[key] = int(value)
			except ValueError:
				raise ValueError("invalid GI number: {!r}".format(value))
		elif isinstance(key, slice):
			self.gilist[key] = NcbiGilist.fromIterable(value)
		else:
			raise ValueError("cannot set item of NcbiGilist with key of type {!r}".format(key.__class__.__name__))

	def __str__(self):
		'''Convert GI list to string.'''
		return str(self.gilist)
		
	def insert(self, index, value):
		'''Insert GI number in this list.'''
		try:
			self.gilist.insert( index, int(value) )
		except ValueError:
			raise ValueError("invalid GI number: {!r}".format(value))
	
	def read(self, gilist_file):
		'''Read GI list from given file.'''
		
		self.gilist = list()
		
		with open(gilist_file, 'r') as f:
			
			for line in f:
				
				line = line.rstrip()
				
				if line == '':
					continue
				
				try:
					gi = int(line)
				except ValueError:
					raise ValueError("invalid GI number ({!r}) in text file {!r}".format(line, gilist_file))
			
				self.gilist.append(gi)  

	def write(self, gilist_file):
		'''Write GI list to given file.'''
		
		if gilist_file.endswith('.n.gil') or gilist_file.endswith('.p.gil'):
					
			with NamedTemporaryFile(mode='w') as fh:
						
				for gi in self.gilist:
					fh.write("{}\n".format(str(gi)))
			
				with open(os.devnull, 'w') as nul:
					check_call([ 'blastdb_aliastool', '-gi_file_in', fh.name, 
						'-gi_file_out', gilist_file], stdout=nul, stderr=nul)
		else:
			
			with open(gilist_file, 'w') as f:
			
				for gi in self.gilist:
					f.write("{}\n".format(str(gi)))

################################################################################
