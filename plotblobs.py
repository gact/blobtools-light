#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File	: plotblobs.py
Version : 0.1
Author	: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs	: ?
To do	:	- Add read proportion histogram 
			- legend placement at bottom
			- total reads mapped per taxon
		
"""

# # # # # 
# MODULES										
# # # # # 

from __future__ import division
import argparse
from collections import OrderedDict
import math as math
from matplotlib import cm
from matplotlib.lines import Line2D
import matplotlib as mat
import matplotlib.pyplot as plt
from matplotlib.pyplot import NullLocator
from matplotlib.ticker import NullFormatter
import numpy as np
import os
import sys

from MiscFunctions import n50
from Util import formatBaseCount
from Util import special_plot_colors

mat.rcParams.update({'font.size': 30})
mat.rcParams['xtick.major.pad']='8'
mat.rcParams['ytick.major.pad']='8'
mat.rcParams['lines.antialiased']=True

def set_canvas():
	left, width = 0.1, 0.60
	bottom, height = 0.1, 0.60
	bottom_h = left_h = left+width+0.02
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	rect_legend = [left_h, bottom_h, 0.2, 0.2]
	return rect_scatter, rect_histx, rect_histy, rect_legend

def set_format_scatterplot(axScatter, max_cov=0):
	axScatter.set_xlabel("GC proportion", fontsize=35)
	axScatter.set_ylabel("Coverage", fontsize=35)
	axScatter.grid(True, which="major", lw=2., color=white, linestyle='-') 
	axScatter.set_axisbelow(True)
	axScatter.set_xlim( (0, 1) )
	axScatter.set_ylim( (0.01, max_cov+1000) ) # This sets the max-Coverage so that all libraries + sum are at the same scale
	axScatter.xaxis.labelpad = 20
	axScatter.xaxis.labelpad = 20
	return axScatter

def set_format_hist_x(axHistx, axScatter):
	axHistx.set_xlim( axScatter.get_xlim() )
	axHistx.grid(True, which="major", lw=2., color= white, linestyle='-')
	axHistx.xaxis.set_major_formatter(nullFmt) # no labels since redundant
	axHistx.set_axisbelow(True)
	axHistx.yaxis.labelpad = 20
	return axHistx

def set_format_hist_y(axHisty, axScatter):
	axHisty.set_yscale('log')
	axHisty.yaxis.set_major_formatter(nullFmt) # no labels since redundant
	axHisty.set_ylim( axScatter.get_ylim() )
	axHisty.grid(True, which="major", lw=2., color= white, linestyle='-')
	axHisty.set_axisbelow(True)
	axHisty.xaxis.labelpad = 20
	return axHisty

def plot_ref_legend(axScatter, fontsize=24):
	s = 15
	# markersize in scatter is in "points^2", markersize in Line2D is in "points" ... that's why we need math.sqrt()
	ref_1 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(1000/15),  markerfacecolor=grey))
	ref_2 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(5000/15), markerfacecolor=grey))
	ref_3 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(10000/15), markerfacecolor=grey))
	axScatter.legend([ref_1,ref_2,ref_3], ["1,000nt", "5,000nt", "10,000nt"], numpoints=1, loc = 4, fontsize=fontsize)

def plot(data, cov_data, outfile, title, multi_plot=False, ignore_contig_len=False,
	hist_span=True, min_cov=None, max_cov=None, tax_groups=None, color_dict=None, 
	label_dict=None, fig_format='png', classifier=None):
	""" Plotting function which gets masked data and plots to outfile"""
	
	legend_fontsize = 24

	assert tax_groups is not None, "'tax_groups' not found"
	assert color_dict is not None, "'color_dict' not found"
	assert label_dict is not None, "'label_dict' not found"
	
	# If not specified explicitly, set coverage bounds from data.
	if min_cov is None:
		min_cov = np.amin(cov_data)
	if max_cov is None:
		max_cov = np.amax(cov_data)
	
	rect_scatter, rect_histx, rect_histy, rect_legend = set_canvas()
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Setting up plots and axes
	plt.figure(1, figsize=(35,35), dpi=400)

	axScatter = plt.axes(rect_scatter, axisbg=background_grey, yscale = 'log')
	axScatter = set_format_scatterplot(axScatter, max_cov=max_cov)
	axHistx = plt.axes(rect_histx, axisbg=background_grey)
	axHistx = set_format_hist_x(axHistx, axScatter)
	axHisty = plt.axes(rect_histy, axisbg=background_grey)
	axHisty = set_format_hist_y(axHisty, axScatter)
	axScatter.yaxis.get_major_ticks()[0].label1.set_visible(False)
	axScatter.yaxis.get_major_ticks()[1].label1.set_visible(False)
	if (title):
		plt.suptitle(title, fontsize=35, verticalalignment='top')
	#plt.suptitle(out_file, fontsize=25, verticalalignment='bottom')
	
	axLegend = plt.axes(rect_legend, axisbg=white)
	axLegend.xaxis.set_major_locator(nullLoc)
	axLegend.xaxis.set_major_formatter(nullFmt)  # no labels on axes
	axLegend.yaxis.set_major_locator(nullLoc)
	axLegend.yaxis.set_major_formatter(nullFmt)  # no labels on axes
	#
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Setting bins for histograms
	top_bins = np.arange(0, 1.01, 0.01)
	right_bins = np.logspace(-2, (int(math.log(max_cov)) + 1), 200, base=10.0)

	# empty handles for big legend
	legend_handles = []
	legend_labels = []

	# change file name if span (should be in input parsing function)
	if hist_span:
		outfile += ".hist_span"
	else:
		outfile += ".hist_count"

	# counter necessary for multiplot so that PNGs are in order when sorted by name
	i = 0

	# initiate variables for plotting
	s, lw, alpha, color = 0, 0, 0, ''
	
	# Maybe make an STDOUT printing func?
	print "[STATUS] Plotting : " + outfile

	# for each tax-group in order
	for g in tax_groups:
	
		i += 1

		# get indices for those rows in data where the taxon is in this group
		index_for_group = np.in1d(data[:,3].astype(str), tax_groups[g])
		# count of contigs ... total number of contigs comes from previous step?
		number_of_contigs_for_group = np.sum(index_for_group)
		
		# uses number_of_contigs for checking whether plotting should be carried out ... maybe there is a better place for this ...
		if number_of_contigs_for_group == 0:
			pass
		else:
			
			# create np_arrays for length, gc and cov for all contigs in tax-group 
			len_array = data[index_for_group][:,1].astype(int)
			gc_array = data[index_for_group][:,2].astype(float)
			cov_array = cov_data[index_for_group].astype(float)
			# sets label from label_dict
			label = label_dict[g]

			# another status message
			print "\t" + label 
			s_array = []
			# ignore contig length ... maybe do this in input and set these params for plotting there ...
			if (ignore_contig_len):
				if g == 'no-hit':
					s, lw, alpha = 15, 0.5, 0.5
				else:
					s, lw, alpha = 65, 0.5, 1
				s_array = [s for contig_length in len_array]
			else:
				if g == 'no-hit':
					s, lw, alpha = 15, 0.5, 0.5
				else:
					s, lw, alpha = 15, 0.5, 1
				# these are the sizes for plotting with contig sizes
				s_array = [contig_length/s for contig_length in len_array]
			# sets colour from color_dict
			color = color_dict[g]
						
			# making copies of gc/cov_array
			gc_hist_array = gc_array
			cov_hist_array = cov_array

			#######
			# if hist span ... 
			#	 make a new array ...
			#	 add to the array : (gc * len/1000) - 1
			# substitute old array with new array

			# set histogram labels depending on type ... can be set before ... 
			weights_array = len_array/1000 
			if (hist_span):
				axHistx.set_ylabel("Span (kb)")
				axHisty.set_xlabel("Span (kb)", rotation='horizontal')
			else:
				axHistx.set_ylabel("Count")
				axHisty.set_xlabel("Count", rotation='horizontal')
		
			# this should be set before ... or after? ... but only once 
			for xtick in axHisty.get_xticklabels(): # rotate text for ticks in cov histogram 
				xtick.set_rotation(270)

			# add text to legend ... label was build before ... could be a function
			legend_handles.append(Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=24, markerfacecolor=color))
			legend_labels.append(label)
			if (number_of_contigs_for_group):
				if (hist_span):
					axHistx.hist(gc_hist_array, weights=weights_array , color = color, bins = top_bins, histtype='step', lw = 3)
					axHisty.hist(cov_hist_array, weights=weights_array , color = color, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
				else:			
					axHistx.hist(gc_hist_array, color = color, bins = top_bins, histtype='step', lw = 3)
					axHisty.hist(cov_hist_array , color = color, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
		
			# If classifier specified, draw classifier regions and boundary.
			if classifier is not None:
			
				# Get GC-content range.
				gc_range = np.linspace( *axScatter.get_xlim() )
				
				# Get log-scaled range of coverage depth.
				cov_range = np.logspace( *np.log10(axScatter.get_ylim()).tolist() )
				
				# Get mesh grid of plot space.
				X, Y = np.meshgrid(gc_range, cov_range)
				
				# Get classifier value at each point in mesh grid.
				Z = classifier.decision_function( np.c_[ X.ravel(), Y.ravel() ] )
				Z = Z.reshape(X.shape)
				
				# Set contour levels (positive for target, negative for foreign).
				contour_levels = (np.amin(Z), 0, np.amax(Z))
				
				# Draw classifier regions and boundary in background.
				axScatter.contourf(X, Y, Z, levels=contour_levels, 
					cmap=plt.cm.coolwarm_r, alpha=0.1, zorder=0)
				axScatter.contour(X, Y, Z, levels=[0], colors='black', 
					linestyles='dashed', lw=0.5, zorder=0.1)
				
			axScatter.scatter(gc_array, cov_array, color = color, s = s_array, lw = lw, alpha=alpha, edgecolor=black, label=label)
		
			axLegend.axis('off')

			if (multi_plot): # MULTI-PLOT!!!
				axLegend.legend(legend_handles, legend_labels, loc=6, numpoints=1, fontsize=legend_fontsize, frameon=True)
				plot_ref_legend(axScatter, fontsize=legend_fontsize)
				#plt.savefig(outfile + "." + str(i) + "_"+ g.replace("/","") + "." + fig_format, format=fig_format)
				plt.savefig(outfile + "." + str(i) + "_"+ g.replace("/","") + "." + fig_format, format=fig_format)
	
	if ignore_contig_len:
		pass
	else: # print scale-legend
		plot_ref_legend(axScatter, fontsize=legend_fontsize)

	axLegend.legend(legend_handles, legend_labels, numpoints=1, fontsize=legend_fontsize, frameon=True, loc=6 )		
	sys.stdout.write("Saving file " + outfile)
	plt.savefig(outfile + "." + fig_format, format=fig_format)
	plt.close()
	print " [Done]\n" 

def getInput():

	parser = argparse.ArgumentParser(
		prog='plotblobs.py',
		usage = '%(prog)s infile [-p] [-f] [-t] [-e] [-n] [-o] [-l] [-c] [-s] [-h]',
		add_help=True)
	parser.add_argument('i', metavar = 'infile', help='Input file (blobplot.txt)')
	parser.add_argument('-p', metavar = 'max_taxa_plot', default=7, type = int, help='Maximum number of taxa to plot (Default = 7)')
	#parser.add_argument('-t', metavar = 'tax_level', default=2, type = int, help='Taxonomic level on which to plot.  Species = 0, Order = 1, Phylum = 2, Superkingdom = 3 (Default = 2)')
	#parser.add_argument('-e', metavar = 'eval_cutoffs' , default=[1.0], type = float, nargs='+', help='Set maximal e-value(s) (Default = 1.0)') 
	parser.add_argument('-c', metavar = 'len_cutoffs' , default=[100], type = int, nargs='+', help='Set minium contig length(s) (Default = 100)') 
	parser.add_argument('-s', action='store_true' , help='Ignore contig length for plotting.') 
	parser.add_argument('-n', action='store_true', help='Hides "no-hit" contigs') 
	parser.add_argument('-o', metavar ='out_prefix', default='' , help='Set output file prefix.') 
	parser.add_argument('-m', action='store_true' , help='Multi-plot. Print PNG after each tax-addition.') 
	parser.add_argument('-sort', action='store_false' , help='Sort by number of contigs per taxon (Default: Sort by span of taxon)') 
	parser.add_argument('-hist', action='store_false' , help='Make histograms based on contig counts. (Default: Span-Weighted histograms)') 
	parser.add_argument('-title', metavar='title' , help='Add title to plot') 
	parser.add_argument('-v', action='version', version='%(prog)s version 0.1')
	args = parser.parse_args()
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	out_prefix, sort_by_span, multi_plot, title = args.o, args.sort, args.m, args.title

	len_cutoffs, ignore_contig_len, hide_not_annotated, hist_span = args.c, args.s, args.n, args.hist

	infile = args.i

	if not os.path.exists(infile):
		parser.exit("[ERROR] : Input file {} does not exist!".format(infile))

	max_taxa_plot = args.p
	
	if (max_taxa_plot < 1):
		parser.exit("[ERROR] : 'max_taxa_plot' must be a positive integer!")
	
	return infile, out_prefix, max_taxa_plot, len_cutoffs, ignore_contig_len, sort_by_span, multi_plot, hist_span, title 

def parseInfile(data):
	
	contig_data_list = []
	cov_data_dict = {}
	tax_dict = {}
	with open(infile) as fh:
		for line in fh:
			if line.startswith('#'):
				pass
			else:
				line_data = line.rstrip("\n").split("\t")
				contig_id = line_data[0]
				length = int(line_data[1])
				gc = float(line_data[2])
				cov_dict = dict(string.split('=') for string in line_data[3].split(";"))
				cov_dict = {k: (float(v) if float(v) >= 0.1 else 0.1) for k, v in cov_dict.items()} # set coverages below 0.1 to 0.1
				blast_dict = dict(string.split('=') for string in line_data[4].split(";"))
				blast_dict = {k: v.split(":")[0] for k, v in blast_dict.items()} # removing any text after colon (if present), since not needed
				tax = blast_dict['tax']
				
				contig_data_list.append([(contig_id), (length), (gc), (tax)])
				for cov_lib, cov in cov_dict.items():
					if not cov_lib in cov_data_dict:
						cov_data_dict[cov_lib] = []
					cov_data_dict[cov_lib].append(cov)

				if not tax in tax_dict:
					tax_dict[tax] = {}
				tax_dict[tax]['count'] = tax_dict[tax].get('count', 0) + 1
				tax_dict[tax]['span'] = tax_dict[tax].get('span', 0) + int(length)

	contig_data_array = np.array(contig_data_list)
	cov_data_array = {}
	for cov_lib, cov_data in cov_data_dict.items():
		cov_array = np.array(cov_data, dtype=np.float)
		cov_data_array[cov_lib]=cov_array

	return contig_data_array, tax_dict, cov_data_array

def getMasks(data, len_cutoffs):
	""" Returns dict with parameters as keys and mask-array as values. """
	mask_dict = {}
	for len_cutoff in len_cutoffs:
		key = str(len_cutoff) #+ "_" + str(eval_cutoff)
		mask = np.where(data[:,1].astype(int) >= len_cutoff) #& (data[:,4].astype(float) <= eval_cutoff))	
		mask_dict[key]=mask
	return mask_dict

def getTaxGroups(tax_dict, sort_by_span, max_taxa_plot):
	"""Create OrderedDict of tax groups sorted by span or count.
	
	Each tax group maps to a list of taxonomic annotations (including special 
	annotations such as 'ambig-hit' and 'no-hit'). The tax groups are ordered 
	by the span or count of their constituent taxa. All tax groups will consist
	of one annotation, except the tax group 'other-taxa', which groups taxa of
	low span/count together so as to ensure all taxa are plotted while keeping 
	within the maximum number of taxa to plot.
	"""
	
	tax_groups = OrderedDict()
	
	k = 'span' if sort_by_span else 'count'
	
	tax_array = sorted(tax_dict, key=lambda x: tax_dict[x][k], reverse=True)

	special_groups = ('ambig-hit', 'no-hit')
	standard_indices = np.where(np.in1d(tax_array, special_groups, invert=True) )[0]
	special_indices = np.where(np.in1d(tax_array, special_groups) )[0]
	
	effective_max_taxa = max_taxa_plot - len(special_indices)

	if len(standard_indices) > effective_max_taxa:
		effective_max_groups = effective_max_taxa - 1
	else:
		effective_max_groups = effective_max_taxa
		
	for i in standard_indices[:effective_max_groups]:
		tax_groups[ tax_array[i] ] = [ tax_array[i] ]
		
	for j in special_indices:	
		tax_groups[ tax_array[j] ] = [ tax_array[j] ]
	
	if effective_max_groups < len(standard_indices):
		tax_groups['other-taxa'] = [ tax_array[i] 
			for i in standard_indices[effective_max_groups:] ]

	return tax_groups
	
def getColorDict(tax_groups, colormap="Set2"):
	""" Returns colour dict, Not annotated is always grey. """
	colors = cm.get_cmap(name=colormap)
	color_index = 1
	color_dict = {}
	
	num_colors = len(tax_groups) - sum( t in tax_groups 
		for t in ('ambig-hit', 'no-hit', 'other-taxa') )
	
	for i, g in enumerate(tax_groups):
		if g in ('ambig-hit', 'no-hit', 'other-taxa'):
			color_dict[g] = special_plot_colors[g]
		else:
			color_dict[g] = mat.colors.rgb2hex(colors(1.0 * (color_index/num_colors)))
			color_index += 1
			
	return color_dict

def getLabelDict(tax_groups, data):
	"""Returns label dict."""	
	
	label_dict = {}
	
	for g in tax_groups:
	
		m = np.in1d(data[:,3].astype(str), tax_groups[g])
		number_of_contigs_for_group = np.sum(m)
	
		# sums span for taxon group
		span_of_contigs_for_group = np.sum(data[m][:,1].astype(int))
	
		# create np_arrays for length for all contigs in taxon group
		len_array = data[m][:,1].astype(int)
	
		label_dict[g] = "{} ({:,}; {}; {})".format(g, 
			int(number_of_contigs_for_group), 
			formatBaseCount(span_of_contigs_for_group), 
			formatBaseCount(n50(len_array)))

	return label_dict

def getMinMaxCov(cov_dict):
	max_cov, min_cov = 100.00, 100.00
	for lib in cov_dict:
		lib_max_cov, lib_min_cov = np.amax(cov_dict[lib].astype(float)), np.amin(cov_dict[lib].astype(float))
		if lib_max_cov > max_cov:
			max_cov = lib_max_cov
		if lib_min_cov < min_cov:
			min_cov = lib_min_cov
	print "[STATUS] - Max.cov = " + str(max_cov) + " / Min.cov = " + str(min_cov)
	return max_cov, min_cov

black, grey, background_grey, white = '#262626', '#d3d3d3', '#F0F0F5', '#ffffff'
nullFmt = NullFormatter()  # no labels on axes
nullLoc = NullLocator()  # no ticks on axes

if __name__ == "__main__":
	colormap = "Set2" # "Paired"
	
	infile, out_prefix, max_taxa_plot, len_cutoffs, ignore_contig_len, sort_by_span, multi_plot, hist_span, title = getInput()

	data, tax_dict, cov_dict = parseInfile(infile) # data is a numpy array, tax_dict is a dict, cov_dict is a dict of numpy arrays

	mask_dict = getMasks(data, len_cutoffs) # allows filtering of blobs by length 

	# Get OrderedDict mapping of tax-groups to taxa.
	tax_groups = getTaxGroups(tax_dict, sort_by_span, max_taxa_plot)
	
	# Get mapping of tax-groups to plot colour.
	color_dict = getColorDict(tax_groups, colormap=colormap)
	
	# Get mapping of tax-groups to plot label.
	label_dict = getLabelDict(tax_groups, data)

	max_cov, min_cov = getMinMaxCov(cov_dict)

	for lib in sorted(cov_dict):
	# sanitation of names occurs in cov_dict 
		for key in mask_dict:
			mask = mask_dict[key]
			cov = cov_dict[lib]
			if out_prefix: 
				outfile = out_prefix + "." + lib + "." + key 
			else: 
				outfile =  infile + "." + lib + "." + key 
			plot(data[mask], cov[mask], outfile, title, 
				tax_groups=tax_groups, color_dict=color_dict, label_dict=label_dict, 
				ignore_contig_len=ignore_contig_len, min_cov=min_cov, max_cov=max_cov, 
				hist_span=hist_span, multi_plot=multi_plot, fig_format='png')

################################################################################
