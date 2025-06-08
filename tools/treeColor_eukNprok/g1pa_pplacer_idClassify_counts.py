#!/usr/bin/env python

# count_pplacer_csv_by_taxonomy.py

"""
Ryan Groussman
Armbrust Lab, University of Washington 2017
This script will take as input the 'guppy to_csv' CSV output from pplacer jplace
files. It will go through the list of pqueries (designed for the 'keep
at most 1' option) and will count the incidence in each 'classification' field
against a list of NCBI tax_ids for groups of particular taxonomy (e.g. Bacillariophyta)
supplied by the --treecolor_csv flag.
The normalization factors are given by the --normfactors_csv flag.
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_csv", help="CSV-formatted input file")
parser.add_argument("-g", "--high_level_groups", help="Print high-level group results", action="store_true")
parser.add_argument("-t", "--top_taxa", help="Print top taxa per group", action="store_true")
parser.add_argument("-i", "--ingroup_edges", help="Only count edge numbers in this file", type=str)
parser.add_argument("-c", "--treecolor_csv", help="Specify path to an alternate treecolor.csv-like file.", type=str)
parser.add_argument("-e", "--print_header", help="Print header for csv and quit", action="store_true")
parser.add_argument("-n", "--normfactors_csv", help="Specify path to an csv file with norm factors for your samples.", type=str)
parser.add_argument("-w", "--write_unknown", help="Write out all the tax_ids in the 'unknown' group and quit", action="store_true")
parser.add_argument("-x", "--write_edge_counts", help="Write the normalized counts for each edge to specified file", type=str)

# Parse the arguments given at the command line
args = parser.parse_args()
other_counter = 0
outgroup_counter = 0
other_tax_set = set([])

def parse_treecolor_info(TCInfo):
	"""
	This function loads a 'treecolor.csv' type file, which tells us the names of the
	taxonomic groups to count and the corresponding files with tax_ids.
	Returns a dictionary with this information.
	"""

	TreecolorDict = {} # treecolor information dictionary

	for line in TCInfo:
		line_elts = line.split(",")
		ListFile = line_elts[0]	# the group name (e.g., Bacillariophyta)
		ListFilePath = '/'.join(args.treecolor_csv.split('/')[0:-1]) + '/' + ListFile.strip()
		group_name = ListFile.split(".")[0]	# pop off the prefix, this is the clade name
		tax_file = open(ListFilePath, 'r') # open the list file
		tax_set = set([]) # make a set to store the tax_ids
		for tax_id in tax_file:
			tax_id = tax_id.strip()
			tax_set.add(tax_id)
		TreecolorDict[group_name] = tax_set

	return TreecolorDict

def test_treecolor_dict(TreecolorDict):
	"""
	Tests whether or not there are overlapping tax_id between the different
	groups. Prints a warning to stdout if yes; no message if not."""

	# Test that there are no overlapping tax_ids between the different groups in TreecolorDict:
	tax_group_list = TreecolorDict.keys()
	while len(tax_group_list) > 0:
		any_group = tax_group_list.pop()
		for other_group in tax_group_list:
			overlap = TreecolorDict[any_group].intersection(TreecolorDict[other_group])
			if len(overlap) > 0:
				print "!!! WARNING !!!"
				print "Overlap between ", any_group, "and", other_group
				print overlap
				quit()

def initialize_norm_counts_dict():
	"""
	This function sets up a dictionary with sample normalization factors.
	The	default is diel1_standard_counts_SUMS.txt (via Bryn Durham, 13Jan17).
	Norm factor is the sum of BPDSmall1, BPDSmall2, BPDMedium1, BPDMedium2,
	BPDLarge1 and BPDLarge2 from bowtie2 counts (first summed by for each standard
	across machine runs and then summed in total), then averaged by # of standards (divided by 6),
	and used to divide total 'added' counts along with volume:
	NORM FACTOR = ADDED / AVG / VOLUME
	If called with -n or --normfactors_csv, will load a CSV with norm factors
	in this format specified to your samples, with sample names and machine
	reps in the first field separated by a period (no header). Without this call,
	this script will load the default diel1 norm factors.
		S14C1_0_2umC.NS7,4779.004092
		S14C1_0_2umC.NS8,4891.099897
		S14C1_3umA.NS1,51153.58824
		S14C1_3umA.NS2,17896.11284
	"""

	if args.normfactors_csv != None:
		NormFactorsDict = {}
		normfactors_csv = open(args.normfactors_csv, 'r')
		for line in normfactors_csv:
			line_elts = line.split(",")
			sample = line_elts[0]
			norm_factor = float(line_elts[1].strip())
			NormFactorsDict[fix_sample_name(sample)] = norm_factor
		return NormFactorsDict

	elif args.normfactors_csv == None:
		# ( to maintain back compatibility - these are the 'default' norm factors if not called with --normfactors_csv )
		# Values taken from diel1_standard_counts_SUMS.txt (via Bryn Durham, 13Jan2017)
		NormFactorsDict = {
		'S11C1_A_1800':2393,
		'S11C1_C_1800':2254,
		'S14C1_B_2200':2978,
		'S14C1_C_2200':2365,
		'S15C1_B_200':2674,
		'S15C1_C_200':1970,
		'S16C1_A_600':2597,
		'S16C1_B_600':2520,
		'S17C1_A_1000':2270,
		'S17C1_B_1000':1904,
		'S18C1_A_1400':1339,
		'S18C1_C_1400':2287,
		'S19C1_A_1800':2485,
		'S19C1_C_1800':2097,
		'S20C1_B_2200':2506,
		'S20C1_C_2200':2666,
		'S21C1_B_200':2107,
		'S21C1_C_200':2172,
		'S22C1_A_600':2059,
		'S22C1_C_600':1654,
		'S23C1_B_1000':3005,
		'S23C1_C_1000':1577,
		'S24C1_A_1400':1917,
		'S24C1_B_1400':1466,
		'S26C1_A_1800':2833,
		'S26C1_C_1800':2327,
		'S28C1_B_2200':1247,
		'S28C1_C_2200':2669,
		'S29C1_A_200':2240,
		'S29C1_C_200':1348,
		'S30C1_A_600':1831,
		'S30C1_C_600':2598,
		'S31C1_A_1000':2192,
		'S31C1_C_1000':1834,
		'S32C1_B_1400':1791,
		'S32C1_C_1400':2055,
		'S33C1_A_1800':1949,
		'S33C1_C_1800':2500,
		'S34C1_B_2200':1858,
		'S34C1_C_2200':2641,
		'S35C1_A_200':2625,
		'S35C1_C_200':2630,
		'S06C1_A_600':1396,
		'S06C1_C_600':3435,
		'S07C1_A_1000':1777,
		'S07C1_B_1000':1560,
		'S08C1_B_1400':1887,
		'S08C1_C_1400':4185
		}

	return NormFactorsDict

def normalize_counts(sample, raw_sum_count):
	"""
	Normalizes against the combined machine-run sums of counts from
	diel1_standard_counts_SUMS.txt (via Bryn Durham, 13Jan17)
	Norm factor is the sum of BPDSmall1, BPDSmall2, BPDMedium1, BPDMedium2,
	BPDLarge1 and BPDLarge2 from bowtie2 counts, then averaged, and used to
	divide total 'added' counts along with volume:
	NORM FACTOR = ADDED / AVG / VOLUME
	Takes as input the sample_id and raw_count, returns the normalized count.
	"""

	if sample in NormFactorsDict.keys():
		norm_count = raw_sum_count * NormFactorsDict[sample]
	else:
		print "Sample ID not found!"
	return norm_count

def initialize_taxa_counts_dict(TCInfo):
	"""Build a blank dictionary like default_taxa_counts_dict
	but using specified into in alternate treecolor.csv-like file."""

	TaxaCountsDict = {}
	for line in open(TCInfoPath, 'r'):
		base_file = line.split(",")[0]
		group_name = base_file.split(".")[0]
		TaxaCountsDict[group_name] = 0
	return TaxaCountsDict

def default_taxa_counts_dict():
	"""Initialize TaxaCountsDict starting with 0 for desired taxons.
	"""

	# Initialize the taxa counts dictionary
	TaxaCountsDict = {
	'Amoebozoa' : 0,
	'Bacillariophyta' : 0,
	'Cryptophyta' : 0,
	'Dinophyceae' : 0,
	'Excavata' : 0,
	'Glaucophyta' : 0,
	'Haptophyta' : 0,
	'NonDiatomStramenopile' : 0,
	'NonDinoAlveolata' : 0,
	'Opisthokonta' : 0,
	'Rhizaria' : 0,
	'Rhodophyta' : 0,
	'Viridiplantae' : 0
	}

	return TaxaCountsDict

def print_high_level_results():
	"""Prints the high-level taxonomic results to stdout"""

	# Spit out results!
	counts_results = [fix_sample_name(summed_sample_name)]
	# Normalize for all the groups in TaxaCountsDict and append to counts_results list
	for group_name in sorted(TaxaCountsDict.keys()):
		norm_count = normalize_counts(summed_sample_name, TaxaCountsDict[group_name])
		counts_results.append(str(norm_count))
	# Normalize counts for 'other' and 'outgroup' counters
	norm_count_other = normalize_counts(summed_sample_name, other_counter)
	counts_results.append(str(norm_count_other))
	norm_count_outgroup = normalize_counts(summed_sample_name, outgroup_counter)
	counts_results.append(str(norm_count_outgroup))
	# line below prints the header
	if args.print_header == True:
		print 'sample_name' + "," + ",".join(sorted(TaxaCountsDict.keys())) + ",Other,Outgroups"
	elif args.write_unknown == True:
		for other in other_tax_set:
			print other + "\n"
	elif args.print_header == False:
		print ",".join(counts_results)

def print_top_taxa_per_group(group):
	"""Given a group, will run through all of the tax_id that run under its
	umbrella. Prints out a list of tax_ids belonging to group along with
	the proportion of hits in the file (a quick way to normalize)"""

	# header
	# top_taxa_header = "sample_name,group,tax_id,fraction_of_all_hits"

	# for our taxa in TreecolorDict
	for taxa in TreecolorDict[group]:
		if taxa in TaxIDDict.keys():
			hit_fraction = TaxIDDict[taxa] / count_denominator
			print fix_sample_name(sample_name) + "," + group + "," + str(taxa) + "," + str(hit_fraction)

def add_classification_to_taxiddict(classification):
	"""Given an NCBI tax_id classification, will add it to a dictionary if
	it's not there with value 1, otherwise will increment counter"""

	if classification in TaxIDDict.keys():
		TaxIDDict[classification] += 1
	elif classification not in TaxIDDict.keys():
		TaxIDDict[classification] = 1

def process_ingroup_list(ingroup_file_path):
	""" Process a list of edge numbers for the ingroup and their inclusive
	internal edge numbers, as noted in a .jplace file. If given as
	argument, this script will only count these edges. Example:
	72
	73
	74
	...
	"""

	ingroup_edges_file = open(args.ingroup_edges, 'r')
	ingroup_edges_set = set([])
	for edge in ingroup_edges_file:
		edge = edge.strip()
		ingroup_edges_set.add(edge)
	return ingroup_edges_set

def search_and_count_treecolor_dict(tax_id):
	"""go through each tax list in TreecolorDict
	if it matches a list, increment tally in TaxaCountsDict
	"""
	global other_counter
	global other_tax_set
	global outgroup_counter

	# taxID starts off not found:
	taxid_found = False
	# If it's found in one of the sets of the TreecolorDict, then increment a counter:
	# test_treecolor_dict(TreecolorDict) prevents overlap between the TreecolorDict sets
	for group_name in TreecolorDict:
		if tax_id in TreecolorDict[group_name]:
			TaxaCountsDict[group_name] += 1
			taxid_found = True
	# Otherwise, add this taxID to a list of others and increment a counter
	if taxid_found == False:
		other_counter += 1
		if args.write_unknown == True:
			other_tax_set.add(tax_id)

def fix_sample_name(sample_name):
	"""Fixes single-digit station sample names to standardize format.
	E.g. S8C1 > S08C1, S7C1 > S07C1, etc """

	for i in range(1,10):
		short_id = "S" + str(i) + "C1"
		long_id = "S0" + str(i) + "C1"
		if sample_name.startswith(short_id):
			sample_name = sample_name.replace(short_id, long_id)
	return sample_name

def main():
	"""Execute the main sequence of commands for counting pplacer csv files"""

	global TCInfoPath
	global TaxIDDict
	global TreecolorDict
	global TaxaCountsDict
	global summed_sample_name
	global NormFactorsDict
	global outgroup_counter

	# load csv file containing list and color information
	if args.treecolor_csv != None:
		TCInfoPath = args.treecolor_csv
	elif args.treecolor_csv == None:
		TCInfoPath = "/blue/b.durham/rebeccakey/4_boysenTrees/tools/euk_tools/treecolor_files_euk/treecolors_w_euks.csv"

	# parse treecolor info: custom high-level classification
	TCInfo = open(TCInfoPath, 'r')
	TreecolorDict = parse_treecolor_info(TCInfo)
	# test the dictionary for any overlapping sets of tax_ids:
	test_treecolor_dict(TreecolorDict)

	# if ingroup list file is provided, load the file:
	if args.ingroup_edges != None:
		ingroup_edges_set = process_ingroup_list(args.ingroup_edges)

	# Initialize TaxaCountsDict and Other count and Outgroup count
	if args.treecolor_csv != None:
		TaxaCountsDict = initialize_taxa_counts_dict(TCInfo)
	elif args.treecolor_csv == None:
		TaxaCountsDict = default_taxa_counts_dict()

	# Initialize NormFactorsDict
	NormFactorsDict = initialize_norm_counts_dict()

	# now open up the guppy csv file
	input_csv_path = args.input_csv
	in_csv = open(input_csv_path, 'r')

	# Counts per edge: this dict stores raw counts for each taxid.
	TaxIDDict = {}

	# Start off summed_sample_name as False until it's been given an assignment
	summed_sample_name = False
	line_counter = 0
	in_csv.readline()   # skip the first line (header)
	# iterate through each subsequent line of the CSV file:
	for line in in_csv:
		line_counter += 1 # increment line counter
		line_elts = line.split(",") # comma-separated values
		classification = line_elts[10].strip() # classification is in '10th' field
		edge_num = line_elts[3] # this is used to determine whether placement is in an 'ingroup'
###		
		origin = line_elts[0] # sample ID, may include machine runs
		line_sample_name = origin.split(".")[1]
#		line_sample_name = origin.split(".NS")[0]
#		line_sample_name = sample_name.split("atrr.")[1]
###		
		# 5 following lines assume that the CSV is from one discrete sample (although multipe machine runs may be present)
		if summed_sample_name == False:
			summed_sample_name = line_sample_name
		elif summed_sample_name != line_sample_name:
			print "Multiple samples present in CSV file!"
			quit()
		# If ingroup flag is set, will classify counts outside of this group as 'outgroup'
		if args.ingroup_edges != None:
			if edge_num in ingroup_edges_set:
				search_and_count_treecolor_dict(classification)
				add_classification_to_taxiddict(classification)
			elif edge_num not in ingroup_edges_set:
				outgroup_counter += 1
		# With no ingroups (default), classify through following functions
		elif args.ingroup_edges == None:
			add_classification_to_taxiddict(classification)
			search_and_count_treecolor_dict(classification)

	# This flag outputs results for the high level groups given by the supplied treecolor file:
	if args.high_level_groups == True:
		print_high_level_results()

	# top_taxa will print out the most significant taxa per high level 'treecolor' group
	elif args.top_taxa == True:
		if args.ingroup_edges == None:
			count_denominator = float(line_counter)
		elif args.ingroup_edges != None:
			count_denominator = float(line_counter - outgroup_counter)
		for group in TreecolorDict.keys():
			print_top_taxa_per_group(group)
		# also 'other' and 'outgroup'
		# removed division by count_denominator for other_counter and outgroup_counter (ex: (other_counter / count_denominator))
		print fix_sample_name(sample_name) + "," + "Other" + "," + "1" + "," + str(other_counter / count_denominator)
		if args.ingroup_edges != None:
			print fix_sample_name(sample_name) + "," + "Outgroups" + "," + "1" + "," + str(outgroup_counter / count_denominator)

	in_csv.close()

if __name__ == "__main__":
	main()
