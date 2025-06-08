#!/usr/bin/env python

import re
import xml.etree.ElementTree as ET
from Bio import Phylo
import argparse

def build_output_handle(infile_path):
	handle_elts = infile_path.split(".")
	handle_elts.insert(-1,"col")
	out_xml_path = ".".join(handle_elts)
	return out_xml_path

def build_color_dict(TCInfoPath):
	"""
	Builds a dictionary of the color scheme from treecolors.csv
	"""

	treecolor_file = "treeColor_key.csv"
	if args.treecolor_csv != None:
		treecolor_file = args.treecolor_csv

	TCInfo = open((TCInfoPath + treecolor_file), 'r')

	TreecolorDict = {} # treecolor information dictionary
	for line in TCInfo:
		line_elts = line.split(",")	# tab delimited file
		ListFile = line_elts[0]	# the group name
		ListFilePath = TCInfoPath + ListFile.strip()
		group_name = ListFile.split(".")[0]	# pop off the prefix, this is the clade name
		tax_file = open(ListFilePath, 'r') # open the list file
		tax_set = set([])
		for tax_id in tax_file: # create a list of tax_id from each list file
			tax_id = tax_id.strip()
			tax_set.add(tax_id.strip())
		RedVal = line_elts[2] # Collect RGB color values.
		GreenVal = line_elts[3]
		BlueVal = line_elts[4]
		rgb_values = (RedVal.strip(), GreenVal.strip(), BlueVal.strip())
		TreecolorDict[group_name] = (rgb_values,tax_set) # Color values and MMETSP# assigned as values to group names

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
			overlap = TreecolorDict[any_group][1].intersection(TreecolorDict[other_group][1])
			if len(overlap) > 0:
				print "Overlap between ", any_group, "and", other_group
				print overlap

def retrieve_tax_id(defline):
	"""
	Given a sequence defline; will retrieve the NCBI taxonomy ID given
	somewhere in the defline as '_taxXXXXXX'
	e.g.: Stauroneis_constricta_0119572860_MMETSP1352_tax265584/330-455
	Returns the NCBI taxi_id (e.g. 265584)
	"""

	# Check if there is a tax_id marker to be found:
	if re.search("(_tax[0-9]+)", defline):
		# if it's found, retrieve the tax_id number and return it
		tax_id = re.sub(r'.*_tax([0-9]+).*', r'\1', defline)
		return tax_id
	else:
		print "No tax_id found for", defline
		return None

def get_colors_from_color_dict(tax_id):
	"""
	Given an NCBI tax_id, looks up the tax_id in TreecolorDict and
	returns a tuple with the R, G, B values if present. Otherwise, returns none.
	"""

	for tax_group in TreecolorDict:
		if tax_id in TreecolorDict[tax_group][1]:
			return TreecolorDict[tax_group][0]

def change_taxa_colors(parent, rgb_values):
	"""
	Changes the colors in the XML file for a given tax_id in accordance with
	the user parameters in the color dictionary.
	"""

	if parent.find('{http://www.phyloxml.org}color') != None:
		color = parent.find('{http://www.phyloxml.org}color')
		color.find('{http://www.phyloxml.org}red').text = rgb_values[0]
		color.find('{http://www.phyloxml.org}green').text = rgb_values[1]
		color.find('{http://www.phyloxml.org}blue').text = rgb_values[2]
	elif parent.find('{http://www.phyloxml.org}color') == None and args.all_leaves == True:
		ET.SubElement(parent, '{http://www.phyloxml.org}color')
		color = parent.find('{http://www.phyloxml.org}color')
		ET.SubElement(color, '{http://www.phyloxml.org}red')
		ET.SubElement(color, '{http://www.phyloxml.org}green')
		ET.SubElement(color, '{http://www.phyloxml.org}blue')
		color.find('{http://www.phyloxml.org}red').text = rgb_values[0]
		color.find('{http://www.phyloxml.org}green').text = rgb_values[1]
		color.find('{http://www.phyloxml.org}blue').text = rgb_values[2]

parser = argparse.ArgumentParser()
parser.add_argument("input_xml", help="XML-formatted phylogenetic tree")
parser.add_argument("-a", "--all_leaves", help="Color all leaves", action="store_true")
parser.add_argument("-f", "--fat_tree", help="Color leaves with any given width", action="store_true")
parser.add_argument("-m", "--min_width", help="Used with --fat_tree; specify the minimum width to color the leaf", type=float)
parser.add_argument("-c", "--treecolor_csv", help="Specify path to an alternate treecolor.csv-like file.", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile", type=str)

args = parser.parse_args()

# Default colors for internal nodes:
default_rgb = ('150','150','150')

# load tab delimited file containing list and color information
TCInfoPath = "/blue/b.durham/rebeccakey/4_boysenTrees/tools/treeColor_eukNprok/treeColor_ids/"

input_xml_path = args.input_xml
TreecolorDict = build_color_dict(TCInfoPath)
test_treecolor_dict(TreecolorDict)

if args.out_file != None:
	out_xml_path = args.out_file
else:
	out_xml_path = build_output_handle(input_xml_path)


### test code ### elementtree
in_xml = ET.parse(input_xml_path)
root = in_xml.getroot()


for parent in root.getiterator():
	# print parent.tag
	# if the element is a clade and does not have a name (not a leaf):
	if parent.tag == '{http://www.phyloxml.org}clade' and parent.find('{http://www.phyloxml.org}name') == None:
		change_taxa_colors(parent, default_rgb)
	# if the element is a 'clade' and contains a 'name' then:
	if parent.tag == '{http://www.phyloxml.org}clade' and parent.find('{http://www.phyloxml.org}name') != None:
		# for child in parent:
		# 	print child.tag
		name_tag = parent.find('{http://www.phyloxml.org}name').text
		# print name_tag
		# print parent.text
		tax_id = retrieve_tax_id(parent.find('{http://www.phyloxml.org}name').text)
		if tax_id != None:
			rgb_values = get_colors_from_color_dict(tax_id)
			if rgb_values != None:
				change_taxa_colors(parent, rgb_values)
			elif rgb_values == None:
				print "### tax_id not in a group:", tax_id

# Iteratively run through the parent nodes; if all of their children
# share the same color scheme, then give the parent the same color.

# Write back to a file
in_xml.write(out_xml_path, xml_declaration=True)
#################

# Need to 'reset' the xml for strange formatting reasons:
tree = Phylo.read(out_xml_path, 'phyloxml')
Phylo.write(tree, out_xml_path, "phyloxml")
