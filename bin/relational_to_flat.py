#!/usr/local/bin/python

import argparse
from collections import defaultdict
from openpyxl.utils.cell import get_column_letter

def main(args):
	
	control_key, comparison_key = read_comparisons(args.comparisons)
	sample_groups = read_groups(args.groups)
	write_output(args.outfile, control_key, comparison_key, sample_groups)

def read_comparisons(comparison_infile):

	controls = set()
	comparison_key = defaultdict(list)
	with open(comparison_infile) as infile:
		for line in infile:
			cur = line.rstrip().split('\t')
			if len(cur) == 3:
				control, treatment, name = line.rstrip().split('\t')
				column = 'group'
			else:
				control, treatment, name, column = line.rstrip().split('\t')
			controls.add((control, column))
			comparison_key[(control, column)].append((treatment, name, column))
	
	control_key = {}
	for i, sample in enumerate(list(controls)):
		control_key[sample] = get_column_letter(i+1)

	return(control_key, comparison_key)

def read_groups(group_infile):

	sample_groups = defaultdict(list)
	with open(group_infile) as infile:
		header = infile.readline().rstrip().split('\t')
		for line in infile:
			cur = line.rstrip().split('\t')
			for i, column in enumerate(header):

				sample_groups[(cur[i], column)].append(cur[0])
	return(sample_groups)

def write_output(output_file, control_key, comparison_key, sample_groups):

	output = []

	for comparison in comparison_key:
		letter = control_key[comparison]

		for sample in sample_groups[comparison]:
			output.append('%s\t%s0\t' % (sample, letter))

		for i, group in enumerate(comparison_key[comparison]):
			for sample in sample_groups[(group[0], group[2])]:
				output.append('%s\t%s%s\t%s' % (sample, letter, i+1, group[1]))

	with open(output_file, 'w') as out:
		out.write('\n'.join(output))

def parseArguments():
	parser = argparse.ArgumentParser(prog="Convert from groups and comparisons files into flat file.", description='', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-g', '--groups', required=True, help='Groups file.', metavar='', dest='groups')
	required.add_argument('-c', '--comparisons', required=True, help='Comparisons file.', metavar='', dest='comparisons')
	optional = parser.add_argument_group('optional arguments')
	optional.add_argument('-o', '--output', default='groups.tsv', help='Output file.', metavar='', dest='outfile')

	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)