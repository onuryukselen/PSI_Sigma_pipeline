#!/usr/bin/python3

import sys, re, argparse
from collections import defaultdict

def main(args):

	sample_key, columns, names, comparisons = parse_input(args.infile)
	write_groups_file(sample_key, columns, args.groups)
	write_comparisons_file(comparisons, names, args.comparisons)

def parse_input(input_file):

	sample_key = {}
	columns = set()
	names = {}
	comparisons = defaultdict(set)

	with open(input_file) as infile:
		for line in infile:
			sample, code, name = line.rstrip().split('\t')

			if sample not in sample_key:
				sample_key[sample] = defaultdict(list)

			letter = re.findall('^[A-Z]+', code)[0]
			number = re.findall('\d+$', code)[0]
			
			columns.add(letter)
			sample_key[sample][letter] = number
			names[code] = name

			if number != '0':
				comparisons[letter].add(number)

	return sample_key, columns, names, comparisons

def write_groups_file(sample_key, columns, outfile):
	final_columns = sorted(list(columns))
	output = ['sample_name\t%s' % ('\t'.join(column for column in final_columns))]
	for sample in sample_key:
		output.append('%s\t%s' % (sample, '\t'.join([sample_key[sample][column] if len(sample_key[sample][column]) > 0 else '' for column in final_columns])))
	with open(outfile, 'w') as out:
		out.write('\n'.join(output))

def write_comparisons_file(comparisons, names, outfile):
	output = ['controls\ttreats\tnames\tcolumn']
	for comparison_group in comparisons:
		for i in comparisons[comparison_group]:
			output.append('%s\t%s\t%s\t%s' % ('0', i, names['%s%s' % (comparison_group, i)], comparison_group))
	with open(outfile, 'w') as out:
		out.write('\n'.join(output))

def parseArguments():
	parser = argparse.ArgumentParser(prog="Convert from flat file into groups and comparisons files.", description='', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', '--input', required=True, help='Input file [sample] [Code] [Name]', metavar='', dest='infile')
	optional = parser.add_argument_group('optional arguments')
	optional.add_argument('-g', '--groups-outfile',  default='groups.tsv', help='Output name for groups file.', metavar='', dest='groups')
	optional.add_argument('-c', '--comparisons-outfile', default='comparisons.tsv', help='Output name for comparisons file.', metavar='', dest='comparisons')

	return parser.parse_args()

args = parseArguments()

main(args)