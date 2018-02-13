#!/usr/bin/env python3

# usage: python make_snpbac_infile.py </path/to/directory/with/reads/> <name_of_outfile.txt> <"_forward_reads_suffix.fastq.gz"> <"_reverse_reads_suffix.fastq.gz">

# example: python make_snpbac_infile.py /path/to/my/reads/ /path/to/my/reads/read_list.txt "_R1.fastq.gz" "_R2.fastq.gz"

# another example: python make_snpbac_infile.py /a/different/path/to/my/reads/ /directory/where/I/want/my/list_of_reads_named_anything_I_want.txt "_1.fastq.gz" "_2.fastq.gz"

import sys, os, glob
import argparse

def main():

	parser = argparse.ArgumentParser(usage = 'python make_snpbac_infile.py --input </path/to/directory/with/reads/> --out <name_of_outfile.txt> --forward <"_forward_reads_suffix.fastq.gz"> --reverse <"_reverse_reads_suffix.fastq.gz">')

	parser.add_argument("--input", help = "Path to directory with all reads in fastq.gz format", nargs = 1, required = True)

	parser.add_argument("--out", help = "Desired name of your list of sample file paths; if you want the list deposited in a directory other than the current one, make sure to include a path", nargs = 1, required = True)

	parser.add_argument("--forward", help = 'Suffix shared by all of your forward reads, in quotation marks; common suffixes might include "_1.fastq.gz" or "_R1.fastq.gz"', nargs = 1, required = True)
	parser.add_argument("--reverse", help = 'Suffix shared by all of your reverse reads, in quotation marks; common suffixes might include "_2.fastq.gz" or "_R2.fastq.gz"', nargs = 1, required = True)

	args = parser.parse_args()
	
	reads = args.input[0]
	
	outfile = open(args.out[0],"a")
	
	forward = args.forward[0]
	
	reverse = args.reverse[0]

	if not reads.endswith("/"):
		reads = reads + "/"

	for f in glob.glob(reads+"*"+forward):
		prefix = f.split(forward)[0]
		if os.path.exists(prefix.strip()+reverse):
			print >> outfile, os.path.abspath(f)+","+os.path.abspath(prefix.strip()+reverse)
	outfile.close()

if __name__ == "__main__":
	main()
