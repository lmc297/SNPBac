from collections import OrderedDict

class CoresnpMatrix:
	"""
	Create a matrix of core SNPs
	"""
	
	def __init__(self, final_results_directory, recomb_filter):
		self.final_results_diretory = final_results_directory
		self.recomb_filter = recomb_filter
	
	def create_core_snp_matrix(self, final_results_directory, recomb_filter):
		bad_snps = [] # contains non-core SNPs and SNPs in recombined regions, if applicable
		if recomb_filter == "gubbins":
			print "Collecting recombination positions from gubbins..."
			gubbins = open(final_results_directory+"gubbins/snpbac_gubbins.recombination_predictions.embl","r")
			for line in gubbins:
				if "misc_feature" in line:
					snp_range=line.split()[2]
					snp_begin = int(snp_range.split("..")[0].strip())
					snp_end = int(snp_range.split("..")[1].strip())
					realrange = range(snp_begin, snp_end + 1)
					bad_snps = bad_snps + [int(r) for r in realrange]
			gubbins.close()
	
		print "Finished with gubbins."
		in_vcf = open(final_results_directory+"snpbac_snps.recode.vcf", "r")
		out_vcf = open(final_results_directory+"snpbac_core_snps.vcf", "a")
		snpmat = {}
		for line in in_vcf:
			if "#" not in line:
				position = int(line.split("\t")[1].strip())	
				alleles = []
				ref = line.split("\t")[3]
				alt = line.split("\t")[4]
				alt = alt.split(",")
				alleles.append(ref.strip())
				alleles.extend(alt)
				if position not in bad_snps:
					if position not in snpmat.keys():
						snpmat[position] = []
					splits = line.split("\t")[9:]
					snps = [s.split(":")[0].strip() for s in splits]
					if "." not in snps and len(set(snps)) > 1:
						newsnps = [alleles[int(s)] for s in snps]
						snpmat[position].append([s for s in newsnps])
						print >> out_vcf, line.strip()
			else:
				print >> out_vcf, line.strip()
				if "##" not in line:
					samples = line.split("\t")[9:]
		in_vcf.close()
		out_vcf.close()

		out_snpmatrix = open(final_results_directory+"snpbac_core_snps.txt", "a")
		print >> out_snpmatrix, "position" + "\t" + "\t".join([s.strip() for s in samples])
		print "Printing core SNP matrix..."
		for key in sorted(snpmat.keys()):
			mysnps = snpmat[key]
			for m in mysnps:
				print >> out_snpmatrix, str(key) + "\t" + "\t".join([mm.strip() for mm in m])
		out_snpmatrix.close()
					

	def create_snp_fasta(self, final_results_directory):

		snp_dict = OrderedDict()
		
		infile = open(final_results_directory+"snpbac_core_snps.txt", "r")
		for line in infile:
			splits = line.split("\t")
			if "position" in line:
				for s in splits:
					if s.strip() != "position":
						snp_dict[s.strip()] = ""
			else:
				counter = 1
				for key in snp_dict.keys():
					snp_dict[key] = snp_dict[key] + splits[counter].strip()
					counter = counter + 1
		infile.close()
		
		outfile = open(final_results_directory+"snpbac_core_snps.fasta", "a")
		for key in snp_dict.keys():
			sequence = snp_dict[key]
			print >> outfile, ">" + key
			print >> outfile, sequence.strip()
		outfile.close()	
