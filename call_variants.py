import subprocess
import os
import sys
import glob
import signal
from Bio import SeqIO

class VariantCaller:
	"""
	Call variants using either samtools/bcftools or freebayes
	"""

	def __init__(self, final_results_directory, reference, threads, minq):
		self.final_results_directory = final_results_directory
		self.reference = reference
		self.threads = threads
		self.minq = minq
	
	def gather_bamfiles(self, final_results_directory):
		bam_names = []
		for f in glob.glob(final_results_directory+"*_snpbac.bam"):
			bam_names.append(f)
		return(bam_names) 

	def samtools_variant_caller(self, final_results_directory, reference, threads, bam_names):
		additional_threads = str(int(threads) - 1)
		with open(final_results_directory+"samtools_error_log.txt", "w") as fb_error:
			with open(final_results_directory+"snpbac_raw.vcf", "w") as raw_vcf:
				samtools_mpileup = ["samtools", "mpileup", "-uf", reference]
				samtools_mpileup.extend(bam_names)
				samtools_mpileup_call = subprocess.Popen(samtools_mpileup, stdout = subprocess.PIPE)
				bcftools_call = ["bcftools", "call", "-mv", "--ploidy", "1", "--threads", additional_threads]
				bcftools_call_call = subprocess.Popen(bcftools_call, stdin = samtools_mpileup_call.stdout, stdout = raw_vcf, stderr = fb_error)
				samtools_mpileup_call.stdout.close()
				bcftools_call_call.communicate()

	def freebayes_variant_caller(self, final_results_directory, reference, threads, bam_names):
		with open(final_results_directory+"freebayes_error_log.txt", "w") as fb_error:
			with open(final_results_directory+"snpbac_raw.vcf", "w") as raw_vcf:
				freebayes = ["freebayes", "-f", reference, "-p", "1"]	
				bam_command = []
				for b in bam_names:
					bam_command.append("-b")
					bam_command.append(b)	
				freebayes.extend(bam_command)
				freebayes_call = subprocess.Popen(freebayes, stdout = raw_vcf, stderr = fb_error).wait()	

	def filter_variants(self, final_results_directory, minq):
		minq = str(minq)
		snps = ["vcftools", "--vcf", final_results_directory+"snpbac_raw.vcf", "--minQ", minq, "--out", final_results_directory+"snpbac_snps", "--recode", "--remove-indels"]
		snps_call = subprocess.Popen(snps).wait()
		indels = ["vcftools", "--vcf", final_results_directory+"snpbac_raw.vcf", "--minQ", minq, "--out", final_results_directory+"snpbac_indels", "--recode", "--keep-only-indels"]
		indels_call = subprocess.Popen(indels).wait()

	def subset_vcf(self, final_results_directory):
		sample_list = []
		vcf = open(final_results_directory+"snpbac_snps.recode.vcf", "r")
		for line in vcf:
			if "##" not in line and "#" in line:
				splits = line.split("\t")[9:]
				print splits
				sample_list.extend([s.strip() for s in splits])
		print "sample_list"
		print sample_list
		vcf.close()
		bgzip_vcf = ["bgzip", final_results_directory+"snpbac_snps.recode.vcf"]
		bgzip_vcf_call = subprocess.Popen(bgzip_vcf).wait()
		tabix = ["tabix", "-p", "vcf", final_results_directory+"snpbac_snps.recode.vcf.gz"]
		tabix_vcf_call = subprocess.Popen(tabix).wait()
		for sample in sample_list:
			print sample
			with open(final_results_directory+sample+".vcf.gz", "w") as sample_vcf:
				vcf_subset = ["vcf-subset", "-c", sample, "-e", final_results_directory+"snpbac_snps.recode.vcf.gz"]
				vcf_subset_call = subprocess.Popen(vcf_subset, stdout = subprocess.PIPE)
				bgzip_sample = ["bgzip"]
				bgzip_sample_call = subprocess.Popen(bgzip_sample, stdin = vcf_subset_call.stdout, stdout = sample_vcf)
				vcf_subset_call.stdout.close()
				bgzip_sample_call.communicate()
				tabix_sample = ["tabix", "-p", "vcf", final_results_directory+sample+".vcf.gz"]
				print "Tabixing..."
				tabix_sample_call = subprocess.Popen(tabix_sample).wait()
		gunzip_main_vcf = ["gunzip", final_results_directory+"snpbac_snps.recode.vcf.gz"]
		gunzip_main_vcf_call = subprocess.Popen(gunzip_main_vcf).wait()

	def default_sigpipe(self):
		signal.signal(signal.SIGPIPE, signal.SIG_DFL)

	def vcf_consensus(self, final_results_directory, reference):
		for f in glob.glob(final_results_directory+"*.vcf.gz"):
			print "Making consensus sequence for "+f.split(".vcf.gz")[0]
			with open(f.split(".vcf.gz")[0]+"_temporary.fasta", "w") as sample_consensus:
				cat = ["cat", reference]
				cat_call = subprocess.Popen(cat, stdout = subprocess.PIPE, preexec_fn=self.default_sigpipe)
				consensus = ["vcf-consensus", f]
				consensus_call = subprocess.Popen(consensus, stdin = cat_call.stdout, stdout = sample_consensus, preexec_fn=self.default_sigpipe)
				cat_call.stdout.close()
				consensus_call.communicate()

	def rename_concat_consensus(self, final_results_directory):
		for f in glob.glob(final_results_directory+"*_temporary.fasta"):
			outfile = open(f.split("_temporary.fasta")[0]+"_consensus.fasta", "a")
			seq_concat = ""
			infile = open(f, "r")
			for record in SeqIO.parse(infile, "fasta"):
				seqseq = str(record.seq).strip()
				seq_concat += seqseq
			infile.close()
			seqid = f.split("/")[-1]
			seqid = seqid.split("_temporary.fasta")[0]
			print >> outfile, ">"+seqid.strip()
			print >> outfile, seq_concat
			outfile.close()
			os.remove(f)

