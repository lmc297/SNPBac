import subprocess
import shlex

class ReadMapper:
	"""
	Map reads to a reference using bwa mem or bowtie2, and process using samtools
	"""

	def __init__(self, prefix, seq, reference, threads, final_results_directory):
		self.prefix = prefix
		self.seq = seq
		self.reference = reference
		self.threads = threads
		self.final_results_directory = final_results_directory
		#self.params = params
	
	def concatenate_reference_contigs(self, reference, final_results_directory):
		ref_name = reference.split(".")[0]
		if "/" in ref_name:
			ref_name = ref_name.split("/")[-1]
		ref_in = open(reference, "r")
		ref_out = open(final_results_directory+ref_name.strip()+"_pseudochrom.fasta", "a")
		print >> ref_out, ">"+ref_name.strip()
		print >> ref_out, "".join(line.strip() for line in ref_in if ">" not in line)
		ref_out.close()
		ref_in.close()
		reference = final_results_directory+ref_name.strip()+"_pseudochrom.fasta"
		return(reference)
		

	def  bwa_setup(self, reference, threads, seq, prefix):
		bwa_index = ["bwa", "index", reference]# shlex.split(params)
		index = subprocess.Popen(bwa_index).wait()
		read_group = "@RG\tID:"+prefix.split("/")[-1].strip()+"\tSM:"+prefix.split("/")[-1].strip()
		bwa_mem = ["bwa", "mem", "-t", threads, "-R", read_group.strip(), reference]
		bwa_mem.extend(seq)
		return(bwa_mem)

	def bowtie2_setup(self, reference, threads, seq, prefix):
		bowtie_index = ["bowtie2-build", "-f", reference, reference]
		index = subprocess.Popen(bowtie_index).wait()
		read_group = prefix.split("/")[-1].strip()
		if len(seq) > 1:
			bowtie2 = ["bowtie2", "-x", reference, "--rg-id", read_group, "--rg", "SM:"+read_group, "-1", seq[0], "-2", seq[1], "-p", threads]
		else:
			bowtie2 = ["bowtie2", "-x", reference, "--rg-id", read_group, "--rg", "SM:"+read_group, "-U", seq[0], "-p", threads]
		return(bowtie2)

	def map_and_process(self, map_command, prefix, threads, final_results_directory):	
		additional_threads = str(int(threads) - 1)	
		samtools_view = ["samtools", "view", "-b", "-@", additional_threads, "-"]	
		mapreads = subprocess.Popen(map_command, stdout = subprocess.PIPE)
		samtools_view_call = subprocess.Popen(samtools_view, stdin = mapreads.stdout, stdout = subprocess.PIPE)
		mapreads.stdout.close()
		samtools_namesort = ["samtools", "sort", "-n", "-@", additional_threads, "-"]
		samtools_namesort_call = subprocess.Popen(samtools_namesort, stdin = samtools_view_call.stdout, stdout = subprocess.PIPE)
		samtools_view_call.stdout.close()
		samtools_fixmate = ["samtools", "fixmate", "-m", "-@", additional_threads, "-", final_results_directory+prefix+".bam"]
		samtools_fixmate_call = subprocess.Popen(samtools_fixmate, stdin = samtools_namesort_call.stdout)
		samtools_namesort_call.stdout.close()
		samtools_fixmate_call.communicate()
		samtools_sort = ["samtools", "sort", "-@", additional_threads, final_results_directory+prefix+".bam"]	
		samtools_sort_call = subprocess.Popen(samtools_sort, stdout = subprocess.PIPE)
		samtools_rmdup = ["samtools", "markdup", "-r", "-@", additional_threads, "-", final_results_directory+prefix+"_snpbac.bam"]
		samtools_rmdup_call = subprocess.Popen(samtools_rmdup, stdin = samtools_sort_call.stdout)
		samtools_sort_call.stdout.close()
		samtools_rmdup_call.communicate()		
		samtools_index = ["samtools", "index", final_results_directory+prefix+"_snpbac.bam"]
		samtools_index_call = subprocess.Popen(samtools_index).wait()
		rm_intermediate = ["rm", final_results_directory+prefix+".bam"]
		rm_intermediate_call = subprocess.Popen(rm_intermediate)
		print("Mapping and processing complete for {0}").format(prefix)
		bamout = final_results_directory+prefix+"_snpbac.bam"
		return(bamout)

