import subprocess
import os

class RecombinationFilter:
	"""
	Detect and remove regions of recombination using Gubbins
	"""

	def __init__(self, final_results_directory, threads):
		self.final_results_directory = final_results_directory
		self.threads = threads

	def run_gubbins(self, final_results_directory, threads):
		os.mkdir(final_results_directory+"gubbins")	
		os.system("cat "+final_results_directory+"*_consensus.fasta > "+final_results_directory+"gubbins/snpbac.fna")		
		gubbins = ["run_gubbins.py", "--prefix", final_results_directory+"gubbins/snpbac_gubbins", "--threads", threads, final_results_directory+"gubbins/snpbac.fna"]
		print gubbins
		gubbins_call = subprocess.Popen(gubbins).wait()
