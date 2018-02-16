# SNPBac
## SNP and short indel calling pipeline for bacteria

## Overview

SNPBac is an easy-to-use command-line tool for single nucleotide polymorphism (SNP) and short insertion/deletion (indel)
calling using bacterial whole-genome sequencing (WGS) data. 

The pipeline has 3 major steps:

1. Reads are mapped to a reference genome using either <a href="https://github.com/lh3/bwa">bwa mem</a> or 
<a href="https://github.com/BenLangmead/bowtie2">bowtie2</a>

2. Variants are called using <a href="http://www.htslib.org/doc/samtools.html">samtools/bcftools</a> or 
<a href="https://github.com/ekg/freebayes">freebayes</a>

3. Variants meeting various thresholds (e.g. appropriate minimum quality, not located in a region of recombination identified
using <a href="https://sanger-pathogens.github.io/gubbins/">Gubbins</a>) are retained

Output files include a matrix and fasta file of core SNPs, isolate consensus sequences in fasta format, and SNPs and indels
in vcf format, among others.

Post issues at https://github.com/lmc297/SNPBac/issues

### SNPBac Citation

#### If you found SNPBac or its source code useful, please cite the following:

Carroll, Laura M. SNPBac: SNP and short indel calling pipelines for bacteria. Version 1.0.0. https://github.com/lmc297/SNPBac.

### Additional Citations

#### If you use SNPBac, please also cite:

Cock, Peter J. A., Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg,
Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, and Michiel J. L. de Hoon. 2009.
<a href="https://www.ncbi.nlm.nih.gov/pubmed/19304878"> Biopython: freely available Python tools for 
computational molecular biology and bioinformatics.</a> *Bioinformatics* 25(11): 1422-1423.

Danecek, Petr, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert E. Handsaker,
Gerton Lunter, Gabor T. Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin, and 1000 Genomes Project Analysis Group. 
2011. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137218/">The variant call format and VCFtools. </a>
*Bioinformatics* 27(15): 2156–2158.

Heng, Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin,
and 1000 Genome Project Data Processing Subgroup. 2009. <a href="https://www.ncbi.nlm.nih.gov/pubmed/19505943">
The Sequence Alignment/Map format and SAMtools.</a> *Bioinformatics* 25(16): 2078–2079.

#### If you use SNPBac's default pipeline (read mapping using bwa mem, variant calling using samtools/bcftools, recombination identification and filtration using Gubbins), please also cite the following:

Croucher, Nicholas J., Andrew J. Page, Thomas R. Connor, Aidan J. Delaney, Jacqueline A. Keane, Stephen D. Bentley,
Julian Parkhill, and Simon R. Harris. 2015. <a href="https://www.ncbi.nlm.nih.gov/pubmed/25414349">
Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins.</a>
*Nucleic Acids Research* 43(3): e15.\*

Li, Heng. 2013. <a href="https://arxiv.org/abs/1303.3997">Aligning sequence reads, clone sequences and assembly contigs with
BWA-MEM.</a> arXiv:1303.3997v2 [q-bio.GN]. 

Stamatakis, Alexandros. 2014. <a href="https://www.ncbi.nlm.nih.gov/pubmed/24451623"> RAxML version 8: a tool for 
phylogenetic analysis and post-analysis of large phylogenies.</a> *Bioinformatics* 30(9): 1312–1313.\*

\* If you set -\-remove_recombination to False, you can omit the citations for Gubbins and RAxML

#### If you use the -\-pipeline freebayes option, please also cite:

Garrison E, Marth G. 2012. <a href="https://arxiv.org/abs/1207.3907">Haplotype-based variant detection from short-read 
sequencing.</a> arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012

#### If you use the -\-aligner bowtie2 option, please cite the following in lieu of the BWA-MEM paper (Li, Heng 2013):

Langmead, Ben and Steven L Salzberg. 2012. <a href="https://www.ncbi.nlm.nih.gov/pubmed/22388286">
Fast gapped-read alignment with Bowtie 2. </a> *Nature Methods* 9(4): 357–359.


------------------------------------------------------------------------


## Quick Start

#### Quick Tips

* Make sure you don't have any spaces, commas, or characters that need to be escaped in your file paths and file names

* Make sure the fastq.gz files for your samples have an extension (i.e. all file names end in .fastq or .fastq.gz) 
and a unique sample name before the first period (example: for 2 different bacterial isolates, "sample1" and "sample2",
```sample1.fastq.gz``` and ```sample2.fastq.gz``` are okay because sample1 &#8800; sample2;
```sample.1.fastq.gz``` and ```sample.2.fastq.gz``` are not okay because sample == sample)

* Make sure your output directory doesn't have a snpbac_final_results directory in it from a previous run
  
#### Command structure:
  
```
snpbac -i [list of paths to fastq.gz files] -o [output directory] -r [reference genome fasta] [options...]
```

For help, type `snpbac -h` or `snpbac --help`

For your current version, type `snpbac --version`

#### Sample commands and pipelines:

**Call variants using the default pipeline (simplest command):**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory -r reference.fasta
```

**Call variants using the default pipeline (verbose):**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory -r reference.fasta --aligner bwa --pipeline samtools \
--quality 20 --remove_recombination True --threads 1
```

**Map reads using bowtie2 instead of bwa mem:**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory -r reference.fasta --aligner bowtie2
```

**Call variants with freebayes instead of samtools/bcftools:**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory -r reference.fasta --pipeline freebayes
```

**Map reads with bowtie2, call variants using freebayes, and do not remove recombination using Gubbins:**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory -r reference.fasta --aligner bowtie2 --pipeline freebayes \
--remove_recombination False
```

**Map reads with bwa mem, call variants using samtools/bcftools, and raise vcftools minimum quality threshold for SNPs
and indels to 30:**

```
snpbac -i sample_list.txt -o /path/to/output_directory -r reference.fasta --quality 30
```


------------------------------------------------------------------------


## Installation
### Install SNPBac using Homebrew (macOS users)
  
SNPBac and its dependencies can be installed using <a href="https://brew.sh/">Homebrew</a>.

1. First, install Homebrew, if necessary, by running the following command from your terminal:
  
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

2. Install pip, if necessary, by running the following command from your terminal:

```
sudo easy_install pip
```

3. Install Biopython, if necessary, by running the following command from your terminal:

```
pip install biopython
```
Note: if you don't have permissions, you may need to use sudo:

```
sudo pip install biopython
```

4. Tap brewsci/science, if necessary, by running the following command from your terminal:

```
brew tap brewsci/science
```

5. Tap SNPBac by running the following command from your terminal:
  
```
brew tap lmc297/homebrew-snpbac
```

6. Install SNPBac and its dependencies by running the following command from your terminal:
  
```
brew install snpbac
```

### Download and run SNPBac using source file (macOS and Ubuntu)

1. To run SNPBac, please download and install the following dependencies, if necessary:
  
  <a href="https://www.python.org/downloads/"> Python v. 2.7</a>
  
  <a href="http://biopython.org/DIST/docs/install/Installation.html"> Biopython v. 1.6.9 or higher</a>
  
  <a href="https://github.com/samtools/samtools">SAMtools v. 1.6</a>
  
  <a href="https://github.com/ekg/freebayes">Freebayes v. 1.1.0-dirty</a>
  
  <a href="https://vcftools.github.io/downloads.html">VCFtools v. 0.1.15</a>
  
  <a href="https://github.com/samtools/bcftoolsl">BCFtools v. 0.1.6</a>
  
  <a href="https://sanger-pathogens.github.io/gubbins/">Gubbins v. 2.3.1</a>
  
  <a href="https://github.com/lh3/bwa">BWA v. 0.7.17-r1188</a>
  
  <a href="https://github.com/BenLangmead/bowtie2"> Bowtie 2 v. 2.3.4.1 </a>
  
  <a href="https://github.com/stamatak/standard-RAxML"> RAxML v. 8.2.11</a>
  
2. Download the SNPBac tar.gz file, and store it in your directory of choice:

https://github.com/lmc297/SNPBac/raw/master/archive/snpbac-1.0.0.tar.gz

3. Extract SNPBac

```
tar -xzvf snpbac-1.0.0.tar.gz
```

Note: to ensure that SNPBac works correctly, make sure the files remain in the same directory.

4. To run SNPBac, call Python 2.7 and supply the full path to the SNPBac executable:

```
python /path/to/executable/snpbac [options...]
```

Note: In the examples below, snpbac commands are shown as ```snpbac [options...]```. If you are calling snpbac from 
the source file (i.e. you didn't install snpbac using Homebrew), keep in mind that you may have to call python and supply 
the path to snpbac to execute the program: ```python snpbac [options...]```.


------------------------------------------------------------------------
  
  
## Usage and Options

### Input
  
To call variants using SNPBac, you need the following as input:

1. **Reads from the samples in which you would like to call variants**
These files will be in **fastq** or gzipped fastq (**fastq.gz**) format. Single end reads should be in a single fastq/fastq.gz
file (one per sample), and paired end reads should be separated into two separate fastq/fastq.gz files, one for forward reads 
and one for reverse (two files per sample). 

Importantly, before running SNPBac, **you should make sure your reads are compatible with your selected alignment
method (bwa mem or bowtie2).** If you are unsure whether your sequencing method is compatible with either bwa mem
or bowtie2, check the appropriate manual.

2. **A text file containing the paths to your reads, one sample per line**
This is a simple text file that you will need to generate to tell SNPBac where your samples are located. Each line of this
file should have the absolute path to the fastq file(s) for a single sample. Single-end reads, which have one file per sample,
should have one absolute path per line. Paired-end reads, which have two files per sample, should have two absolute paths per
line: the absolute path to the forward reads for a sample, and the absolute path to the reverse reads for a sample, separated
by a comma (,). **Make sure none of your file paths or file names have spaces in them!**

Example: I sequenced 5 *Bacillus cereus* isolates using an Illumina MiSeq, and I named them sample1, sample2, sample3, sample4, 
and sample5. 3 isolates (sample1, sample2, and sample3) were sequenced with paired-end reads, giving me 2 fastq.gz files per
sample (sample1_R1.fastq.gz, sample1_R2.fastq.gz; sample2_R1.fastq.gz, sample2_R2.fastq.gz; sample3_R1.fastq.gz, 
sample3_R2.fastq.gz). The remaining 2 isolates (sample4 and sample5) were sequenced using single-end reads, giving me
1 fastq.gz file per sample (sample4.fastq.gz; sample5.fastq.gz). All of my fastq.gz files are stored in a directory called
/path/to/reads/. My input file will look like this:
```
/path/to/reads/sample1_R1.fastq.gz,/path/to/reads/sample1_R2.fastq.gz
/path/to/reads/sample2_R1.fastq.gz,/path/to/reads/sample2_R2.fastq.gz
/path/to/reads/sample3_R1.fastq.gz,/path/to/reads/sample3_R2.fastq.gz
/path/to/reads/sample4.fastq.gz
/path/to/reads/sample5.fastq.gz
```

If you don't want to type this up yourself, SNPBac has a handy helper function, *make_snpbac_infile.py*, to help you!
See the **Making a SNPBac Input File** section below for more information.
  
3. **A reference genome** 
This file will be in **fasta** format, and can contain either a single chromosome or multiple sequences. If the 
reference genome contains multiple sequences, as is the case with a draft genome in the form of contigs or scaffolds, the 
sequences will be concatenated in the order in which they appear in the fasta file to form a closed pseudochromosome.


### Required Arguments

SNPBac can be run from your terminal with the following command line:
  
```
snpbac -i [list of paths to fastq.gz files] -o [output directory] -r reference_genome.fasta [options...]
```

Required arguments are:
  
  
**-i/-\-input [string]**
Path to text file containing absolute paths to fastq.gz files, one sample per line. For paired-end reads, absolute paths to
forward read files and absolute paths to reverse read files should be separated by a comma:
```
/absolute/path/to/sample_forward_1.fastq.gz,/absolute/path/to/sample_reverse_2.fastq.gz
```
Absolute paths to single-end read files can be simply printed (do not include commas):
```
/absolute/path/to/sample.fastq.gz
```
See the **Making a SNPBac Input File** section below for more information.

**-o/-\-output [string]**
Path to desired output directory. Specify the path to the output directory where a results directory (snpbac_final_results)
containing output files will be created.

**-r/-\-reference [string]**
Path to bacterial reference genome in fasta format.


### Optional Arguments

Options that can be specified in SNPBac include the following:
  
**-\-aligner [bwa or bowtie2]**
Align reads to reference genome using bwa mem or bowtie2. Default is set to bwa.

**-\-pipeline [samtools or freebayes]**
Call variants using samtools (which calls samtools/bcftools) or freebayes. Both options call samtools for processing 
sam/bam files, but the actual variant calling is done by samtools/bcftools when "samtools" is called, or freebayes when
"freebayes" is called. Default is set to samtools.

**-\-quality [integer]**
Minimum quality threshold of SNPs and indels to retain when filtering variants with vcftools. SNPs and indels above 
this threshold, which corresponds to the estimated probability of polymorphism based on the Phred score,
will be retained and used downstream. Default is set to 20 (99% chance of polymorphism).

**-\-remove_recombination [True or False]**
Remove SNPs found in regions of recombination identified using Gubbins. If set to True, Gubbins will be used to identify 
regions of recombination in the consensus genomes, and SNPBac will filter out SNPs found in these regions when producing
a final matrix and fasta of core SNPs. Setting this to False will bypass Gubbins. Default is set to True.

**-t/-\-threads [integer]**
Number of threads to use. Default is set to 1.


------------------------------------------------------------------------


## Output Directories and Files
  
A single SNPBac run will deposit the following in your specified output directory (-\-output):
  
**snpbac_final_results**
*directory*
Final results directory in which SNPBac deposits all of its output files. SNPBac creates this directory in your specified
output directory (-\-output). 

Note that if you want to re-run an analysis in the sample output directory (-\-output) as a previous analysis, you should
move the snpbac_final_results directory from your previous run to a new directory (not your specified output directory)
to avoid over-writing your previous results.

***your_genome_final_results.txt***
*file*
Final results text file, 1 per input genome. seq2mlst creates this final results text file, which contains a tab-separated 
line, containing the isolate's sequence type (ST), if available, and the corresponding allelic types (AT). The best-matching 
allele is reported at each locus; an allele that does not match with 100\% identity or coverage is denoted by an asterisk 
(\*), while an allele that is not detected in the genome at the given e-value threshold is denoted by "?". If a sequence type
cannot be determined using the best-matching allelic types, a "?" is listed in its place. A ST that is detemined using any 
best-matching alleles that did not match with 100\% identity or coverage is denoted by \*, regardless of whether all 7 
alleles could be associated with a ST or not.

**genefiles**
*directory*
Directory in which seq2mlst deposits genefiles, (multi)fasta files which contain the sequences of all genes detected in 
a run. seq2mlst creates this directory within the seq2mlst_final_results directory within your specified output directory 
(output_directory/seq2mlst_final_results/genefiles).

***some_gene_genefile.fasta***
*file*
seq2mlst genefiles, (multi)fasta files which contain the sequences of all genes detected in a run. 1 file per each allele 
is created files are created, and if seq2mlst is run using more than 1 genome as input (either in multifasta format, or if 
seq2mlst is run in a loop), genes from each genome are aggregated together in each genefile. These files are formatted so you
can easily input them into your favorite aligner, phylogenetic tree construction program, the NCBI BLAST server, etc.

**isolatefiles**
*directory*
Directory in which seq2mlst deposits results directories for individual genomes. seq2mlst creates this directory within the 
seq2mlst_final_results directory within your specified output directory 
(output_directory/seq2mlst_final_results/isolatefiles).

***your_genome_results***
*directory*
Directory in which seq2mlst deposits additional results files for each input genome. seq2mlst creates this directory within 
the isolatefiles directory (output_directory/seq2mlst_final_results/isolatefiles/*your_genome*_results). Within this 
directory, you'll find detailed tab-separated results files for each typing analysis performed, as well as fasta files 
containing genes extracted from the genome in question. 


------------------------------------------------------------------------
