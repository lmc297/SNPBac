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

Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin,
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

* If you get an error regarding Vcf.Pm ("Can't locate Vcf.pm in @INC..") when you try to run SNPBac, make sure you <a href="http://vcftools.sourceforge.net/examples.html">set the environment variable PERL5LIB to include the Vcf.pm module </a> by typing ```export PERL5LIB=/path/to/your/installation/perl``` in your terminal prior to running SNPBac or vcftools (replace "/path/to/your/installation/perl" with the appropriate path for your installation).

For example: if you are a Mac user who installed vcftools using Homebrew, running the following command in your terminal prior to running SNPBac or vcftools should alleviate this problem:

```
export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB}
```
  
#### Command structure:
  
```
snpbac -i [list of paths to fastq.gz files] -o [output directory] -r [reference genome fasta] [options...]
```

For help, type `snpbac -h` or `snpbac --help`

For your current version, type `snpbac --version`

#### Sample commands and pipelines:

**Call variants using the default pipeline (simplest command):**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory/ -r reference.fasta
```

**Call variants using the default pipeline (verbose):**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory/ -r reference.fasta --aligner bwa --pipeline samtools \
--quality 20 --remove_recombination True --threads 1
```

**Map reads using bowtie2 instead of bwa mem, and use 8 threads:**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory/ -r reference.fasta --aligner bowtie2 -t 8
```

**Call variants with freebayes instead of samtools/bcftools:**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory/ -r reference.fasta --pipeline freebayes
```

**Map reads with bowtie2, call variants using freebayes, and do not remove recombination using Gubbins:**
  
```
snpbac -i sample_list.txt -o /path/to/output_directory/ -r reference.fasta --aligner bowtie2 --pipeline freebayes \
--remove_recombination False
```

**Map reads with bwa mem, call variants using samtools/bcftools, and raise vcftools minimum quality threshold for SNPs
and indels to 30:**

```
snpbac -i sample_list.txt -o /path/to/output_directory/ -r reference.fasta --quality 30
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

Note: if you are a Mac user who installed SNPBac or vcftools using Homebrew (as described above), and you get an error regarding Vcf.Pm ("Can't locate Vcf.pm in @INC..") when you try to run SNPBac or vcftools, try running the following command in your terminal prior to running SNPBac or vcftools:

```
export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB}
```

### Download and run SNPBac using source file (macOS and Ubuntu)

1. To run SNPBac, please download and install:

A. All of the following dependencies, if necessary:
  
<a href="https://www.python.org/downloads/"> Python v. 2.7</a>
  
<a href="http://biopython.org/DIST/docs/install/Installation.html"> Biopython v. 1.6.9 or higher</a>
  
<a href="https://github.com/samtools/samtools">SAMtools v. 1.6</a>
  
<a href="https://vcftools.github.io/downloads.html">VCFtools v. 0.1.15</a>

B. One or more of the following dependencies for read alignment (SNPBac's defalt is BWA):

<a href="https://github.com/lh3/bwa">BWA v. 0.7.17-r1188</a>

<a href="https://github.com/BenLangmead/bowtie2"> Bowtie 2 v. 2.3.4.1 </a>

C. One or more of the following dependencies for variant calling, if necessary (SNPBac's default is BCFtools, used in combination with SAMtools):

<a href="https://github.com/samtools/bcftoolsl">BCFtools v. 0.1.6</a>

<a href="https://github.com/ekg/freebayes">Freebayes v. 1.1.0-dirty</a>
  
D. Both of the following dependencies, if you would like to filter out SNPs identified in regions of recombination (SNPBac will query these programs by defaut):

<a href="https://github.com/stamatak/standard-RAxML"> RAxML v. 8.2.11</a>

<a href="https://sanger-pathogens.github.io/gubbins/">Gubbins v. 2.3.1</a>

To ensure that SNPBac runs properly, **make sure these programs are in your path.**
  
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

Additional note: If you get an error regarding Vcf.Pm ("Can't locate Vcf.pm in @INC..") when you try to run SNPBac, make sure you <a href="http://vcftools.sourceforge.net/examples.html">set the environment variable PERL5LIB to include the Vcf.pm module </a> by typing ```export PERL5LIB=/path/to/your/installation/perl``` in your terminal prior to running SNPBac or vcftools (replace "/path/to/your/installation/perl" with the appropriate path for your installation).

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
or bowtie2, consult the appropriate manual.

2. **A text file containing the paths to your reads, one sample per line**
This is a simple text file that you will need to generate to tell SNPBac where your samples are located. Each line of this
file should have the absolute path to the fastq/fastq.gz file(s) for a single sample. Single-end reads, which have one file per sample, should have one absolute path per line. Paired-end reads, which have two files per sample, should have two absolute paths per line: the absolute path to the forward reads for a sample, and the absolute path to the reverse reads for a sample, separated by a comma (,). **Make sure none of your file paths or file names have spaces, commas, or other special characters in them!**

Example: I sequenced 5 *Bacillus cereus* isolates using an Illumina MiSeq, and I named them sample1, sample2, sample3, sample4, and sample5. 3 isolates (sample1, sample2, and sample3) were sequenced with paired-end reads, giving me 2 fastq.gz files per sample, with forward reads ending in "_R1.fastq.gz" and reverse reads ending in "_R2.fastq.gz" (sample1_R1.fastq.gz, sample1_R2.fastq.gz; sample2_R1.fastq.gz, sample2_R2.fastq.gz; sample3_R1.fastq.gz, sample3_R2.fastq.gz). The remaining 2 isolates (sample4 and sample5) were sequenced using single-end reads, giving me 1 fastq.gz file per sample (sample4.fastq.gz; sample5.fastq.gz). All of my fastq.gz files are stored in a directory called /path/to/reads/. My input file will look like this:
```
/path/to/reads/sample1_R1.fastq.gz,/path/to/reads/sample1_R2.fastq.gz
/path/to/reads/sample2_R1.fastq.gz,/path/to/reads/sample2_R2.fastq.gz
/path/to/reads/sample3_R1.fastq.gz,/path/to/reads/sample3_R2.fastq.gz
/path/to/reads/sample4.fastq.gz
/path/to/reads/sample5.fastq.gz
```

If you don't want to type this up yourself, you may be able to use SNPBac's has a handy helper function, *make_snpbac_infile.py*. See the **Making a SNPBac Input File** section below for more information.
  
3. **A reference genome** 
This file will be in **fasta** format, and can contain either a single chromosome or multiple sequences. If the 
reference genome contains multiple sequences, as is the case with a draft genome in the form of contigs or scaffolds, the 
sequences will be concatenated in the order in which they appear in the fasta file to form a closed pseudochromosome.


### Required Arguments

SNPBac can be run from your terminal with the following command line:
  
```
snpbac -i [list of fastq.gz paths] -o [output directory] -r reference_genome.fasta [options...]
```

Required arguments are:
  
  
**-i/-\-input [string]**
Path to text file containing absolute paths to fastq.gz files, one sample per line. For paired-end reads, absolute paths to
forward read files and absolute paths to reverse read files should be separated by a comma:
```
/absolute/path/to/sample_forward_1.fastq.gz,/absolute/path/to/sample_reverse_2.fastq.gz
```
Absolute paths to single-end read files can be simply printed (make sure there are no commas in the path or file name):
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
regions of recombination in consensus genomes for all samples, and SNPBac will filter out SNPs found in these regions when producing a final matrix and fasta of core SNPs. Setting this to False will bypass Gubbins. Default is set to True (Note that Gubbins requires at least 4 sequences; if you are calling variants in fewer than 4 samples, set this option to False).

**-t/-\-threads [integer]**
Number of threads to use. Default is set to 1.


------------------------------------------------------------------------


## Output Directories and Files
  
A single SNPBac run will deposit the following in your specified output directory (-\-output):
  
**snpbac_final_results**

*1 directory*

Final results directory in which SNPBac deposits all of its output files. SNPBac creates this directory in your specified
output directory (-\-output). 

Note that if you want to re-run an analysis in the same output directory (-\-output) as a previous analysis (i.e. there is already a directory called "snpbac_final_results" in your specified -\-output directory), you should
rename the old directory to something else, or move the snpbac_final_results directory from your previous run to a new directory.

***snpbac_core_snps.fasta***

*1 file*

A fasta file containing core SNPs (i.e. SNPs at sites present in all of your input samples), with one sequence per sample (the reference is excluded). 

Core SNPs in this file are above the minimum quality threshold specified by the -\-quality argument (default is 20). If ```--remove_recombination True``` is used (as is done by default), SNPs in regions of recombination identified using Gubbins are excluded from this file. In addition, SNPs that are called relative to the reference genome, but not within the samples themselves (e.g. the reference genome has a T, and all of the samples have a C; the site is a SNP relative to the reference, but the site is invariant when the samples are compared to each other) are also excluded from this file, as these would be considered invariant when just the samples are taken into consideration and not the reference.

SNPs in this file are the same as those in the *snpbac_core_snps.txt* and *snpbac_core_snps.vcf* files (below), just in fasta format.

This file can be used to build a phylogenetic tree with your favorite software, or to get pairwise SNP differences between samples (all sites are variant sites).

***snpbac_core_snps.txt***

*1 file*

A tab-separated file containing a matrix of core SNPs (i.e. SNPs at sites present in all of your input samples), with one column per sample (reference excluded) and one row per position in the reference genome at which the core SNP was detected. If the reference genome used contained multiple sequences, as is the case with a draft genome in the form of contigs or scaffolds, the sequences were concatenated in the order in which they appeared in the fasta file to form a closed pseudochromosome, and the "position" column corresponds to the position in this pseudochromosome at which a core SNP was found.

Core SNPs in this file are above the minimum quality threshold specified by the -\-quality argument (default is 20). If ```--remove_recombination True``` is used (as is done by default), SNPs in regions of recombination identified using Gubbins are excluded from this file. In addition, SNPs that are called relative to the reference genome, but not within the samples themselves (e.g. the reference genome has a T, and all of the samples have a C; the site is a SNP relative to the reference, but the site is invariant when the samples are compared to each other) are also excluded from this file, as these would be considered invariant when just the samples are taken into consideration and not the reference.

SNPs in this file are the same as those in the *snpbac_core_snps.fasta* file (above) and the *snpbac_core_snps.vcf* file (below), just in matrix format, with information on their position included.

This file can be easily viewed in Excel or loaded into R.

***snpbac_core_snps.vcf***

*1 file*

A variant call format (VCF) file containing core SNPs (i.e. SNPs at sites present in all of your input samples). If the reference genome used contained multiple sequences, as is the case with a draft genome in the form of contigs or scaffolds, the sequences were concatenated in the order in which they appeared in the fasta file to form a closed pseudochromosome, and the "POS" column corresponds to the position in this pseudochromosome at which a core SNP was found.

Core SNPs in this file are above the minimum quality threshold specified by the -\-quality argument (default is 20). If ```--remove_recombination True``` is used (as is done by default), SNPs in regions of recombination identified using Gubbins are excluded from this file. In addition, SNPs that are called relative to the reference genome, but not within the samples themselves (e.g. the reference genome has a T, and all of the samples have a C; the site is a SNP relative to the reference, but the site is invariant when the samples are compared to each other) are also excluded from this file, as these would be considered invariant when just the samples are taken into consideration and not the reference.

SNPs in this file are the same as those in the *snpbac_core_snps.fasta* and *snpbac_core_snps.txt* files (above), just in VCF format.

This file can be used with vcftools/bcftools, or viewed in a text editor or Excel.

***snpbac_snps.recode.vcf***

*1 file*

A variant call format (VCF) file containing SNPs above the minimum quality threshold specified by the -\-quality argument (default is 20). In addition to containing core SNPs, this VCF file also contains SNPs present at sites in the genome that are not present in all samples (i.e. it contains SNPs that are not part of the core genome of all of the samples), as well as sites that are SNPs relative to the reference genome, but not within the samples themselves (e.g. the reference genome has a T, and all of the samples have a C; the site is a SNP relative to the reference, and is treated as such in this file). This file also contains SNPs that *ARE* present in regions of recombination (this step is performed prior to running Gubbins, so regions of recombination have not been filtered out).

If the reference genome used contained multiple sequences, as is the case with a draft genome in the form of contigs or scaffolds, the sequences were concatenated in the order in which they appeared in the fasta file to form a closed pseudochromosome, and the "POS" column corresponds to the position in this pseudochromosome at which a SNP was found.

This file can be used with vcftools/bcftools, or viewed in a text editor or Excel.

This file has an assoicated log file, *snpbac_snps.log*, which explains the command used by vcftools to filter SNPs.

This file also has an associated index file, *snpbac_snps.recode.vcf.gz.tbi*, which can be deleted if desired. 

***snpbac_indels.recode.vcf***

*1 file*

A variant call format (VCF) file containing insertions and deletions (indels) above the minimum quality threshold specified by the -\-quality argument (default is 20). 

If the reference genome used contained multiple sequences, as is the case with a draft genome in the form of contigs or scaffolds, the sequences were concatenated in the order in which they appeared in the fasta file to form a closed pseudochromosome, and the "POS" column corresponds to the position in this pseudochromosome at which an indel was found.

This file can be used with vcftools/bcftools, or viewed in a text editor or Excel.

This file has an assoicated log file, *snpbac_indels.log*, which explains the command used by vcftools to filter indels.

***snpbac_raw.vcf***

*1 file*

A variant call format (VCF) file containing all SNPs and indels called using the specified variant caller (samtools/bcftools or freebayes; default is samtools/bcftools). This file contains all raw variants and has undergone no filtering (this file is created prior to running Gubbins, as well as prior to quality filtering using vcftools; it is the raw VCF file produced by samtools/bcftools or freebayes).

If the reference genome used contained multiple sequences, as is the case with a draft genome in the form of contigs or scaffolds, the sequences were concatenated in the order in which they appeared in the fasta file to form a closed pseudochromosome, and the "POS" column corresponds to the position in this pseudochromosome at which a SNP or indel was found.

This file can be used with vcftools/bcftools, or viewed in a text editor or Excel.

***your_sample_consensus.fasta***

*1 file per sample*

A fasta file containing the consensus sequence of a sample. The consensus sequence is generated by applying the SNPs found in *snpbac_snps.recode.vcf* (i.e. all SNPs above the quality threshold, including non-core SNPs) to the reference genome (this step is performed prior to running Gubbins, so regions of recombination have not been filtered out). These files are concatenated and used as input to Gubbins.

The "your_sample" prefix for your_sample_consensus.fasta is named by using everything prior to the first "." in a sample's .fastq/.fastq.gz file. If a sample uses paired-end reads, the filename prefix is generated using the name of the file containing the forward reads, so "your_sample" will probably be something like "your_sample_1" or "your_sample_R1", depending on the name of the forward read file (the data in the file itself is created using both forward and reverse reads, though, so don't be alarmed!).

In addition to a consensus fasta, each sample will have (i) an individual gzipped vcf file *your_sample.vcf.gz*, and (ii) an associated index file *your_sample.vcf.gz.tbi*. These are produced from *snpbac_snps.recode.vcf* using the vcf-subset command found in vcftools. These files can be deleted, if desired.

***your_sample_snpbac.bam***

*1 file per sample*

A sorted bam file used for variant calling, 1 per input sample. This bam file is created by mapping reads to your reference genome using either bwa mem or bowtie2, filtering out PCR duplicates, and sorting, and serves as the final input into the selected variant caller (samtools/bcftools or freebayes). 

These files may be useful for looking at mapping statistics, and can be easily converted to the human-readable SAM format using samtools (note that SAM files can get very large):
```
samtools view -h -o your_sample.sam your_sample_snpbac.bam
```

If you are not planning on using these files further, they may be deleted, as they can get large. 

In addition to a sorted bam, each sample will have an index file, *your_sample_snpbac.bam.bai*, which may be deleted if you're not planning on manipulating the bam files any further.

***samtools_error_log.txt or freebayes_error_log.txt***

*1 file*

If something goes wrong during variant calling, it will be logged here. 

***your_reference_genome_pseudochrom.fasta***

*1 file*

If the reference genome you supply to SNPBac contains more than one sequence (e.g. it is a draft genome made up of contigs or scaffolds), SNPBac concatenates the sequences into a pseudochromosome in the order of their appearance. This is the fasta file for the concatenated pseudochromosome. Any position coordinates given by vcf files or snp matrices output by SNPBac will refer to positions in this file.

In addition to the actual fasta, several index files are output, depending on which pipeline you run.

**gubbins**

*directory*

If ```--remove_recombination True``` is used (as is done by default), this directory is created in the snpbac_final_results directory within the specified output directory (output_directory/snpbac_final_results/gubbins). This directory contains the input file that SNPBac passes to Gubbins, titled *snpbac.fna*, which is a multifasta with all of the concatenated consensus sequences output by SNPBac (feel free to delete this file and/or the snpbac consensus fasta files for each of your samples if you do not plan on using them any more, as they are redundant).

This directory contains many useful and interesting files that Gubbins outputs; if you are interested in all SNPs above a certain quality threshold (rather than just the core SNPs that SNPBac outputs), but you still want to exclude SNPs in regions of recombination, check the files in this directory. Some examples of output files that Gubbins produces include a phylogenetic tree constructed using RAxML, recombination predictions, and a fasta file constructed using all SNPs above the quality threshold that are not present in regions of recombination (including the non-core ones that SNPBac filters out). You can read more about the types of files Gubbins produces <a href="https://sanger-pathogens.github.io/gubbins/"> here. </a>


------------------------------------------------------------------------


## Making a SNPBac Input File

As mentioned above, SNPBac requires you to produce a text file containing the absolute paths to your reads, one sample per line. Obviously, you can open up your favorite text editor and type/copy/paste this out, but you may find the following tricks to be useful.

### Single-End Reads
For single-end reads, this is simply a list of absolute paths. You can generate this simply by moving all of your single-end read files to one directory (we'll call it "my_directory"), moving to that directory using ```cd /path/to/my_directory/```, and typing the following command in your terminal (this command works if all of your read files end in ".fastq.gz"; if your read files end in ".fastq", just omit the ".gz" part):

```
ls -d -1 $PWD/*.fastq.gz > se_list.txt
```
This will produce an input file named "se_list.txt" within "my_directory" (but feel free to name it anything you want!). If "my_directory" contains 5 single-end read files ending in .fastq.gz (sample1.fastq.gz ... sample5.fastq.gz), our input file "se_list.txt" will look something like this:

```
/path/to/my_directory/sample1.fastq.gz
/path/to/my_directory/sample2.fastq.gz
/path/to/my_directory/sample3.fastq.gz
/path/to/my_directory/sample4.fastq.gz
/path/to/my_directory/sample5.fastq.gz
```
If we want to run SNPBac with these samples as input, we can now give this to SNPBac's ```-i``` argument:
```
snpbac -i /path/to/my_directory/se_list.txt -o /path/to/any/output/directory/ -r /path/to/reference_genome.fasta
```

### Paired-End Reads
For paired-end reads, which have two files per sample, each line of the input file should have the absolute path to the forward reads for a sample, and the absolute path to the reverse reads for a sample, separated by a comma (,). This is done to ensure that SNPBac knows exactly which forward and reverse read pairs belong together. This can get a little tricky/tedious if your read pairs have different prefixes or suffixes; however, if your read pairs for each sample have identical sample names, as well as the same extensions indicating forward or reverse reads, you can use the handy ```make_snpbac_infile.py``` command.

**Example files names that WILL work with make_snpbac_infile.py:**

sample1_1.fastq.gz, sample1_2.fastq.gz, sample42_1.fastq.gz, sample42_2.fastq.gz, sampleABCD_1.fastq.gz, sampleABCD_2.fastq.gz

All of these files have (i) a shared suffix for forward reads ("_1.fastq.gz"), (ii) a shared suffix for reverse reads ("_2.fastq.gz"), (iii) matching prefixes (sample1==sample1, sample42==sample42, sampleABCD==sampleABCD)

**Example files names that WILL work with make_snpbac_infile.py:**

spam_R1_001.fastq.gz, spam_R2_001.fastq.gz, ham_R1_001.fastq.gz, ham_R2_001.fastq.gz, eggs_R1_001.fastq.gz, eggs_R2_001.fastq.gz

All of these files have (i) a shared suffix for forward reads ("_R1_001.fastq.gz"), (ii) a shared suffix for reverse reads ("_R2_001.fastq.gz"), (iii) matching prefixes (spam==spam, ham==ham, eggs==eggs)

**Example files names that will NOT work with make_snpbac_infile.py:**

sample1_1.fastq.gz, sample1_2.fastq.gz, sample2_R1_001.fastq.gz, sample2_R2_001.fastq.gz

This samples do not have a shared suffix for forward reads ("_1.fastq.gz" &#8800; "_R1_001.fastq.gz") or reverse reads ("_2.fastq.gz &#8800; "_2_001.fastq.gz")

**Example files names that will NOT work with make_snpbac_infile.py:**

sample1_S14_R1.fastq.gz, sample1_S3_R2.fastq.gz, sample2_H4_R1.fastq.gz, sample2_H7_R2.fastq.gz

While these samples have a shared suffix for both forward reads ("_R1.fastq.gz"=="_R1.fastq.gz") and reverse reads ("_R2.fastq.gz"=="_R2.fastq.gz"), each sample does not have a shared prefix ("sample1_S14" &#8800; "sample1_S3"; "sample2_H4" &#8800; "sample2_H7")

If your file names meet these criteria, make sure all of your paired-end read files are in a directory together (we'll call it "my_directory_pe"), and run make_snpbac_infile.py, substituting "_forward_suffix.fastq.gz" and "_reverse_suffix.fastq.gz" with the appropriate forward and reverse suffixes for your sequences, respectively (you can type ```make_snpbac_infile.py -h``` or ```make_snpbac_infile.py --help``` for a help message in your terminal):

```
make_snpbac_infile.py --input /path/to/my_directory_pe/ --out /path/to/my_directory_pe/pe_list.txt \
--forward "_forward_suffix.fastq.gz" --reverse "_reverse_suffix.fastq.gz"
```
This will output our desired input file for paired end reads, titled "pe_list.txt", in our "my_directory_pe" (obviously, you can change the name/path to the output file!). Common forward/reverse suffix pairs you might run into are "_1.fastq.gz"/"_2.fastq.gz", "_R1.fastq.gz"/"_R2.fastq.gz"

If we had 3 paired-end samples in our "my_directory_pe" directory (sampleA_R1.fastq.gz, sampleA_R2.fastq.gz, sampleB_R1.fastq.gz, sampleB_R2.fastq.gz, sampleC_R1.fastq.gz, sampleC_R2.fastq.gz), our output file should look something like this:
```
/path/to/my_directory_pe/sampleA_R1.fastq.gz,/path/to/my_directory_pe/sampleA_R2.fastq.gz
/path/to/my_directory_pe/sampleB_R1.fastq.gz,/path/to/my_directory_pe/sampleB_R2.fastq.gz
/path/to/my_directory_pe/sampleC_R1.fastq.gz,/path/to/my_directory_pe/sampleC_R2.fastq.gz
```
If we want to run SNPBac with these samples as input, we can now give this to SNPBac's ```-i``` argument:
```
snpbac -i /path/to/my_directory_pe/pe_list.txt -o /path/to/any/output/directory/ -r /path/to/reference_genome.fasta
```

### Mixture of Single- and Paired-End Samples

If we have some samples that have single-end reads and others that have paired-end reads, we can just concatenate our lists (se_list.txt and pe_list.txt) from the two examples above using ```cat``` in our terminal:

```
cat /path/to/my_directory/se_list.txt /path/to/my_directory_pe/pe_list.txt > /path/to/any_directory/combined_list.txt
```
This produces a concatenated list called combined_list.txt in a directory called "any_directory" that looks like this (as always, you can change the path and name of the output list):
```
/path/to/my_directory/sample1.fastq.gz
/path/to/my_directory/sample2.fastq.gz
/path/to/my_directory/sample3.fastq.gz
/path/to/my_directory/sample4.fastq.gz
/path/to/my_directory/sample5.fastq.gz
/path/to/my_directory_pe/sampleA_R1.fastq.gz,/path/to/my_directory_pe/sampleA_R2.fastq.gz
/path/to/my_directory_pe/sampleB_R1.fastq.gz,/path/to/my_directory_pe/sampleB_R2.fastq.gz
/path/to/my_directory_pe/sampleC_R1.fastq.gz,/path/to/my_directory_pe/sampleC_R2.fastq.gz
```

We can now give this to SNPBac's ```-i``` argument when we run the command:
```
snpbac -i /path/to/any_directory/combined_list.txt -o /path/to/any/output/directory/ -r /path/to/reference_genome.fasta
```

------------------------------------------------------------------------


## References

Cock, Peter J. A., Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg,
Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, and Michiel J. L. de Hoon. 2009.
<a href="https://www.ncbi.nlm.nih.gov/pubmed/19304878"> Biopython: freely available Python tools for 
computational molecular biology and bioinformatics.</a> *Bioinformatics* 25(11): 1422-1423.

Croucher, Nicholas J., Andrew J. Page, Thomas R. Connor, Aidan J. Delaney, Jacqueline A. Keane, Stephen D. Bentley,
Julian Parkhill, and Simon R. Harris. 2015. <a href="https://www.ncbi.nlm.nih.gov/pubmed/25414349">
Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins.</a>
*Nucleic Acids Research* 43(3): e15.

Danecek, Petr, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert E. Handsaker,
Gerton Lunter, Gabor T. Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin, and 1000 Genomes Project Analysis Group. 
2011. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137218/">The variant call format and VCFtools. </a>
*Bioinformatics* 27(15): 2156–2158.

Garrison E, Marth G. 2012. <a href="https://arxiv.org/abs/1207.3907">Haplotype-based variant detection from short-read 
sequencing.</a> arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012

Langmead, Ben and Steven L Salzberg. 2012. <a href="https://www.ncbi.nlm.nih.gov/pubmed/22388286">
Fast gapped-read alignment with Bowtie 2. </a> *Nature Methods* 9(4): 357–359.

Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin,
and 1000 Genome Project Data Processing Subgroup. 2009. <a href="https://www.ncbi.nlm.nih.gov/pubmed/19505943">
The Sequence Alignment/Map format and SAMtools.</a> *Bioinformatics* 25(16): 2078–2079.

Li, Heng. 2013. <a href="https://arxiv.org/abs/1303.3997">Aligning sequence reads, clone sequences and assembly contigs with
BWA-MEM.</a> arXiv:1303.3997v2 [q-bio.GN]. 

Stamatakis, Alexandros. 2014. <a href="https://www.ncbi.nlm.nih.gov/pubmed/24451623"> RAxML version 8: a tool for 
phylogenetic analysis and post-analysis of large phylogenies.</a> *Bioinformatics* 30(9): 1312–1313.


------------------------------------------------------------------------


Disclaimer: SNPBac is...pretty neat! However, no tool is perfect, and SNPBac and its dependencies are not guaranteed to predict variants with perfect accuracy. As always, interpret your results with caution. We are not responsible for taxonomic misclassifications, misclassifications of an isolate's pathogenic or antimicrobial resistance potential, and/or misinterpretations (biological, statistical, or otherwise) of SNPBac results.




