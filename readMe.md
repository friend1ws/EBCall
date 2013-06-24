EBCall (Empirical Baysian mutation Calling)
===========

EBCall is a software package for somatic mutation detection (including InDels). EBCall uses not only paired tumor/normal sequence data of a target sample, but also multiple non-paired normal reference samples for evaluating distribution of sequencing errors, which leads to an accurate mutaiton detection even in case of low sequencing depths and low allele frequencies.

Please mail friend1ws@gmail.com for problems or questions.



Paper
----------

"[An empirical Bayesian framework for mutation detection from cancer genome sequencing data](http://nar.oxfordjournals.org/content/early/2013/03/06/nar.gkt126.abstract)", Shiraishi et al.,  
Nucleic Acids Research (First published online: March 6, 2013)

Dependecy
----------

* [samtools](http://samtools.sourceforge.net/)
* [R](http://www.r-project.org/)
* [The VGAM package for R](http://cran.r-project.org/web/packages/VGAM/index.html)


SetUp
----------

1. Download the EBCall package to any directory.

2. Dowload external tools  
 **samtools**: Download and extract and install samtools (we tested ver. 0.1.18) to any directory.  
 **R**: Install R and then locate the R and Rscript executables. Install the VGAM package.

3. Onen config.sh and set each entry.  
 **a reference genome**: the path to the reference genome (.fasta file) to which your sequence data is aligned.(we just test on the hg19 human genome reference)  
 **samtools path**: the path to the samtools executable  
 **R path**: the path to the R script  
 **min_tumor_variant_read**: Mutations whose numbers of variant reads are less than this value (4 for defult) in the target tumor sample are filtered.  
 **min_tumor_allele_freq**: Mutations whose allele frequencies are less than this value (0.1 for default) in the target tumor sample are filtered.  
 **max_normal_allele_freq**: Mutations whose allele frequencies are more than this value (0.1 for defualt) in the paired normal sample are filtered.  
 **min_minus_log10(p-value)**: The threshold of the P-value generated in our method (3 = -log10(0.001) for default).
		


Prepare Input Data
----------

All the input should be indexed bam files.  
**target tumor sample**: the .bam file of the target tumor sample  
**target normal sample**: the .bam file of the paired normal sample. This is used for filtering germline mutations  
**list of normal reference samples**: the list of paths to .bam files for non-paired normal reference samples. Please name the text file as you like (e.g., myNormalRef.txt), and list the paths of .bam files as follows:  

	/home/yshira/ngs/data/sequence/normalreference1.bam
	/home/yshira/ngs/data/sequence/normalreference2.bam
	...
	/home/yshira/ngs/data/sequence/normalreference10.bam
	
We recommend that at least 10 normal samples which are not contaminated by cancer cells are prepared. The paired normal sample may be included to the set of normal reference samples. 
	
How to run
---

Compile C++ programs

	$ make

Just type the following command

	$ sh ebCall_v2.sh <path to the target tumor bam file> <path to the target normal bam file> <path to the output directory> <path to the text file of the list for normal reference samples>

Then you will get the output under the specified output directory.


Output
---

The format of the result is suitable for adding annotation by annovar.  
**Chr, Start, End**: the position of the candidate mutation  
**Ref**: the reference base for that position ('-' for insertions)  	
**Obs**: the alternated sequence for the mutation candidate ('-' for deletions)  
**bases_tumor, base_normal**: (sequencing depth for positive strand: the number of variant reads for positive strand: sequencing depth for negative strand: the number of variant read for negative strand) for both tumor and normal samples  
**misRate_tumor,misRate_normal**: the mismatch rates computed in the tumor and normal samples  
**strandRatio_tumor, strandRatio_normal**: the ratio of variant reads aligned to positive strand for the tumor and normal samples  
**depth_tumor depth_normal**: sequencing depths for that position for the tumor and normal samples  
**variantNum_tumor, variantNum_normal**: the number of supporting variant read in the tumor and normal samples  	
**p-value**: the minus logarithm of p-value of the EBCall. This is a combined value from the two p-values caluculated in positive and negative strands  
**p-value (+strand), p-value (-strand)**: the minus logarithm of p-value of the EBCall for positive and negative strand, respectively  
**p-value(Fisher)**: the minus logarithm of the p-value by Fisher's exact test  
**alpha (+starnd),beta (+strand),alpha (-strand),beta(-strand)**: the estimated parameter values of Beta-Binominal sequencing model for that variant.

Test run
----------
We provide a set of test data files in the EBCall-master/testdata directory and the result in the EBCall-master/testresult directory.   

Edit EBCall/testdata/list_normal_sample.txt to adjust the paths to the EBCall-master directory.

	/home/your_username/EBCall-master/testdata/normalreference1.bam
	/home/your_username/EBCall-master/testdata/normalreference2.bam
	...
	/home/your_username/EBCall-master/testdata/normalreference10.bam

Type the following command after setup EBCall and compiling C++ programs. 

	sh ebCall_v2.sh testdata/tumor.bam testdata/normal.bam testout testdata/list_normal_sample.txt

Result is stored under the testout directory.

When EBCall goes into error, please download chr11.fa from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr11.fa.gz and decompress it.
Then, onen config.sh and set the path to the chr11.fa.

	# path to the reference genome
	PATH_TO_REF=/home/your_username/ref/chr11.fa

Copyright
----------
Copyright (c) 2012, Yuichi Shiraishi, Kenichi Chiba

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
  * We ask you to cite one of the following papers using this software.
  	** "Shiraishi et al., An empirical Bayesian framework for mutation detection from cancer genome sequencing data, Nucleic Acids Research, First published online: March 6, 2013"

"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

