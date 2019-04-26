# VATK
[Virus Analysis Tool Kit](https://github.com/centre-for-virus-research/VATK)



1. Evaluating your raw sequencing data
  * counting sequences: fastq and fastq.gz
  * quality assurance: FASTQC and samstats
  * Determining the quality and complexity of next-generation sequencing data without a reference genome 
2. trimming and de-multiplexing: Trimmomatic, Fastx_toolkit, fastx_trimmer, cutadapt, trim_galore, condetri, flexor (https://github.com/seqan/flexbar)
3. genotype detection, reference selection 
  * [MIRNA\_SEARCH](https://github.com/centre-for-virus-research/VATK/tree/master/GenotypingTools)
4. alignment: bowtie2, bwa, tanoti
5. post alignment stats: 
  * [weeSAM](https://github.com/centre-for-virus-research/weeSAM) - python script for various coverage statistics (legacy versions in perl)
  * [SamRemoveIndels](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/SamRemoveIndels.awk) - awk script to remove reads with indels from a SAM file
  * [UniqSamPE](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/UniqSamPE.awk) - awk script to remove paired-end fragments that start and end at the exact same position
6. post alignment analysis: haplotype window, remove primer reads, detecting chimeras
7. consensus calling 
8. variant calling (link to ValVs)
