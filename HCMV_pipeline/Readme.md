# VATK
[Virus Analysis Tool Kit](https://github.com/centre-for-virus-research/VATK)

## HCMV pipeline

**alignmentStats_ReCVR.sh** is a bash script written to align multiple fastq files
to a reference sequence. The script expect 4 arguments.

```
alignmentStats_ReCVR.sh <input_dir> <ref.fa> <signature.fa> <lib_name>
```

**Input**

**\<input_dir\>** is a directory of paired fastq files. The fastq files need to be in a directory 
and have the extension \_R1\_001.fastq and \_R2\_001.fastq (or \_R1\_001.fq and \_R2\_001.fq).

**\<ref.fa\>** is the name of the reference fasta file (e.g., merlin.fa)

**\<signature.fa\>** is the fasta file with the signatures of different HCMV strains.

**\<lib_name\>** is the name of the reference library when building the bowtie2 index

### Dependencies

[cutadapt](https://cutadapt.readthedocs.io/en/stable/)

[trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[samtools](https://sourceforge.net/projects/samtools/files/)

[weeSAMv1.4](https://github.com/centre-for-virus-research/weeSAM/blob/master/legacy_versions/weeSAMv1.4)

[gawk](https://www.gnu.org/software/gawk/) 

[SamRemoveIndels.awk](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/SamRemoveIndels.awk) - hash-bang may need to be changed depending on your gawk installation

[UniqSamPE.awk](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/UniqSamPE.awk) - hash-bang may need to be changed depending on your gawk installation

[miRNA_Search](https://github.com/centre-for-virus-research/VATK/tree/master/GenotypingTools)
 
 
### Pipeline

#### Step 1:

Each fastq file in the folder is trimmed using [trim_galore](https://github.com/FelixKrueger/TrimGalore)
 with the following settings (--paired --length 21 --quality 10 --stringency 3).

#### Step 2:

The processed reads are subsequently aligned against the reference provided using bowtie2 
allowing for a maximum fragment length of 1200 (-X 1200)

#### Step 3:
 
The assembly statistics are generated using [weeSAMv1.4](https://github.com/centre-for-virus-research/weeSAM/blob/master/legacy_versions/weeSAMv1.4).
A newer version of weeSAM is available [here](https://github.com/centre-for-virus-research/weeSAM) 
if you wish to have more comprehensive statistics.
A number of assembly statistics are also printed to the terminal.
 
#### Step 4:

The library diversity is also estimated, first using  [SamRemoveIndels.awk](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/SamRemoveIndels.awk) 
and then with [UniqSamPE.awk](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/UniqSamPE.awk).
These provide additional statistics which enable the calculation of the Ratio of total to unique coverage.

The diversity of genotypes in the sample is also estimated using [miRNA_Search](https://github.com/centre-for-virus-research/VATK/tree/master/GenotypingTools), 
which used the signature motifs to determine the number of posible strains in the sample.

#### Step 5:

The final stats are printed to output.csv in the input directory.
