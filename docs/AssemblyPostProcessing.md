# VATK
[Virus Analysis Tool Kit](https://github.com/centre-for-virus-research/VATK)

## Assembly Post-processing

Once a SAM/BAM has been generated, a number of post-processing can be undertaken.

### Extracting statistics from SAM/BAM

An important step for checking the quality of a particular reference assembly is to 
generate depth and breadth of coverage of the reads. A number of tools have been developed 
for this.

#### weeSAM

[weeSAM](https://github.com/centre-for-virus-research/weeSAM) takes either a SAM or BAM file 
as an input and can generate the output in text-tab delimited format or as an html file with 
figures of coverage plots. If a SAM file is provided as input, the BAM file will be automatically produced.

#### SamRemoveIndels

[SamRemoveIndels](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/SamRemoveIndels.awk) is an awk script 
which is used in the HCMV pipeline to remove reads with indels from a SAM file.

```
SamRemoveIndels samfile.sam > samfile_with_without_indels.sam
```


#### UniqSamPE

[UniqSamPE](https://github.com/centre-for-virus-research/VATK/blob/master/AssemblyPostProcessing/UniqSamPE.awk) is awk script 
to remove duplicate paired-end fragments that start and end at the exact same position. This script is used in the HCMV pipeline to 
measure the level of clonality in a given alignment, i.e, it remove reads lacking unique 
fragment coordinates.

```
UniqSamPE.awk samfile.sam > samfile_with_unique_reads.sam
```

