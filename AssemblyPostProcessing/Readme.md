# VATK
[Virus Analysis Tool Kit](https://github.com/centre-for-virus-research/VATK)

Once a SAM/BAM has been generated, a number of post-processing can be undertaken.

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
