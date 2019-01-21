#!/bin/bash
#===================================================================
# Maha Maabar, 27/01/2017
# modified by Joseph Hughes, 03/12/2018
#
# A script to calculate statistical information from 
# the results of aligning FASTQ files to the Merlin 
# reference genome and different HCMV strains.
#
# Usage:
# alignmentStats.sh inputDir  merlinRef.fasta signatureFile.fa
#
# <inputDir> Path to the directory which contains the paired-end 
# FASTQ files, the Merlin complete reference genome <merlinRef.fasta>,
# the list of strain genotype signatures <signatureFile.fa>,
# and the name of library to do the alignment <merlin_ref>
#
# Output: a file called "output.csv" in <inputDir>
#=====================================================

usage='echo -e "\n Usage: alignmentStats_ReCVR.sh <input_dir> <ref.fa> <signature.fa> <lib_name>
		<input_dir> is a directory of paired fastq files. The fastq files need to be in a directory and have the extension _R1_001.fastq and _R2_001.fastq (or _R1_001.fq and _R2_001.fq).
		<ref.fa> is the name of the reference fasta file (e.g., merlin.fa)
		<signature.fa> is the fasta file with the signatures of different HCMV strains.
		<lib_name> is the name of the reference library when building the bowtie2 index \n"';

if [[ ! $4 ]] 
then
	printf "${usage}\n\n";
exit;
fi

type cutadapt >/dev/null 2>&1 || { echo >&2 "I require cutadapt but it's not installed or not in the PATH system variable.  Aborting.";exit 0; }
type trim_galore >/dev/null 2>&1 || { echo >&2 "I require trim_galore but it's not installed or not in the PATH system variable.  Aborting.";exit 0; }
type samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed or not in the PATH system variable.  Aborting.";exit 0; }


#INPUT PARAMETERS 
inputDir="$1"
refName="$2"
sigName="$3"
libName="$4"

#create a temp directory in current directory to hold intermediate files
tmpDir="./tmpDir"
mkdir tmpDir
chmod 775 tmpDir

#create the output file with the headings info
labelFile="$tmpDir"/labels.csv
summaryFile="$tmpDir"/output.csv

### preparing the output summary file

#headings for alignment with Merlin reference genome
echo " ,SUMMARY ANALYSIS" > "$labelFile"
echo " ,Dataset name" >> "$labelFile"
echo " ,Total trimmed reads (no.)" >> "$labelFile"
echo " ,HCMV reads (no.)" >> "$labelFile"
echo " ,HCMV reads (%)" >> "$labelFile"
echo " ,Coverage (reads/nt)" >> "$labelFile"
echo "Merlin sequence, Genome size (nt)" >> "$labelFile"
echo " ,Coverage (unique fragment reads/nt)" >> "$labelFile"
echo " ,Coverage total/unique" >> "$labelFile"
echo " ,minimum Coverage" >> "$labelFile"
echo " ,maximun Coverage" >> "$labelFile"
echo " ,Coverage (standard Deviation)" >> "$labelFile"
echo " ,Coverage (> 0.2*Depth)" >> "$labelFile"
echo " ,Coverage (> 1.0*Depth)" >> "$labelFile"
echo " ,Coverage (> 1.8*Depth)" >> "$labelFile"
echo " ,Coefficient of Variation (stdDev/mean)" >> "$labelFile"
echo " ,Number of strains" >> "$labelFile"

#Merlins strains 
echo " ,SIG ANALYSIS" >> "$labelFile"
echo "Genotype,Sig sequence" >> "$labelFile" 

while read line
do 
   if [[ ${line:0:1} == ">" ]] ; then 
      strain=$(echo "$line" |sed 's/^.//')
      read line
      echo $strain","$line  >> "$labelFile"    
   fi
   
done <$sigName
#}


echo $sigName

#Prepare library for alignment
refFile="$inputDir"/"$refName"

#echo bowtie2-build $refFile
echo bowtie2-build $refFile $libName
bowtie2-build $refName $libName &> /dev/null
cp $libName.* $tmpDir

echo "##################################################"
#Main Process 
k=1

for file in "$inputDir"/*R1* #reads files must have .fastq extension
do
  
  echo "File Name ($k): "$file""
  echo "##################################################"
  #get the extension of the file
  fileExtension="${file##*.}"
  
  if [[ "$fileExtension" == "fastq" ]]
  then 
     file2=${file%_R1*}_R2_001.fastq 
     
  fi
  if  [[ "$fileExtension" == "fq" ]] 
  then 
     file2=${file%_R1*}_R2_001.fq 
     
  fi 
  
    
  #copy both files to tmp dir to process
  cp {$file,$file2} "$tmpDir"
  fileName1=${file##*/} 
  fileName2=${file2##*/} 
  echo "Files to process: $fileName1 and $fileName2  "
   
  #get the sampleID from the file names
  basename=${fileName1%.*} 
  sampleID=${basename%_S*} 
  echo "Sample ID: $sampleID "
 


  #Step1: preprocess reads files and concatenate both files
  #trim_galore
  echo "Trimming reads files ..."
  echo "trim_galore --paired --length 21 --quality 10 --stringency 3 "$file" "$file2" &> /dev/null"
  trim_galore --paired --length 21 --quality 10 --stringency 3 "$file" "$file2" &> /dev/null

  trimFile1=${fileName1%.*}"_val_1.fq"
  trimFile2=${fileName2%.*}"_val_2.fq"
  
  mv $trimFile1 $tmpDir/. 
  mv $trimFile2 $tmpDir/.
  rm *_report.txt
  
  catFile="$tmpDir"/"$sampleID""_catfile.fq" #concatenate trimmed files
  cat "$tmpDir"/"$trimFile1" "$tmpDir"/"$trimFile2" > "$catFile"
  echo "Trimmed and concatenated file: $catFile"

  #get the total number of reads 
  numReads=$(wc -l "$catFile" | cut -d ' ' -f1) 
  totalReads=$((numReads / 4))
  echo "Total Number of Reads $totalReads" 

  #Step2: assemble trimmed files to Merlin
  #bowtie2 assembly
  echo "Assembling to Merlin ..."
  assmbFile="$tmpDir"/${sampleID}"_alignmentFile.sam"
  echo bowtie2 -x "$tmpDir/$libName" -1 "$tmpDir"/"$trimFile1" -2 "$tmpDir"/"$trimFile2" -X 1200 -p 8 -S "$assmbFile" 
  bowtie2 -x "$tmpDir"/"$libName" -1 "$tmpDir"/"$trimFile1" -2 "$tmpDir"/"$trimFile2" -X 1200 -p 8 -S "$assmbFile" &> /dev/null

  #convert sam to bam, sort bam and index it
  echo "converting to bam files ..."
  bamFile="$tmpDir"/${sampleID}"_alignmentFile.bam"    
  samtools view -bS $assmbFile | samtools sort -o $bamFile &> /dev/null
  samtools index $bamFile &> /dev/null

  #Step3: get stats from assembly- reference length,#mappedReads, %mapped Reads, average depth
  assmbStats="$tmpDir"/${sampleID}"_assmbStats.txt"
  weeSAMv1.4 -b $bamFile -out $assmbStats  &> /dev/null

  statsLine=$(awk 'NR==2{print;exit}' $assmbStats)
  refLen=$(awk 'NR==2{print;exit}' $assmbStats |cut -f2)
  mappedReads=$(awk 'NR==2{print;exit}' $assmbStats |cut -f3)
  prcMappedReads=$(bc -l <<< "scale = 2; (($mappedReads/$totalReads*100))")
  minDepth=$(awk 'NR==2{print;exit}' $assmbStats |cut -f6)
  maxDepth=$(awk 'NR==2{print;exit}' $assmbStats |cut -f7)
  avrgDepth=$(awk 'NR==2{print;exit}' $assmbStats |cut -f8)
  echo $avrgDepth
  avrgDepthRd=$(echo $avrgDepth |awk '{printf "%.0f", $1}')
  echo "Average Depth =  $avrgDepthRd" 

  stdDev=$(echo $statsLine |cut -d ' ' -f9)
  
  coverage1=$(echo $statsLine |cut -d ' ' -f10) #proportion of num of sites with depth > 0.2*AverageDepth
  coverage2=$(echo $statsLine |cut -d ' ' -f12) #proportion of num of sites with depth > 1.0*AverageDepth
  coverage3=$(echo $statsLine |cut -d ' ' -f13) #proportion of num of sites with depth > 1.8*AverageDepth
  
  varCoef=$(echo $statsLine |cut -d ' ' -f14)  #coefficient of variation

  #Step4: get library diversity info
  #remove reads containing indels
  noIndelsReads="$tmpDir"/${sampleID}"_alignmentFile2.sam"
  ../AssemblyPostProcessing/SamRemoveIndels.awk $assmbFile > $noIndelsReads
  #remove reads lacking unique fragment coordinates
  unqReads="$tmpDir"/${sampleID}"_alignmentFile3.sam"
  ../AssemblyPostProcessing/UniqSamPE.awk $noIndelsReads > $unqReads

  unqBamFile="$tmpDir"/${sampleID}"_alignmentFile_unqReads.bam"
  samtools view -bS $unqReads | samtools sort -o $unqBamFile &> /dev/null
  samtools index $unqBamFile &> /dev/null

  #get average depth of unqReads
  unqReadStats="$tmpDir"/${sampleID}"_unqStats.txt"
  weeSAMv1.4 -b $unqBamFile -out $unqReadStats  &> /dev/null
  
  unqAvrgDepth=$(awk 'NR==2{print $8;exit}' $unqReadStats)
  unqAvrgDepthRd=$(echo $unqAvrgDepth |awk '{printf "%.0f", $unqAvrgDepth}') 
  echo "Unique Reads Average Depth =  $unqAvrgDepthRd"  
    
  ratio=$(bc -l <<< "scale = 2; (($avrgDepthRd/$unqAvrgDepthRd))")
  echo "Ratio Coverage (total/unique)= $ratio "
 

  #get the number of strains in the sample
  # sigFile=${inputDir}/${sigName}
  strainOutput=${tmpDir}/${sampleID}_output.sig
  TMPFILE=$(mktemp ./tmp.XXXXXXXX)
  awk 'NR%4==2' $catFile > ${TMPFILE}
  echo "../GenotypingTools/MIRNA_SEARCH ${TMPFILE} ${sigName} > ${strainOutput}"
  ../GenotypingTools/MIRNA_SEARCH ${TMPFILE} ${sigName} > ${strainOutput}
  rm ${TMPFILE}
  

  #get the maximum number of strains in the sample
  numStrain=$(perl ./HCMV_GenotypeCount.pl -in ${strainOutput})
  echo "Number of strain = $numStrain "

  #Step5: print the result to the output file
  #An array contains the alignment statistical information
  refInfo=($sampleID $totalReads $mappedReads $prcMappedReads $avrgDepthRd $refLen $unqAvrgDepthRd $ratio $minDepth $maxDepth $stdDev $coverage1 $coverage2 $coverage3  $varCoef $numStrain)  
  SampleStats=${tmpDir}/${sampleID}_stats.txt
#  echo "\n$sampleID\n$totalReads\n$mappedReads\n$prcMappedReads\n$avrgDepthRd\n$refLen\n$unqAvrgDepthRd\n$ratio\n$minDepth\n$maxDepth\n$stdDev\n$coverage1\n$coverage2\n$coverage3\n$varCoef\n$numStrain\n" > $SampleStats
  echo "" > $SampleStats
  printf '%s\n' $sampleID $totalReads $mappedReads $prcMappedReads $avrgDepthRd $refLen $unqAvrgDepthRd $ratio $minDepth $maxDepth $stdDev $coverage1 $coverage2 $coverage3  $varCoef $numStrain >> $SampleStats
  echo -e "\n" >> $SampleStats
  awk '{print $1}' $strainOutput >> $SampleStats

  # pasting the results to each other
  if [[ $k == 1 ]]; then
    paste -d "," $labelFile  $SampleStats > $summaryFile     
  else
    # echo "No labels added"
    TMPFILE=$(mktemp ./tmp.XXXXXXXX)
    paste -d "," ${summaryFile} ${SampleStats} > ${TMPFILE}
    mv ${TMPFILE} ${summaryFile}
#    rm ${TMPFILE}
  fi



  
  
  echo "##################################################"
  ((k++))  
  done

mv $summaryFile ${inputDir}/output.csv  
#remove tmpDir
rm -rf $tmpDir 
echo "Done!"
  
