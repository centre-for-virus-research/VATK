#!/bin/bash
#===================================================================
# Maha Maabar, 27/01/2017
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

#INPUT PARAMETERS 
inputDir="$1"
refName="$2"
sigName="$3"
libName="$4"

#create the output file with the headings info
outputFile="$inputDir"/output.csv
formatOutput.sh $inputDir/$sigName $outputFile



#create a temp directory in current directory to hold intermediate files
tmpDir="./tmpDir"
mkdir tmpDir
chmod 775 tmpDir


#Prepare library for alignment
refFile="$inputDir"/"$refName"

bowtie2-build "$refFile" "$libName" &> /dev/null
mv $libName.* $tmpDir

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
  trim_galore --paired --length 21 --quality 10 --stringency 3 "$file" "$file2" &> /dev/null

  trimFile1=${fileName1%.*}"_val_1.fq"
  trimFile2=${fileName2%.*}"_val_2.fq"
  
  mv *.fq $tmpDir 
  rm *_report.txt
  
  catFile="$tmpDir"/"$sampleID""_catfile.fq" #concatenate trimmed files
  cat "$tmpDir"/"$trimFile1" "$tmpDir"/"$trimFile2" > "$catFile"
  echo "Trimmed and concatenated file: $catFile"

  #get the total number of reads 
  numReads=$(wc -l "$catFile" | cut -d' ' -f1) 
  totalReads=$((numReads / 4))
  echo "Total Number of Reads $totalReads" 

  #Step2: assemble trimmed files to Merlin
  #bowtie2 assembly
  echo "Assembling to Merlin ..."
  assmbFile="$tmpDir"/${sampleID}"_alignmentFile.sam"
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
  refLen=$(echo $statsLine |cut -d ' ' -f2)
  mappedReads=$(echo $statsLine |cut -d ' ' -f3)
  prcMappedReads=$(bc -l <<< "scale = 2; (($mappedReads/$totalReads*100))")
  minDepth=$(echo $statsLine |cut -d ' ' -f6)
  maxDepth=$(echo $statsLine |cut -d ' ' -f7)
  avrgDepth=$(echo $statsLine |cut -d ' ' -f8)
  avrgDepthRd=$(echo $avrgDepth |awk '{printf "%.0f", $avrgDepth}')
  echo "Average Depth =  $avrgDepthRd" 

  stdDev=$(echo $statsLine |cut -d ' ' -f9)
  
  coverage1=$(echo $statsLine |cut -d ' ' -f10) #proportion of num of sites with depth > 0.2*AverageDepth
  coverage2=$(echo $statsLine |cut -d ' ' -f12) #proportion of num of sites with depth > 1.0*AverageDepth
  coverage3=$(echo $statsLine |cut -d ' ' -f13) #proportion of num of sites with depth > 1.8*AverageDepth
  
  varCoef=$(echo $statsLine |cut -d ' ' -f14)  #coefficient of variation

  #Step4: get library diversity info
  #remove reads containing indels
  noIndelsReads="$tmpDir"/${sampleID}"_alignmentFile2.sam"
  SamRemoveIndels.awk $assmbFile > $noIndelsReads
  #remove reads lacking unique fragment coordinates
  unqReads="$tmpDir"/${sampleID}"_alignmentFile3.sam"
  UniqSamPE.awk $noIndelsReads > $unqReads

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
  sigFile="$inputDir"/"$sigName"
  strainOutput="$tmpDir"/${sampleID}"_output.sig"
  ReadSearchAndrew.sh $catFile $sigFile > $strainOutput
  genotypeFile="$inputDir"/"genotypes.txt"

  #get the maximum number of strains in the sample
  #javaPath= "java -cp /software/bin Genotype"
  #numStrain=$($javaPath $genotypeFile $strainOutput )
  numStrain=$(java -cp /software/bin Genotype $genotypeFile $strainOutput )
  echo "Number of strain = $numStrain "

  #An array contains the alignment statistical information
  refInfo=($sampleID $totalReads $mappedReads $prcMappedReads $avrgDepthRd $refLen $unqAvrgDepthRd $ratio $minDepth $maxDepth $stdDev $coverage1 $coverage2 $coverage3  $varCoef $numStrain)  

  #An array contains the strain information 
  result=$(cut -f1 $strainOutput )
  read -a strainInfo  <<<$result
      
  #Step5: print the result to the output file
  tmpFile="$tmpDir"/"tmp.csv"
  d=1
  while read line
  do 
    if [[ $d == 1 ]] ; then echo $line >> $tmpFile 
    elif [[ $d == 2 ]]; then 
         for (( i=0; i< ${#refInfo[@]} ; i++)) ; do
             echo $line$','${refInfo[i]} >> $tmpFile 
             read line
        done
        echo $line >> $tmpFile #print the last read line
        d=$i
    elif [[ $d == 17 ]]; then echo $line$','Reads' (no.)' >> $tmpFile # 17= num of refInfo + 1
    else
        for (( i=0; i< ${#strainInfo[@]} ; i++)) ; do
           echo $line$','${strainInfo[i]} >> $tmpFile 
           read line
        done
        d=$i
    fi
    ((d++))
  done <"$outputFile" 

  mv $tmpFile $outputFile   
  
  echo "##################################################"
  ((k++))  
  done


  #remove tmpDir
  rm -rf $tmpDir 
  echo "Done!"
  
