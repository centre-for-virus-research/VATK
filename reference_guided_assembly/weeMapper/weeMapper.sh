#!/bin/bash

################################------- This is a Mapping Script --------########################################
# This script is written by Sejal Modha										#
#														#
# This script can be used to map the paired endinput reads to following human DBs				#	
#	~ Human Genome												#
#	~ Human cDNA												#
#	~ Human CDS												#
#	~ Human RepBase												#
#	~ Human HLA												#
#	~ Human Promoters											#
#---------------------------------------------------------------------------------------------------------------#
# Usage:													#
#	./weeMapper <OPTIONS>											#
#	 This script uses bowtie2 for mapping 									#
#														#
# 	-h This helpful message 										#
# 	-1 Input File1.fq 											#
# 	-2 Input File2.fq 											#
# 	-a Clean the input sequences against all human databases(genome, cDNA, CDS, RepBase, HLA, Promoters) 	#
# 	-b Clean against human Genome 										#
# 	-c Clean against human cDNA 										#
# 	-d Clean against human CDS 										#
# 	-e Clean against human RepBase 										#
# 	-f Clean against human HLA										#
#	-g Clean against human Promoters 									#
#														#
#---------------------------------------------------------------------------------------------------------------#

usage=`echo -e "\n Usage: weeMapper <OPTIONS> \n This script uses bowtie2 for mapping \n\n -h This helpful message \n -1 Input File1.fq \n -2 Input File2.fq \n -a Clean the input sequences against all human databases(genome, cDNA, CDS, RepBase, HLA, Promoters) \n -b Clean against human genome \n -c Clean against human cDNA \n -d Clean against human CDS \n -e Clean against human RepBase \n -f Clean against human HLA \n -g Clean against human Promoters"`

if [[ ! $1 ]] 
then
	printf "${usage}\n\n";
exit;
fi

while getopts 1:2:abcdefgh flag; do
  case $flag in

    1)
	file1=`echo "$OPTARG"`;
	f1=${file1%.fq};
	echo $f1;
	;;
    2)
	file2=`echo "$OPTARG"`;
	f2=${file2%.fq};
	echo $f2;
	#echo -e "Input file1 is $file1"; 
	#echo -e "Input file2 is $file2"; 

	#rm -f ${f1}_${f2}_cleaning_log;
	#unlink temp1;
	#unlink temp2;
	ln -s $file1 temp1;
	ln -s $file2 temp2;
	
	input1=`echo temp1`;
	input2=`echo temp2`;
	
	#rm -f ${f1}_${f2}_cleaning_log;
	touch ${f1}_${f2}_cleaning_log;
	num1=`grep "^@" $file1|wc -l`;
	num2=`grep "^@" $file2|wc -l`;
	reads=`expr $num1 + $num2`;
	echo "--------------Total Number of Input Reads----------------" >>  ${f1}_${f2}_cleaning_log;
	echo -e "Total Number of Reads \n$reads " >> ${f1}_${f2}_cleaning_log;
	;;
    a)
	#Align against all DBs
##################################### mapping to Human Genome ##################################
	echo "Aligning against human genome hg38 using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_genome.sam -x /db/wgs/Human/hg38/bowtie2 -1 $input1 -2 $input2 -p 6;
	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_genome.sam -b -o ${input1}_${input2}_human_genome.bam -@ 6;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_genome.bam;

	
	echo "------------------Human Genome DB----------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_genome.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_genome.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_genome.bam >> ${f1}_${f2}_cleaning_log;
	
################################### Mapping to Human cDNA #######################################
	echo "Aligning against human cDNA using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_cdna.sam -x /db/wgs/Human/hg38/bowtie2/Homo_sapiens_cdna -1 ${input1}_${input2}_unmapped_1.fq -2 ${input1}_${input2}_unmapped_2.fq;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_cdna.sam -b -o ${input1}_${input2}_human_cdna.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_cdna.bam;

	echo "------------------Human cDNA DB------------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_cdna.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_cdna.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_cdna.bam >> ${f1}_${f2}_cleaning_log;
	
##################################### mapping to Human CDS ##################################
	echo "Aligning against human CDS using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_cds.sam -x /db/wgs/Human/hg38/bowtie2/Homo_sapiens_cds -1 ${input1}_${input2}_unmapped_1.fq -2 ${input1}_${input2}_unmapped_2.fq;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_cds.sam -b -o ${input1}_${input2}_human_cds.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_cds.bam;

	echo "------------------Human CDS DB------------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_cds.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_cds.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_cds.bam >> ${f1}_${f2}_cleaning_log;
	
##################################### mapping to Human RepBase ################################
	echo "Aligning against human RepBase using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_repbase.sam -x /db/wgs/Human/hg38/bowtie2/human_repeats_repbase -1 ${input1}_${input2}_unmapped_1.fq -2 ${input1}_${input2}_unmapped_2.fq;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_repbase.sam -b -o ${input1}_${input2}_human_repbase.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_repbase.bam;

	echo "------------------Human RepBase DB---------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_repbase.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_repbase.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_repbase.bam >> ${f1}_${f2}_cleaning_log;
	
##################################### mapping to Human HLA ##################################
	echo "Aligning against human HLA using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_hla.sam -x /db/wgs/Human/hg38/bowtie2/human_hla -1 ${input1}_${input2}_unmapped_1.fq -2 ${input1}_${input2}_unmapped_2.fq;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_hla.sam -b -o ${input1}_${input2}_human_hla.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_hla.bam;

	echo "------------------Human HLA DB------------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_hla.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_hla.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_hla.bam >> ${f1}_${f2}_cleaning_log;


################################## mapping to Human Promoters sequences #######################	
	echo "Aligning against human Promoters sequences using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_promoters.sam -x /db/wgs/Human/hg38/bowtie2/EPD_human -1 ${input1}_${input2}_unmapped_1.fq -2 ${input1}_${input2}_unmapped_2.fq;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_promoters.sam -b -o ${input1}_${input2}_human_promoters.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_promoters.bam;

	echo "------------------Human Promoters DB--------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_promoters.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_promoters.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_promoters.bam >> ${f1}_${f2}_cleaning_log;
	
	unlink temp1;
	unlink temp2;

	ln -s ${input1}_${input2}_unmapped_1.fq temp1;
	ln -s ${input1}_${input2}_unmapped_2.fq temp2;
	;;

    b)
	#Align against human genome sequence only
	
	
	#echo -e "$input1 \n $input2";

	echo "Aligning against human genome hg38 using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_genome.sam -x /db/wgs/Human/hg38/bowtie2/hg38 -1 $input1 -2 $input2 -p 6;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_genome.sam -b -o ${input1}_${input2}_human_genome.bam -@ 6;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_genome.bam;

	echo "------------------Human Genome DB----------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_genome.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_genome.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	#echo -e "$mp %">> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_genome.bam >> ${f1}_${f2}_cleaning_log;
	
	unlink temp1;
	unlink temp2;

	ln -s ${input1}_${input2}_unmapped_1.fq temp1;
	ln -s ${input1}_${input2}_unmapped_2.fq temp2;

	;;	
	
    c)
	#Align Against human cDNA sequences
	
	input1=`echo temp1`;
	input2=`echo temp2`;
	
	#echo -e "$input1 \n $input2";
	
	echo "Aligning against human cDNA using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_cdna.sam -x /db/wgs/Human/hg38/bowtie2/Homo_sapiens_cdna -1 $input1 -2 $input2;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_cdna.sam -b -o ${input1}_${input2}_human_cdna.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_cdna.bam;

	echo "------------------Human cDNA DB------------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_cdna.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_cdna.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo -e "$mr ($mp%)">> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_cdna.bam >> ${f1}_${f2}_cleaning_log;
	
	unlink temp1;
	unlink temp2;

	ln -s ${input1}_${input2}_unmapped_1.fq temp1;
	ln -s ${input1}_${input2}_unmapped_2.fq temp2;
	
	;;
	
    d)
	#Align against Human CDS
	input1=`echo temp1`;
	input2=`echo temp2`;
	
	#echo -e "$input1 \n $input2";
	
	echo "Aligning against human CDS using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_cds.sam -x /db/wgs/Human/hg38/bowtie2/Homo_sapiens_cds -1 $input1 -2 $input2;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_cds.sam -b -o ${input1}_${input2}_human_cds.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_cds.bam;

	echo "------------------Human CDS DB------------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_cds.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_cds.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_cds.bam >> ${f1}_${f2}_cleaning_log;
	
	unlink temp1;
	unlink temp2;

	ln -s ${input1}_${input2}_unmapped_1.fq temp1;
	ln -s ${input1}_${input2}_unmapped_2.fq temp2;
	
	;;
	
    e)
	#Align against human RepBase
	
	input1=`echo temp1`;
	input2=`echo temp2`;
	
	#echo -e "$input1 \n $input2";
	
	echo "Aligning against human RepBase using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_repbase.sam -x /db/wgs/Human/hg38/bowtie2/human_repeats_repbase -1 $input1 -2 $input2;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_repbase.sam -b -o ${input1}_${input2}_human_repbase.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_repbase.bam;

	echo "------------------Human RepBase DB---------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_repbase.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_repbase.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_repbase.bam >> ${f1}_${f2}_cleaning_log;
	
	unlink temp1;
	unlink temp2;

	ln -s ${input1}_${input2}_unmapped_1.fq temp1;
	ln -s ${input1}_${input2}_unmapped_2.fq temp2;
	
	;;

    f)
	#Align against human RepBase
	
	input1=`echo temp1`;
	input2=`echo temp2`;
	
	#echo -e "$input1 \n $input2";
	
	echo "Aligning against human HLA using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_hla.sam -x /db/wgs/Human/hg38/bowtie2/human_hla -1 $input1 -2 $input2;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_hla.sam -b -o ${input1}_${input2}_human_hla.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_hla.bam;

	echo "------------------Human HLA DB------------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_hla.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_hla.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_hla.bam >> ${f1}_${f2}_cleaning_log;
	
	unlink temp1;
	unlink temp2;

	ln -s ${input1}_${input2}_unmapped_1.fq temp1;
	ln -s ${input1}_${input2}_unmapped_2.fq temp2;

	;;	
	
    g)
	#Align against human Promoters sequences
	
	input1=`echo temp1`;
	input2=`echo temp2`;
	
	#echo -e "$input1 \n $input2";
	
	echo "Aligning against human Promoters sequences using bowtie2";
	bowtie2 -q  -S ${input1}_${input2}_human_promoters.sam -x /db/wgs/Human/hg38/bowtie2/EPD_human -1 $input1 -2 $input2;

	echo "Converting SAM to BAM";
	samtools view -S ${input1}_${input2}_human_promoters.sam -b -o ${input1}_${input2}_human_promoters.bam;

	echo "Extracting Unmapped Sequences";
	bam2fastq --no-aligned --force --strict -o ${input1}_${input2}_unmapped#.fq ${input1}_${input2}_human_promoters.bam;

	echo "------------------Human Promoters DB--------------------" >> ${f1}_${f2}_cleaning_log;
	echo "Number of mapped reads" >> ${f1}_${f2}_cleaning_log;
	#samtools view -F 4 -c ${input1}_${input2}_human_promoters.bam >> ${f1}_${f2}_cleaning_log;
	mr=`samtools view -F 4 -c ${input1}_${input2}_human_promoters.bam`;
	mp=`expr $mr \* 100 / $reads`;
	echo "$mr ($mp%)" >> ${f1}_${f2}_cleaning_log;
	echo "Number of unmapped reads" >> ${f1}_${f2}_cleaning_log;
	samtools view -f 4 -c ${input1}_${input2}_human_promoters.bam >> ${f1}_${f2}_cleaning_log;
	
	unlink temp1;
	unlink temp2;

	ln -s ${input1}_${input2}_unmapped_1.fq temp1;
	ln -s ${input1}_${input2}_unmapped_2.fq temp2;
	;;
    h)
     	printf "${usage}\n\n";
	#echo "-a Bowtie help ";
	#echo "-h This helpful screen";
	exit;
      ;;

    \?)
	echo -e "\n Option you selected doesn't exist";
	printf "${usage}\n\n";
	exit;
      ;;
	
  esac
done

echo "Removing temporary files";
unlink temp1;
unlink temp2;
rm -f temp1_temp2_*.sam;

echo "Renaming output files"
mv temp1_temp2_unmapped_1.fq ${f1}_${f2}_unmapped_1.fq;
mv temp1_temp2_unmapped_2.fq ${f1}_${f2}_unmapped_2.fq;

if [ -e "temp1_temp2_human_genome.bam" ]
then
	mv temp1_temp2_human_genome.bam ${f1}_${f2}_human_genome.bam;
fi

if [ -e "temp1_temp2_human_cdna.bam" ]
then
	mv temp1_temp2_human_cdna.bam ${f1}_${f2}_human_cdna.bam;
fi

if [ -e "temp1_temp2_human_cds.bam" ]
then
	mv temp1_temp2_human_cds.bam ${f1}_${f2}_human_cds.bam;
fi

if [ -e "temp1_temp2_human_repbase.bam" ]
then
	mv temp1_temp2_human_repbase.bam ${f1}_${f2}_human_repbase.bam;
fi

if [ -e "temp1_temp2_human_hla.bam" ]
then
	mv temp1_temp2_human_hla.bam ${f1}_${f2}_human_hla.bam;
fi

if [ -e "temp1_temp2_human_promoters.bam" ]
then
	mv temp1_temp2_human_promoters.bam ${f1}_${f2}_human_promoters.bam;
fi

echo "Printing log";
echo -e "\n****************Cleaning Log*********************\n";
cat ${f1}_${f2}_cleaning_log;
echo -e "\n";


