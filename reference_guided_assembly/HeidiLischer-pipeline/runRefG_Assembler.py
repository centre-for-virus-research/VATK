import os
import sys
import datetime
from Bio import SeqIO



confFile = open("config.txt")
logFile = open("logFile.log","w",buffering=0)

amosFolder = "/home3/scc20x/Software/bin/"


#Reading configuration file
#Read data
projectName = ((confFile.readline().rstrip()).split("\t"))[1]

read1 = ((confFile.readline().rstrip()).split("\t"))[1]
read2 = ((confFile.readline().rstrip()).split("\t"))[1]
reads2Sample = ((confFile.readline().rstrip()).split("\t"))[1]
referenceSeq = ((confFile.readline().rstrip()).split("\t"))[1]

now = datetime.datetime.now()
logFile.write("Date: "+now.strftime("%Y-%m-%d")+"\n")
logFile.write("Sample name: "+projectName+"\n")
logFile.write("Read1 fastq: "+read1+"\n")
logFile.write("Read2 fastq: "+read2+"\n\n\n")

#***************************************************************
#**************** 0 Preliminary elaboration ********************
#***************************************************************
read1Seq = {}
read2Seq = {}

print "Collecting reads from file",read1
for seq_record in SeqIO.parse(read1,"fastq"):
    if not str(seq_record.id) in read1Seq:
        read1Seq[str(seq_record.id)] = seq_record

print "Collecting reads from file",read2
for seq_record in SeqIO.parse(read2,"fastq"):
    if not str(seq_record.id) in read2Seq:
        read2Seq[str(seq_record.id)] = seq_record





#***************************************************************
#**************** 1 Reads quality filtering ********************
#***************************************************************

confFile.readline() #Read comment
qualityFiltering = ((confFile.readline().rstrip()).split("\t"))[1]
logFile.write("Quality filtering: "+qualityFiltering+"\n")

if qualityFiltering == "yes" or qualityFiltering == "Yes":
    now = datetime.datetime.now()
    logFile.write("Quality filtering started at "+now.strftime("%H:%M")+"\n")

    qualConfFile = open("qualityFiltering.conf","w")
    qualConfFile.write("ProjectName\tallSamples\n")
    qualConfFile.write("Sample_start****************************************"+"\n")
    qualConfFile.write("SampleName\t"+projectName+"\n")
    qualConfFile.write("Read1\t"+read1+"\n")
    qualConfFile.write("Read2\t"+read2+"\n")
    qualConfFile.write("Sample_start***Filtering and trimming options\n")
    for a in range(7):
        qualConfFile.write(confFile.readline())
    qualConfFile.write("SampleEnd*******************************************\n")
    qualConfFile.close()
    os.system("cp /home3/scc20x/Software/mySoftware/VATK/reference_guided_assembly/HeidiLischer-pipeline/* .")
    os.system("python readsFiltering.py qualityFiltering.conf")
    os.system("mkdir -p 1_cleanReads")
    os.system("mv ./allSamples/"+projectName+"_hq_?.fastq ./1_cleanReads/")
    os.system("rm -rf allSamples")

    now = datetime.datetime.now()
    logFile.write("Quality filtering ended at "+now.strftime("%H:%M")+"\n\n")

else:
    for a in range(7):
        confFile.readline()
    if os.path.isfile("./1_cleanReads/"+projectName+"_hq_1.fastq") == False:
        print "File ","./1_cleanReads/",projectName,"_hq_1.fastq does not exist, now exiting...."
        exit()
    if os.path.isfile("./1_cleanReads/"+projectName+"_hq_2.fastq") == False:
        print "File ","./1_cleanReads/",projectName,"_hq_2.fastq does not exist, now exiting...."
        exit()

#Perform subsampling
print "\nPerforming reads fastq subsampling......"
os.system("python get_randomPE.py ./1_cleanReads/"+projectName+"_hq_1.fastq ./1_cleanReads/"+projectName+"_hq_2.fastq 30000")
os.system("mv ./1_cleanReads/"+projectName+"_hq_1.fastq.subset ./1_cleanReads/"+projectName+"_hq_1_subSample.fastq")
os.system("mv ./1_cleanReads/"+projectName+"_hq_2.fastq.subset ./1_cleanReads/"+projectName+"_hq_2_subSample.fastq")


#***************************************************************
#**************** 2 Reads alignment to reference ***************
#***************************************************************
confFile.readline() #Read comment
performAlignment = ((confFile.readline().rstrip()).split("\t"))[1]
if performAlignment == 'Yes' or performAlignment == 'yes':
    os.system("mkdir 2_referenceAlignment")
    os.chdir("2_referenceAlignment")
    os.system("bowtie2-build "+referenceSeq+" reference")
    os.system("bowtie2 --fast-local -p 10 -x reference -1 ../1_cleanReads/"+projectName+"_hq_1.fastq -2 ../1_cleanReads/"+projectName+"_hq_2.fastq -S alignment.sam")
    os.system("samtools view -bS alignment.sam >alignment.bam")
    os.system("samtools sort -o alignment_sorted.bam alignment.bam")
    os.system("samtools view -b -F 4 alignment_sorted.bam >mapped.bam") #Get all the mapped
    os.system("samtools view -b -f 4 alignment_sorted.bam >unmapped.bam") #Get all the unmapped
    os.system("samtools view -b -f 9 unmapped.bam > noMateMapped.bam") # Among the unmapped get for those whose mate did not map, therefore neither the mates mapped
    os.system("samtools view -b -F 8 -f 64 unmapped.bam >mate1NotMapped.bam") # Among the unmapped get the mates for which mate1 did not map
    os.system("samtools view -b -F 8 -f 128 unmapped.bam >mate2NotMapped.bam") # Among the unmapped get the mates for which mate2 did not map
    os.system("samtools view -b -q 10 mapped.bam >mappedHQ.bam")
    os.chdir("..")
else:
    if os.path.isfile("./2_referenceAlignment/mappedHQ.bam")== False:
        print "File ","./2_referenceAlignment/",projectName,"_hq_1.fastq does not exist, now exiting...."
        exit()


#***************************************************************
#********** 3 Blocks and superblocks construction   ************
#***************************************************************
confFile.readline() #Read comment
reconstructBlocks = ((confFile.readline().rstrip()).split("\t"))[1]
if reconstructBlocks == 'yes' or reconstructBlocks == 'Yes':
    os.system("mkdir 3_blocksAndSuperblocks")
    os.chdir("3_blocksAndSuperblocks")
    os.system("samtools faidx "+referenceSeq)
    os.system("bedtools genomecov -ibam ../2_referenceAlignment/mappedHQ.bam  -bga > coverage.txt") #Produce a general coverage file
    #With the following coverage is calculated considering only the proper paired reads
    os.system("samtools view -b -f 2 ../2_referenceAlignment/mappedHQ.bam | samtools sort - -n | bedtools bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | bedtools genomecov -i - -bga -g "+referenceSeq+".fai > coverage_Paired.txt")
    os.system("java -jar ../GetBlocks.jar -i coverage.txt -paired coverage_Paired.txt -o blocks.txt -oSuper superBlocks.txt -mCov 10 -sLength 12000 -sOverlap 300 -maxLength 100000")
    os.chdir("..")

#***************************************************************
#************** 4 Reassemble superblocks de novo   *************
#***************************************************************
confFile.readline() #Read comment
superBlocksDenovo = ((confFile.readline().rstrip()).split("\t"))[1]
if superBlocksDenovo == 'yes' or superBlocksDenovo == 'Yes':
    os.system("mkdir 4_SuperblocksDeNovoAssembly")
    os.chdir("4_SuperblocksDeNovoAssembly")

    superBlocksFile = "../3_blocksAndSuperblocks/superBlocks.txt"
    infile = open(superBlocksFile)
    os.system("samtools index ../2_referenceAlignment/mapped.bam")
    while True:
        interval = infile.readline().rstrip()
        if not interval:
            break
        print "Extracting reads mapping in interval",interval,"...."
        os.system("samtools view  ../2_referenceAlignment/mapped.bam "+interval + " | cut -f 1 | sort -u >mappedReadsNames.txt")
        mappedReadsFile = open("mappedReadsNames.txt")
        temp1fq = open("temp_1.fq","w")
        temp2fq = open("temp_2.fq","w")
        print "Extracting reads"
        while True:
            line = mappedReadsFile.readline().rstrip()
            if not line:
                break
            SeqIO.write(read1Seq[line+"/1"],temp1fq,"fastq")
            SeqIO.write(read2Seq[line+"/2"],temp2fq,"fastq")
        temp1fq.close()
        temp2fq.close()
        mappedReadsFile.close()
        os.system("mv ../getBestAssembly.py ./")
        os.system("python getBestAssembly.py temp_1.fq temp_2.fq assemblyStats.txt")
        os.system("mkdir "+interval)
        os.system("mv scaffolds.fasta ./"+interval)
 


    temp1fq = open("temp_1.fq","w")
    temp2fq = open("temp_2.fq","w")
    tempsfq = open("temp_s.fq","w")
    os.system("samtools view ../2_referenceAlignment/unmapped.bam >unmapped.sam")
    samFile = open("unmapped.sam")
    unmappedReads = {}
    while True:
        line = samFile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if not fields[0] in unmappedReads:
            unmappedReads[fields[0]] = []
        unmappedReads[fields[0]].append(fields[9])
        unmappedReads[fields[0]].append(fields[10])


    for item in unmappedReads:
        if len(unmappedReads[item])>2:
            temp1fq.write("@"+item+"\n"+unmappedReads[item][0]+"\n+\n"+unmappedReads[item][1]+"\n")
            temp2fq.write("@"+item+"\n"+unmappedReads[item][2]+"\n+\n"+unmappedReads[item][3]+"\n")
        else:
            tempsfq.write("@"+item+"\n"+unmappedReads[item][0]+"\n+\n"+unmappedReads[item][1]+"\n")

    temp1fq.close()
    temp2fq.close()
    tempsfq.close()
    os.system("/home3/scc20x/Software/SPAdes-3.12.0-Linux/bin/spades.py -1 temp_1.fq -2  temp_2.fq -s temp_s.fq --cov-cutoff auto --careful -k 51,61,71 -o unmappedSpades")

    #Add the code to collect unmapped paired and unpaired reads from bam/sam file
    os.system("cat */scaffolds.fasta >mergedScaffolds.fasta") #Da modificare con il nome generico
    #if os.path.isfile("./unmappedSpades/scaffolds.fasta")==True:
    #    os.system("cat mergedScaffolds.fasta ./unmappedSpades/scaffolds.fasta >temp.fasta; mv temp.fasta  mergedScaffolds.fasta")
    os.chdir("..")

#***************************************************************
#******* 5 Remove scaffolds redundancy with AMOScmp   **********
#***************************************************************
confFile.readline() #Read comment
superBlocksDenovo = ((confFile.readline().rstrip()).split("\t"))[1]
if superBlocksDenovo == 'yes' or superBlocksDenovo == 'Yes':
    os.system("mkdir 5_RemoveRedundancy")
    os.chdir("5_RemoveRedundancy")
    #Remove scaffolds shorter than 500
    outfile = open("scaffolds.fasta","w")
    for seq_record in SeqIO.parse("../4_SuperblocksDeNovoAssembly/mergedScaffolds.fasta","fasta"):
        if len(str(seq_record.seq))>=500:
            SeqIO.write(seq_record,outfile,"fasta")
    outfile.close()
    os.system("java -jar ../FastaToAmos.jar -i scaffolds.fasta -o scaffolds_amos.fasta")
    os.system(amosFolder+"bank-transact -c -z -b scaffolds.bnk -m scaffolds_amos.fasta")
    os.system(amosFolder+"dumpreads scaffolds.bnk > scaffolds.seq")
    os.system("nucmer --maxmatch  --prefix=scaffolds "+referenceSeq+" scaffolds.seq")
    os.system("casm-layout -t 1000 -U scaffolds.layout -C scaffolds.conflict -b scaffolds.bnk scaffolds.delta")
    os.system("make-consensus -o 10 -B -b scaffolds.bnk")
    os.system("bank2contig scaffolds.bnk > scaffolds.contig")
    os.system("bank2fasta -b scaffolds.bnk > scaffolds.fasta")
    os.system("java -jar ../RemoveShortSeq.jar -i scaffolds.fasta -o scaffolds_unique.fasta -length 1 -u")
    os.system("java -jar ../FastaStats.jar -i  scaffolds_unique.fasta -min 200 >statistics.txt")
    os.chdir("..")

#***************************************************************
# 6 Read mapping on contigs and de novo assembly on unmapped  **
#***************************************************************
confFile.readline() #Read comment
readsRemapping = ((confFile.readline().rstrip()).split("\t"))[1]
if readsRemapping == 'yes' or readsRemapping == 'Yes':
    os.system("mkdir 6_ReadsRemapping")
    os.chdir("6_ReadsRemapping")
    os.system("mv ../5_RemoveRedundancy/scaffolds_unique.fasta .")
    os.system("bowtie2-build scaffolds_unique.fasta supercontigs")
    os.system("bowtie2 --sensitive -p 8 -q   -x supercontigs -1 ../1_cleanReads/"+projectName+"_hq_1.fastq -2 ../1_cleanReads/"+projectName+"_hq_2.fastq  | samtools  view -bS - | samtools  sort - -T shortName -o alignment_sorted.bam")
    os.system("samtools index alignment_sorted.bam")
    os.system("samtools view -b -f 4 alignment_sorted.bam >supercontUnmapped.bam")
    os.system("samtools view -b -f 9 supercontUnmapped.bam >supercontUnmapped_pairs.bam")
    os.system("java -jar /home3/scc20x/Software/picard-tools/picard-tools-1.109/SamToFastq.jar INPUT=supercontUnmapped_pairs.bam FASTQ=unmapped.1.fastq SECOND_END_FASTQ=unmapped.2.fastq")
    os.system("samtools view -b -F 8  supercontUnmapped.bam | bamtools convert -format fastq -out unmapped.s.fastq")
    os.system("bamtools stats -in supercontUnmapped.bam >unmappedBamStatistics.txt")
    os.system("samtools view -b -F 4 -q 10 alignment_sorted.bam >alignment_mapped_sorted_filtered.bam")
    os.system("bamtools stats -in alignment_mapped_sorted_filtered.bam >mappedBamStatistics.txt")
    os.system("/home3/scc20x/Software/SPAdes-3.12.0-Linux/bin/spades.py -1 unmapped.1.fastq -2  unmapped.2.fastq -s unmapped.s.fastq --cov-cutoff auto --careful -k 51,61,71 -o unmappedSpades")
    outfile = open("unmappedReadsContigs_filtered.fasta","w")
    for seq_record in SeqIO.parse("./unmappedSpades/scaffolds.fasta","fasta"):
        if len(str(seq_record.seq)) >= 200:
            SeqIO.write(seq_record,outfile,"fasta")
    outfile.close()
    os.chdir("..")

#***************************************************************
#  ******************* Read correction *************************
#***************************************************************
confFile.readline() #Read comment
readsCorrection = ((confFile.readline().rstrip()).split("\t"))[1]
if readsCorrection == 'yes' or readsCorrection == 'Yes':
    os.system("mkdir 7_readsCorrection")
    os.chdir("7_readsCorrection")
    os.system("cat ../6_ReadsRemapping/scaffolds_unique.fasta ../6_ReadsRemapping/unmappedReadsContigs_filtered.fasta >mergedSequences.fasta")
    os.system("bowtie2-build mergedSequences.fasta merged")
    os.system("bowtie2 --sensitive -p 8 -q   -x merged -1 ../1_cleanReads/"+projectName+"_hq_1.fastq -2 ../1_cleanReads/"+projectName+"_hq_2.fastq  | samtools  view -bS - | samtools  sort - -T shortName -o alignment_sorted.bam")
    os.system("samtools index alignment_sorted.bam")
    os.system("samtools view -b -F 4 -q 10 alignment_sorted.bam >alignment_sorted_mapped_filtered.bam")
    os.system("samtools index alignment_sorted_mapped_filtered.bam")
    os.system("samtools faidx mergedSequences.fasta")
    os.system("createConsensus mergedSequences.fasta alignment_sorted.bam")
    os.system("java -jar ../RemoveShortSeq.jar -i mergedSequences.fasta_con.fasta -o mergedSequences.fasta_con_noN.fasta -length 100 -n -fq")
    os.system("bowtie2-build mergedSequences.fasta_con.fasta consensus")
    os.system("bowtie2 --sensitive -p 8 -q   -x consensus -1 ../1_cleanReads/"+projectName+"_hq_1.fastq -2 ../1_cleanReads/"+projectName+"_hq_2.fastq  | samtools  view -bS - | samtools  sort - -T shortName -o alignment_sorted2.bam")
    os.system("samtools index alignment_sorted2.bam")
    os.system("samtools view -b -F 4 -q 10 alignment_sorted2.bam >alignment_sorted2_mapped_filtered.bam")
    os.system("bedtools genomecov -ibam alignment_sorted2_mapped_filtered.bam -bga> alignment_sorted2_mapped_filtered_Cov.txt")
    os.system("samtools faidx mergedSequences.fasta_con.fasta")
    os.system("samtools view -bf 0x2 alignment_sorted2_mapped_filtered.bam | samtools sort - -n | bedtools bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | bedtools genomecov -i - -bga -g mergedSequences.fasta_con.fasta.fai  > mergedSequences.fasta_con.fasta_filteredPairedCov.txt")
    os.system("java -jar ../SplitSeqLowCov.jar -i alignment_sorted2_mapped_filtered_Cov.txt -paired mergedSequences.fasta_con.fasta_filteredPairedCov.txt -o mergedSequences.fasta_con.fasta_filteredNotCov.txt -mCov 1 -fasta mergedSequences.fasta_con.fasta -fastaOut mergedSequences.fasta_con.fasta_splitFiltered.fa")
    os.chdir("..")
#***************************************************************
#  ********************* gap Closing ***************************
#***************************************************************
confFile.readline() #Read comment
gapClosing = ((confFile.readline().rstrip()).split("\t"))[1]
if gapClosing == 'yes' or gapClosing == 'Yes':
    os.system("mkdir 8_GapClosing")
    os.chdir("8_GapClosing")
    os.system("java -jar ../WriteSoapConfig.jar -insLength 500 -r1 ../1_cleanReads/"+projectName+"_hq_1.fastq -r2 ../1_cleanReads/"+projectName+"_hq_2.fastq -max 300 -ru 2 -rank -o soap.config")
    os.system("../finalFusion -D -c ../7_readsCorrection/mergedSequences.fasta_con.fasta_splitFiltered.fa -K 61 -g hcmv_61  -p 16")
    os.system("SOAPdenovo-127mer map -s soap.config -g hcmv_61 -p 16")
    os.system("SOAPdenovo-127mer scaff -g hcmv_61 -p 16 -F")
    os.system("java -jar ../RemoveShortSeq.jar -i hcmv_61.scafSeq -o hcmv_61.scafSeq.fa -length 200")
    os.system("java -jar ../RemoveShortSeq.jar -i hcmv_61.scafSeq.fa -o hcmv_61.scafSeq_500.fa -length 500")
    os.system("java -jar ../RemoveShortSeq.jar -i hcmv_61.scafSeq.fa -o hcmv_61.scafSeq_1000.fa -length 1000")















  


    
