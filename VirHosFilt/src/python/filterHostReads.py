import sys
import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

inputFolder = sys.argv[1]
hostReference = sys.argv[2]
virusReference = sys.argv[3]
outputFolder = sys.argv[4]




onlyfiles = [f for f in listdir(inputFolder) if isfile(join(inputFolder, f))]

fileNumber = {}
for item in onlyfiles:
    if not item[:-8] in fileNumber:
        fileNumber[item[:-8]] = 1
    else:
        fileNumber[item[:-8]] += 1

paired2filter = []
single2filter = []
for files in fileNumber:
    if fileNumber[files]==2:
        paired2filter.append(files)
    if fileNumber[files]==1:
        single2filter.append(files)

#Print datasets list
#outfile = open("pairedEnd_datasets.txt","w")
#for dataset in paired2filter:
#    outfile.write("-1 "+inputFolder+dataset+"_1.fastq -2 "+inputFolder+dataset+"_2.fastq")
#outfile.close()

#outfile = open("singleEnd_datasets.txt","w")
#for dataset in single2filter:
#    outfile.write("-U "+inputFolder+dataset+"_1.fastq\n")
#outfile.close()

#print "\n"
#print len(paired2filter),"paired end dataset(s) and",len(single2filter),"dataset(s) were detected"
#print "Press a key to close this window and start the filtering process...."
#sys.stdin.read(1)


for dataset in paired2filter:
    if not dataset[0] ==".":
        print "Filtering dataset",dataset
        print "Mapping reads",dataset,"to the host reference genome...."
        os.system("bowtie2 -x "+hostReference+" -1 "+inputFolder+dataset+"_1.fastq -2 "+inputFolder+dataset+"_2.fastq -p 6 -S alignment.sam")
        print "Converting alignment format for reads",dataset+"...."
        os.system("samtools view -bS -h alignment.sam >alignment.bam") 
        print "Extracting unmapped reads for dataset",dataset
        os.system("bam2fastq --no-aligned --force --strict -o unmapped#.fq alignment.bam")
        os.system("mv unmapped_1.fq "+outputFolder+"/"+dataset+"_noHost_1.fastq")
        os.system("mv unmapped_2.fq "+outputFolder+"/"+dataset+"_noHost_2.fastq")



for dataset in single2filter:
    if not dataset[0] ==".":
        print "Filtering dataset",dataset
        print "Mapping reads",dataset,"to the host reference genome...."
        print "bowtie2 -x "+hostReference+" -U "+inputFolder+dataset+"_1.fastq -p 6 -S alignment.sam"

        os.system("bowtie2 -x "+hostReference+" -U "+inputFolder+dataset+"_1.fastq -p 6 -S alignment.sam")
        print "Converting alignment format for reads",dataset+"...."
        os.system("samtools view -bS -h alignment.sam >alignment.bam") 
        print "Extracting unmapped reads for dataset",dataset
        os.system("bam2fastq --no-aligned --force --strict -o unmapped#.fq alignment.bam")
        os.system("mv unmapped_M.fq "+outputFolder+"/"+dataset+"_noHost_1.fastq")


if not virusReference == "None":
    print "Mapping host free reads to the virus reference genome"
    onlyfiles = [f for f in listdir(outputFolder) if isfile(join(outputFolder, f))]
    print "Found ",len(onlyfiles),"in outputFolder"
    fileNumber = {}
    for item in onlyfiles:
        if not item[:-8] in fileNumber:
            fileNumber[item[:-8]] = 1
        else:
            fileNumber[item[:-8]] += 1

    paired2filter = []
    single2filter = []

    print "Found",len(paired2filter),"paired datasets"
    for files in fileNumber:
        if fileNumber[files]==2:
            paired2filter.append(files)
        if fileNumber[files]==1:
            single2filter.append(files)

    for dataset in paired2filter:
        if not dataset[0] ==".":
            print "Mapping reads",dataset,"to the host reference genome...."
            os.system("bowtie2 -x "+virusReference+" -1 "+outputFolder+"/"+dataset+"_1.fastq -2 "+outputFolder+"/"+dataset+"_2.fastq -p 6 -S alignment.sam")
            print "Converting alignment format for reads",dataset+"...."
            os.system("samtools view -bS -h alignment.sam >alignment.bam") 
            print "Sorting bam file for reads",dataset+"...."
            os.system("samtools sort -o alignment_sorted.bam alignment.bam")
            print "Calculating virus coverage for dataset",dataset
            os.system("samtools depth  alignment_sorted.bam >coverage.txt")
            covFile = open("coverage.txt")
            position = []
            coverage = []
            while True:
                line = covFile.readline().rstrip()
                if not line:
                    break
                fields = line.split("\t")
                position.append(int(fields[1]))
                coverage.append(int(fields[2]))
            plt.plot(position,coverage)
            plt.savefig(dataset+"_covPlot.png")
            plt.close()
            covFile.close()
            os.system("mv *_covPlot.png "+outputFolder+"/")
            os.system("mv alignment_sorted.bam "+outputFolder+"/virusAlignment.bam")
            os.system("rm *.sam *.bam")



    for dataset in single2filter:
        if not dataset[0] ==".":
            print "Mapping reads",dataset,"to the host reference genome...."
            os.system("bowtie2 -x "+virusReference+" -U "+outputFolder+"/"+dataset+"_1.fastq -p 6 -S alignment.sam")
            print "Converting alignment format for reads",dataset+"...."
            os.system("samtools view -bS -h alignment.sam >alignment.bam") 
            print "Sorting bam file for reads",dataset+"...."
            os.system("samtools sort -o alignment_sorted.bam alignment.bam")
            print "Calculating virus coverage for dataset",dataset
            os.system("samtools depth  alignment_sorted.bam >coverage.txt")
            covFile = open("coverage.txt")
            position = []
            coverage = []
            while True:
                line = covFile.readline().rstrip()
                if not line:
                    break
                fields = line.split("\t")
                position.append(int(fields[1]))
                coverage.append(int(fields[2]))
            plt.plot(position,coverage)
            plt.savefig(dataset+"_covPlot.png")
            covFile.close()
            os.system("mv *_covPlot.png "+outputFolder+"/")
            os.system("mv alignment_sorted.bam "+outputFolder+"/virusAlignment.bam")
            os.system("rm *.sam *.bam")



    



    





