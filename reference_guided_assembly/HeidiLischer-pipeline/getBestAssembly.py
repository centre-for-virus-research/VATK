import sys
import os
from Bio import SeqIO

read1 = sys.argv[1]
read2 = sys.argv[2]

outfile = open(sys.argv[3],"w")
outfile.write("RunNumber\tSubSampledReads\tN50\tLongestScaffold\n")

#get total number of reads
print "Retriving number of reads...."
os.system("wc -l "+read1+" >numReads.txt")
nr = open("numReads.txt")
numReads = float(((nr.readline().rstrip()).split(" "))[0])/4
nr.close()
print "The dataset contains ",numReads,"reads"


bestAssemblyN50 = ""
bestN50 = 0
bestLength = 0
bestAssembly_length = 0
for b in range(2):
    for a in range(int(numReads*0.1),int(numReads*0.9),int(numReads*0.2)):
        print "Subsampling ",a,"reads...."
        os.system("seqtk sample -s"+str(a*b)+" "+read1+" "+str(a)+" >subsample_1.fq")
        os.system("seqtk sample -s"+str(a*b)+" "+read2+" "+str(a)+" >subsample_2.fq")
        print "Performing de novo assembly...."
        os.system("/home3/scc20x/Software/SPAdes-3.12.0-Linux/bin/spades.py -1 subsample_1.fq -2 subsample_2.fq --cov-cutoff auto --careful -k 53,63,73,83 -o outputSpades_"+str(a)+"_"+str(b)+" >null 2>&1")
        longestScaffold = 0
        if os.path.isfile("./outputSpades_"+str(a)+"_"+str(b)+"/scaffolds.fasta")==True:
            scaffoldsFile = open("./outputSpades_"+str(a)+"_"+str(b)+"/scaffolds.fasta")
            line = scaffoldsFile.readline()
            line = scaffoldsFile.readline()
            longestScaffold += len(line)
            while not line[0] == ">":
                line = scaffoldsFile.readline()
                if not line:
                    break
                longestScaffold += len(line)

            if longestScaffold > bestLength:
                bestLength = longestScaffold
                bestAssembly_length = "./outputSpades_"+str(a)+"_"+str(b)+"/scaffolds.fasta"

            os.system("getN50 ./outputSpades_"+str(a)+"_"+str(b)+"/scaffolds.fasta >N50.txt")
            os.system("fold -1 outputSpades_"+str(a)+"_"+str(b)+"/scaffolds.fasta | grep -c N >> N50.txt") 
            infile = open("N50.txt")
            infile.readline()
            infile.readline()
            n50 = int(((infile.readline().rstrip()).split(" "))[1])
            numN = infile.readline().rstrip()
            print a,n50,longestScaffold,numN
            outfile.write(str(b)+"\t"+str(a)+"\t"+str(n50)+"\t"+str(longestScaffold)+"\t"+str(numN)+"\n")
            if n50 > bestN50:
                bestN50 = n50
                bestAssemblyN50 = "./outputSpades_"+str(a)+"_"+str(b)+"/scaffolds.fasta"

        

print "The best assembly for N50 is in ",bestAssemblyN50
print "The best assembly for longest scaffold is in ",bestAssembly_length
os.system("cp "+bestAssemblyN50+" .")
os.system("rm -rf outputSpades*")

