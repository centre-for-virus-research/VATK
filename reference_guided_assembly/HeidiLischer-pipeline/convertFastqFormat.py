import sys

read1 = sys.argv[1]
read2 = sys.argv[2]

#converting reads1
print "Converting reads file 1 ...."
infile = open(read1)
outfile = open("reads_1.fq","w")
while True:
    line1 = infile.readline().rstrip()
    if not line1:
        break
    fields = line1.split(" ")
    outfile.write(fields[0]+"/1\n")
    for a in range(3):
        line1 = infile.readline().rstrip()
        outfile.write(line1+"\n")

infile.close()
outfile.close()

#converting reads2
print "Converting reads file 2 ...."
infile = open(read2)
outfile = open("reads_2.fq","w")
while True:
    line1 = infile.readline().rstrip()
    if not line1:
        break
    fields = line1.split(" ")
    outfile.write(fields[0]+"/2\n")
    for a in range(3):
        line1 = infile.readline().rstrip()
        outfile.write(line1+"\n")

infile.close()
outfile.close()

