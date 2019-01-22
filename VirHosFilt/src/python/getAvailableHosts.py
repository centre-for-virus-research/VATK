import os
import sys

hostFolder = sys.argv[1]



os.system("ls "+hostFolder+"/hostsFolder/*_virhosfilt >hostlistfile")
infile = open("hostlistfile")
outfile = open("hostlistfile.txt","w")
while True:
    line = ((infile.readline().rstrip()).split("_virhosfilt"))[0]
    
    if not line:
        break
    line2 = (line.split("/"))[-1]
    outfile.write(line2+"\n")

infile.close()
outfile.close()
#os.system("rm hostlistfile.txt")
#os.system("rm hostlistfile")

    