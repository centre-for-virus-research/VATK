import os
import sys

os.chdir("referencesFolder")
os.system("ls *_virhosfilt >viruslistfile")
infile = open("viruslistfile")
outfile = open("viruslistfile.txt","w")
while True:
    line = ((infile.readline().rstrip()).split("_virhosfilt"))[0]
    if not line:
        break
    outfile.write(line+"\n")

infile.close()
outfile.close()
os.system("mv viruslistfile.txt ../")
os.system("rm viruslistfile")
os.chdir("..")