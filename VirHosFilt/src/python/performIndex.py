import os
import sys

fastaFile = sys.argv[1]
indexSuffix = sys.argv[2]
os = sys.argv[3]
folder = sys.argv[4]



os.system("./binaries/"+os+"/bowtie2-build "+fastaFile+" "+indexSuffix)
os.system("mv "+indexSuffix+"* ./"+folder+"Folder/")




print "\n\nIndex created. Press any key to close this window...."
sys.stdin.read(1)