from pathlib import Path
import argparse
import re
import pandas as pd 
import os
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess

parser = argparse.ArgumentParser(description='Plot the number of reads for each trimming program for R1, R2 and the unpaired.')
parser.add_argument('-d', '--directory', help='the name or path of a directory', required=True)
parser.add_argument('-c', '--countout', help='a name for the output barplot for the read count', required=True)
parser.add_argument('-b', '--baseout', help='a name for the output barplot for the base count', required=True)
args = parser.parse_args()
#print(args.directory)

# proc = subprocess.Popen(["cat", "/etc/services"], stdout=subprocess.PIPE, shell=True)
# (out, err) = proc.communicate()
# print "program output:", out


data = []
pathlist = Path(args.directory).glob('**/*.fastq')
for path in pathlist:
  # because path is object not string
  path_in_str = str(path)
  filename=os.path.basename(path)
  name1,name2=os.path.splitext(filename)
  # print(name1,name2)
  proc = subprocess.Popen(["kseq_cnt "+path_in_str], stdout=subprocess.PIPE, shell=True)
  (out, err) = proc.communicate()
  #print("program output:", out.decode('ascii'))
  filepath, count, bases = out.strip().decode('ascii').split(",")
  if len(name1.split("_", 1))>1:
    name, sample, trimmer, paired = name1.split("_")
    # print(name,sample,trimmer,paired)
    if paired == "unpaired":
      data.append([name,trimmer,"singles",count,bases])
    elif paired == "R1":
      data.append([name,trimmer,"R1",count,bases])
#      elif paired == "R2":
#        data.append([name,trimmer,"R2",count])


df = pd.DataFrame(data, columns = ['Sample', 'Program', 'ReadType','Count','Bases']) 
df.Count=pd.to_numeric(df.Count)
df.Program=pd.Categorical(df.Program)
df = df.sort_values(by=['Program'])

#print(df)
g = sns.set(font_scale=0.7)
g = sns.FacetGrid(df, col="Sample",  row="ReadType", hue="Program")
g = (g.map(plt.bar, "Program", "Count").add_legend())
for ax in g.axes.flat:
    for label in ax.get_xticklabels():
        label.set_rotation(90)
g.savefig(args.countout)

df.Bases=pd.to_numeric(df.Bases)
g2 = sns.set(font_scale=0.7)
g2 = sns.FacetGrid(df, col="Sample",  row="ReadType", hue="Program")
g2 = (g2.map(plt.bar, "Program", "Bases").add_legend())
for ax in g2.axes.flat:
    for label in ax.get_xticklabels():
        label.set_rotation(90)
g2.savefig(args.baseout)


