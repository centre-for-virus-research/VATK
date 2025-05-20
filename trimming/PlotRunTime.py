from pathlib import Path
import argparse
import re
import pandas as pd 
import os
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot the runtime as a barplot for a set of benchmarks in a folder.')
parser.add_argument('-d', '--directory', help='the name or path of a directory', required=True)
parser.add_argument('-o', '--output', help='a name for the output barplot', required=True)
args = parser.parse_args()
#print(args.directory)

data = []
pathlist = Path(args.directory).glob('**/*.txt')
for path in pathlist:
  # because path is object not string
  path_in_str = str(path)
  filename=os.path.basename(path)
  name1,name2=os.path.splitext(filename)
  #print(path_in_str)
  f = open(path_in_str)
  lines = f.readlines()
  #print(lines[1])
  vals = re.split(r'\t+', lines[1])
  #print(vals[0])
  #benchmarks/Rana-07_S3.txt
  if len(name1.split(".", 1))>1:
    name, trimmer = name1.split(".", 1)
    data.append([name,trimmer,vals[0]])

#print(data)
 
#Create the pandas DataFrame 
df = pd.DataFrame(data, columns = ['Sample', 'Program', 'Time']) 
df.Program=pd.Categorical(df.Program)
df = df.sort_values(by=['Program'])
#print dataframe. 
df.Time=pd.to_numeric(df.Time)
sns.set(font_scale=0.7) 
sns_plot = sns.catplot(x="Program", y="Time", hue="Sample", kind="bar", data=df)
for ax in sns_plot.axes.flat:
    for label in ax.get_xticklabels():
        label.set_rotation(90)

sns_plot.savefig(args.output)

