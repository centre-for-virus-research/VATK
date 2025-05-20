from pathlib import Path
import argparse
import re
import pandas as pd 
import os
import matplotlib.pyplot as plt
import seaborn as sns

# These are the per-position base quality averages, calculated with different methods: 
# the _linear is just the average of the Phred read qualities; the _log first translate 
# the quality score to probability, then average of the probabilities, finally translate 
# back to Phred quality score. I think the _log should better reflect the overall base 
# quality than the linear, specially for data sets with large variations in qualities at 
# the same position.

# now using the log by default and linear with the option --linear

parser = argparse.ArgumentParser(description='Plot the quality as lines for each benchmarked trimming program in a folder.')
parser.add_argument('-d', '--directory', help='the name or path of a directory', required=True)
parser.add_argument('-o', '--output', help='a name for the output lineplot', required=True)
parser.add_argument('-l', '--linear', help='An argument that specifies if the linear quality or log quality should be used', required=False, action="store_true")
args = parser.parse_args()
#print(args.directory)

if args.linear:
  colR1=1
  colR2=3
  cols=1
else:
  colR1=2
  colR2=4
  cols=2


data = []
pathlist = Path(args.directory).glob('**/*paired.txt')
for path in pathlist:
  # because path is object not string
  path_in_str = str(path)
  filename=os.path.basename(path)
  name1,name2=os.path.splitext(filename)
  with open(path_in_str) as f:
    next(f)
    for line in f:
      vals = re.split(r'\t+', line)
      if len(name1.split("_", 1))>1:
        name, sample, trimmer, paired = name1.split("_")
        if paired == "unpaired":
          data.append([name,trimmer,"singles",vals[0],vals[cols]])
        else:
          data.append([name,trimmer,"R1",vals[0],vals[colR1]])
          data.append([name,trimmer,"R2",vals[0],vals[colR2]])


df = pd.DataFrame(data, columns = ['Sample', 'Program', 'ReadType','Position','Qual']) 
df.Position=pd.to_numeric(df.Position)
df.Qual=pd.to_numeric(df.Qual)
df.Program=pd.Categorical(df.Program)
df = df.sort_values(by=['Program','Position'])
#print(df)

g = sns.FacetGrid(df, col="Sample",  row="ReadType",hue="Program")
g = (g.map(plt.plot, "Position", "Qual").add_legend())
g.savefig(args.output)

