import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='sans-serif'
rcParams['font.sans-serif']=['Helvetica-normal']
import scipy.stats
from scipy.stats import gaussian_kde
import numpy as np

files = []
matrix={}
dir = "./data/"  # Data directory

for file in os.listdir(dir):  # loop through directory with files
   files.append(dir + file)

for i in range(len(files)):
   print('---- Adding file ' + files[i] + ' to output')
   matrix[files[i]] = []
   with open(files[i]) as tsvfile:
       for line in tsvfile:
        (chr, start, end, read_counts) = line.strip().split("\t")  # columns in files
        matrix[files[i]].append(read_counts)

df = pd.DataFrame.from_dict(matrix)
df = df.iloc[1:]  # if header is present
for column in df:
   df[column] = pd.to_numeric(df[column])

print(df)

# Median Absolute Deviation filtering
mad = df.mad(axis=0)
median= df.median(axis=0)
threshold=6  # set threshold
upperLimit= median + (threshold*mad)
df=df[df <= upperLimit]
df=df.dropna()
df=df[(df != 0).all(1)]  # remove values of 0

for column in df:
   max= df[column].max()
   df[column]=df[column]/max

# Choose columns to plot
x=df.ix[:,2]
y=df.ix[:,0]
xy= np.vstack([x,y])
z= gaussian_kde(xy)(xy)  # Gaussian kernel

# Statistics
spearmanCorr=scipy.stats.spearmanr(x,y)
pValue=scipy.stats.wilcoxon(x,y)

# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
# r2='$R^2$= ' + str(round(float(r_value**2),3))

# Plot
plt.style.use('ggplot')
plt.scatter(x=x,y=y,c=z,edgecolors='face')
#plt.text(x=0.8,y=0.95,s="rho= "+str(round(spearmanCorr[0],2)), fontsize=14)
plt.title("rho= "+str(round(spearmanCorr[0],2)), fontweight='bold')
plt.gca().set_xlim(left=0,right=1.1)
plt.gca().set_ylim(bottom=0, top=1.1)
plt.xlabel('G4_06', fontsize=16)
plt.ylabel('G4_07', fontsize=16)
plt.show()
