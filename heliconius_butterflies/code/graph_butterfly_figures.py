import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

##Parameters
parser = argparse.ArgumentParser()
#Specify file names
parser.add_argument("-o","--outfile",type=str,help="Path to outfile directory.")
parser.add_argument("-i","--infile",type=str,help="Path to results file.")

args=parser.parse_args()

ds=pd.read_csv(args.infile,sep="\t")
min_number_of_sites=3000
df=ds.loc[ds.number_of_sites>=min_number_of_sites]

hmb_scaffold="HE670865";hmb_start=300000;hmb_stop=450000
hmyb_scaffold="HE667780";hmyb_start=650000;hmyb_stop=900000

hmb = df.loc[(df.scaffold==hmb_scaffold)&((df.start >= hmb_start)&(df.stop <= hmb_stop))]
hmyb = df.loc[(df.scaffold==hmyb_scaffold)&((df.start >= hmyb_start)&(df.stop <= hmyb_stop))]
df = df[~((df.scaffold==hmb_scaffold)&((df.start >= hmb_start)&(df.stop <= hmb_stop)))
        & ~((df.scaffold==hmyb_scaffold)&((df.start >= hmyb_start)&(df.stop <= hmyb_stop)))]

#Axes specifications
font = {'family': 'Arial', 'color':  'black','weight': 'normal','size': 14}
x_title='Nucleotide Diversity {}'.format(r"$(\pi)$")

#D figure
stat="D"
d_outfile="{}/D_vs_nucleotide_diversity_butterflies".format(args.outfile)
fig, ax = plt.subplots(figsize = ( 10 , 5 ))

sns.scatterplot(x="pop2_pi",y=stat,data=df,alpha=0.05,color="black",s=50)
sns.scatterplot(x="pop2_pi",y=stat,data=hmb,alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60)
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb,alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60)
ax.set_xlabel(x_title,fontdict=font)
ax.set_ylabel(r"$D$",fontdict=font)
plt.tight_layout()
plt.savefig("{}.png".format(d_outfile),dpi=400)
ax.clear()

#D+ figure
stat="D+"
dplus_outfile="{}/Dplus_vs_nucleotide_diversity_butterflies".format(args.outfile)
fig, ax = plt.subplots(figsize = ( 10 , 5 ))

sns.scatterplot(x="pop2_pi",y=stat,data=df,alpha=0.05,color="black",s=50)
sns.scatterplot(x="pop2_pi",y=stat,data=hmb,alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60)
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb,alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60)
ax.set_xlabel(x_title,fontdict=font)
ax.set_ylabel(r"$D^+$",fontdict=font)
plt.tight_layout()
plt.savefig("{}.png".format(dplus_outfile),dpi=400)
ax.clear()

d_dplus_outfile="{}/D-Dplus_vs_nucleotide_diversity_butterflies".format(args.outfile)
fig, ax = plt.subplots(ncols=2,figsize = ( 15 , 5 ),sharey=False)
stat="D"
sns.scatterplot(x="pop2_pi",y=stat,data=df,alpha=0.05,color="black",s=50,ax=ax[0])
sns.scatterplot(x="pop2_pi",y=stat,data=hmb,alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[0])
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb,alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[0])
ax[0].set_xlabel(x_title,fontdict=font)
ax[0].set_ylabel(r"$D$",fontdict=font)

stat="D+"
sns.scatterplot(x="pop2_pi",y=stat,data=df,alpha=0.05,color="black",s=50,ax=ax[1])
sns.scatterplot(x="pop2_pi",y=stat,data=hmb,alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[1])
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb,alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[1])
ax[1].set_xlabel(x_title,fontdict=font)
ax[1].set_ylabel(r"$D^+$",fontdict=font)

plt.tight_layout()
plt.savefig("{}.png".format(d_dplus_outfile),dpi=400)
fig.clear()


