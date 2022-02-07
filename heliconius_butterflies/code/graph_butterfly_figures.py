import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import logging

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
font = {'family': 'Arial', 'color':  'black','weight': 'normal','size': 12}
dpi_size=550
x_title='Nucleotide Diversity {}'.format(r"$(\pi)$")

logging.getLogger('matplotlib.font_manager').disabled = True

###Figure 7
figure_seven_outfile="{}fig7".format(args.outfile)
fig, ax = plt.subplots(nrows=2,ncols=2,figsize = ( 15,12),sharey=False)

#D
stat="D";axes_index_first=0;axes_index_second=0
sns.scatterplot(x="pop2_pi",y=stat,data=df,alpha=0.05,color="black",s=50,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmb,alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb,alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
ax[axes_index_first][axes_index_second].set_xlabel(x_title,fontdict=font)
ax[axes_index_first][axes_index_second].set_ylabel(r"$D$",fontdict=font)


#fD
stat="fD";axes_index_first=0;axes_index_second=1
sns.scatterplot(x="pop2_pi",y=stat,data=df.loc[df.D>0],alpha=0.05,color="black",s=50,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmb.loc[hmb.D>0],alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb.loc[hmyb.D>0],alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
ax[axes_index_first][axes_index_second].set_xlabel(x_title,fontdict=font)
ax[axes_index_first][axes_index_second].set_ylabel(r"$f_D$",fontdict=font)

#Df
stat="df";axes_index_first=1;axes_index_second=0
sns.scatterplot(x="pop2_pi",y=stat,data=df,alpha=0.05,color="black",s=50,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmb,alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb,alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
ax[axes_index_first][axes_index_second].set_xlabel(x_title,fontdict=font)
ax[axes_index_first][axes_index_second].set_ylabel(r"$d_{f}$",fontdict=font)

#D+
stat="D+";axes_index_first=1;axes_index_second=1
sns.scatterplot(x="pop2_pi",y=stat,data=df,alpha=0.05,color="black",s=50,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmb,alpha=1,markers='o',edgecolor='red',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
sns.scatterplot(x="pop2_pi",y=stat,data=hmyb,alpha=1,markers='o',edgecolor='yellow',marker="$\circ$"
               ,facecolors='none',s=60,ax=ax[axes_index_first][axes_index_second])
ax[axes_index_first][axes_index_second].set_xlabel(x_title,fontdict=font)
ax[axes_index_first][axes_index_second].set_ylabel(r"$D^+$",fontdict=font)

plt.tight_layout()
plt.savefig("{}.svg".format(figure_seven_outfile),dpi=dpi_size)
plt.savefig("{}.png".format(figure_seven_outfile),dpi=dpi_size)
fig.clear()
