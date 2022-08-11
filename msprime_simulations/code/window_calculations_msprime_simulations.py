#!/usr/bin/env python3
import numpy as np
import tskit
from collections import defaultdict
import gzip
import sys
import os
import statistics
import argparse
import re

#Functions
def get_population_dict(string,populations=["p1","p2","p3","p4"]):
  dict={}
  for pop in populations:
    search_string=r'{}_\d+'.format(pop)
    population_names=re.findall(search_string,string)
    dict[pop]=[string.split().index(ind) for ind in population_names]
  return(dict)


def isBiallelic(reference,alternative):
  nucleotides=["A","T","C","G"]
  return((nucleotides.count(reference)==1)&(nucleotides.count(alternative)==1))

def get_major_allele(alleles):
  """Get the major allele of a population. Input a list or array with alleles ["A","T","C","G"].
     Outputs a string with a nucleotide letter "A","T","C" or"G".
  """
  try:
    major_allele=statistics.mode(alleles)
  except statistics.StatisticsError:
    major_allele=np.unique(alleles)[0]
  return(major_allele)

def get_derived_allele(dict):
  """Get derived allele of the populations. Input a dictionary with the alleles ["A","T","C","G"] for the
     biallelic site in each of the four populations. Return a string with a nucleotide letter "A","T","C" or"G".
  """
  #alleles for all populations
  population_alleles=np.asarray(sum(dict.values(),[]))
  alleles=np.unique(population_alleles)
  #Only biallelic alleles
  assert(alleles.shape[0]<=2),"Alleles used are not biallelic!"
  #ancestral allele: fixed allele in outgroup or major allele in all four populations
  if (np.unique(dict["p4"]).shape[0]==1):
    ancestral_allele=np.unique(dict["p4"])[0]
  else:
    ancestral_allele=get_major_allele(population_alleles)
  #derived allele
  derived_allele = alleles[alleles!=ancestral_allele][0]
  return derived_allele

def get_derived_freq_per_population(dict):
  """Get derived allele frequency for each population. Input a dictionary with the alleles ["A","T","C","G"] for the biallelic
    site in each of the four populations. Output a dictionary with the float derived allele frequency for each population.
  """
  derived_allele=get_derived_allele(dict)
  freq_dict={}
  for pop in dict.keys():
    freq_dict[pop]=float(dict[pop].count(derived_allele))/np.shape(dict[pop])[0]
  return freq_dict

def update_introgression_stats_sums(p_one,p_two,p_three,p_four,list):
  dict=defaultdict(int)
  if "D" in list:
    #D statistic
    dict["d_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    dict["d_denominator"]=((1-p_one)*p_two*p_three*(1-p_four))+(p_one*(1-p_two)*p_three*(1-p_four))
  if "D+" in list:
    #D+ statistic
    dict["dplus_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    dict["dplus_numerator"]+=(p_one*(1-p_two)*(1-p_three)*(1-p_four))-((1-p_one)*p_two*(1-p_three)*(1-p_four))
    dict["dplus_denominator"]=((1-p_one)*p_two*p_three*(1-p_four))+(p_one*(1-p_two)*p_three*(1-p_four))
    dict["dplus_denominator"]+=(p_one*(1-p_two)*(1-p_three)*(1-p_four))+((1-p_one)*p_two*(1-p_three)*(1-p_four))
  if "Dancestral" in list:
    #D+ statistic
    dict["dancestral_numerator"]=(p_one*(1-p_two)*(1-p_three)*(1-p_four))-((1-p_one)*p_two*(1-p_three)*(1-p_four))
    dict["dancestral_denominator"]=(p_one*(1-p_two)*(1-p_three)*(1-p_four))+((1-p_one)*p_two*(1-p_three)*(1-p_four))
  if "fD" in list:
    #fD statistic
    dict["fd_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    #fD uses the max(p2,p3) as donor population
    pd=max(p_two,p_three)
    dict["fd_denominator"]=((1-p_one)*pd*pd*(1-p_four))-(p_one*(1-pd)*pd*(1-p_four))
  if "fDM" in list:
    #fDM statistic
    dict["fdm_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    #fDM has the same denominator as fD when p_two is greater than or equal to p3
    if p_one>=p_two:
      pd_m=max(p_two,p_three)
      dict["fdm_denominator"]=((1-p_one)*pd_m*pd_m*(1-p_four))-(p_one*(1-pd_m)*pd_m*(1-p_four))
    elif p_one<p_two:
      pd_m=max(p_one,p_three)
      dict["fdm_denominator"]=(pd_m*(1-p_two)*pd_m*(1-p_four))-((1-pd_m)*p_two*pd_m*(1-p_four))
  if "df" in list:
    #Df statsitic
    dict["df_numerator"]=((1-p_one)*p_two*p_three)-(p_one*(1-p_two)*p_three)
    dict["df_denominator"]=(2*p_one*p_two*p_three)+((1-p_one)*p_two*p_three)+(p_one*(1-p_two)*p_three)
  return(dict)

def allele_site_pattern(pop1_allele,pop2_allele,pop3_allele,pop4_allele):
  #BAAA or ABAA
  if(pop4_allele==pop3_allele):
    #BAAA
    if((pop3_allele==pop2_allele) & (pop1_allele!=pop2_allele)):
      return "BAAA"
    #ABAA
    if((pop3_allele==pop1_allele) & (pop1_allele!=pop2_allele)):
      return "ABAA"
  #ABBA or BABA
  if(pop4_allele!=pop3_allele):
    #ABBA
    if((pop3_allele==pop2_allele) & (pop1_allele!=pop2_allele)):
      return "ABBA"
    #ABAA
    if((pop3_allele==pop1_allele) & (pop1_allele!=pop2_allele)):
      return "BABA"
  return None

def calculate_introgression_stats(dict,statistics):
  results=defaultdict(int)
  if "D" in statistics:
    #D statistic
    try:
      results["D"]=dict["d_numerator"]/dict["d_denominator"]
    except ZeroDivisionError:
      results["D"]=np.nan
  if "D+" in statistics:
    #D+ statistic
    try:
      results["D+"]=dict["dplus_numerator"]/dict["dplus_denominator"]
    except ZeroDivisionError:
      results["D+"]=np.nan
  if "Dancestral" in statistics:
    #Dancestral statistic
    try:
      results["Dancestral"]=dict["dancestral_numerator"]/dict["dancestral_denominator"]
    except ZeroDivisionError:
      results["Dancestral"]=np.nan
  if "fD" in statistics:
    #fD statistic
    try:
      results["fD"]=dict["fd_numerator"]/dict["fd_denominator"]
    except ZeroDivisionError:
      results["fD"]=np.nan
  if "fDM" in statistics:
    #fDM statistic
    try:
      results["fDM"]=dict["fdm_numerator"]/dict["fdm_denominator"]
    except ZeroDivisionError:
      results["fDM"]=np.nan
  if "df" in statistics:
    #Df statsitic
    try:
      results["df"]=dict["df_numerator"]/dict["df_denominator"]
    except ZeroDivisionError:
      results["df"]=np.nan
  return results

def number_intro_bases(window_start,window_stop,intro_start,intro_stop):
  bases=0
  #Introgressed region 100% overlap
  if(window_start >= intro_start and window_stop <= intro_stop):
    bases=window_stop - window_start
  #Introgressed region overlap at beginning
  if(window_start >= intro_start and window_stop >= intro_stop and intro_stop > window_start):
    bases = intro_stop-window_start + 1
  #Introgressed region overlap at end
  if(window_start <= intro_start and window_stop <= intro_stop and intro_start < window_stop):
    bases = window_stop - intro_start + 1
  #Introgressed region within window
  if(window_start <= intro_start and window_stop >= intro_stop):
    bases = intro_stop - intro_start + 1
  return(bases)

##Parameters
default_sequence_length=2e7
parser = argparse.ArgumentParser()

parser.add_argument("--tree_sequence_file",type=str,help="Path to tree sequence file to be output.")
parser.add_argument("--vcf_file",type=str,help="Path to vcf file to be output.")
parser.add_argument("-o","--outfile",type=str,help="Path to results file to store introgression statistics of vcf file.")

parser.add_argument("-f",type=int,help="Admixture proportion f. Write as percentage. For example f = 0.03 is -f 3")
parser.add_argument("-s","--seed",type=int,help="Seed number of simulation.")

parser.add_argument("--report",nargs='?',const=100,type=int,help="When to report progress on windows.")

parser.add_argument("-ws","--window_size",type=int,help="Size of sliding windows.")
parser.add_argument("--sequence_length",type=int,nargs='?',default=default_sequence_length,help="The size of the sequence for each population. Default is 2 MB")

parser.add_argument("--introgressed_bases",dest="count_introgressed",action="store_true",help="Turn on count of introgressed bases.")
parser.add_argument("--haplotype",dest="use_haplotype",action="store_true",help="Calculate statistics on a haplotype per population.")

parser.add_argument("--d_statistic",action="store_true",help="To run the D statistic add option.")
parser.add_argument("--dplus_statistic",action="store_true",help="To run the D+ statistic add option.")
parser.add_argument("--fd_statistic",action="store_true",help="To run the fD statistic add option.")
parser.add_argument("--fdm_statistic",action="store_true",help="To run the fDM statistic add option.")
parser.add_argument("--df_statistic",action="store_true",help="To run the df statistic add option.")

args=parser.parse_args()

statistics_list=[]
sites=[]
#Do D
if args.d_statistic:
  statistics_list.append("D")
  sites.append("ABBA")
  sites.append("BABA")

#Do D+
if args.dplus_statistic:
  statistics_list.append("D+")
  sites.append("BAAA")
  sites.append("ABAA")

#Do fD
if args.fd_statistic:
  statistics_list.append("fD")

#Do fDM
if args.fdm_statistic:
  statistics_list.append("fDM")

#Do df
if args.df_statistic:
  statistics_list.append("df")

#Get introgressed bases
if args.count_introgressed:
  ts=tskit.load(args.tree_sequence_file)
  pop_ids=defaultdict()
  for pop in ts.populations():
    pop_ids[pop.metadata['name']]=pop.id
  migrations=[]
  for migration in ts.migrations():
    if migration.dest==pop_ids["p_3"]:
      migrations.append((migration.left,migration.right))


nextReport=args.report
windows_tested=0

populations=["p1","p2","p3","p4"]
deliminator="|"
start_calculations=True

#Open outfile
outfile_delim="\t"
header = outfile_delim.join(["f","seed","start","stop","window_size","number_of_sites"])
if args.use_haplotype:
  header += "{}{}".format(outfile_delim,outfile_delim.join(sites))
header += "{}{}".format(outfile_delim,outfile_delim.join(statistics_list))
if args.count_introgressed:
  header+="{}{}".format(outfile_delim,outfile_delim.join(["introgressed_bases","introgressed_percentage"]))
if os.path.exists(args.outfile):
  fout=open(args.outfile,"a+")
else:
  fout=open(args.outfile,"w")
  fout.write("{}\n".format(header))


with gzip.open(args.vcf_file,"rt") as file:
  for line in file:
    if (line.startswith("##")):
      continue
    if (line.startswith("#")):
      pop_index_dict=get_population_dict(line)
      continue
    #Read in the file
    spline=np.asarray(line.split())
    pos=int(spline[1])
    #Start first window
    #if not "window_start" in locals():
    if start_calculations:
      start_calculations=False
      window_start=1
      if args.use_haplotype:
        site_dict=defaultdict(int)
      else:
        stats_components=defaultdict(int)
      introgressed_bases=0
      number_of_sites=0
    #Get information from site for current window
    number_of_sites+=1
    #Is a new window needed?
    if (pos > (window_start + args.window_size - 1)):
      #Write new results to outfile
      outline=outfile_delim.join([str(args.f),str(args.seed),str(window_start),str(window_start + args.window_size - 1),str(args.window_size),str(number_of_sites)])
      if args.use_haplotype:
        stats_dict=calculate_introgression_stats(site_dict,statistics_list)
        outline+="{}{}".format(outfile_delim,outfile_delim.join([str(site_dict[key]) for key in sites]))
      else:
        #Calculate the stats for the window
        intro_stats=calculate_introgression_stats(stats_components,statistics_list)
        outline+="{}{}".format(outfile_delim,outfile_delim.join([str(intro_stats[key]) for key in statistics_list]))
      #Count the number of introgressed bases if applicable
      if args.count_introgressed:
        introgressed_bases=np.sum([number_intro_bases(window_start,(window_start+args.window_size-1),migration[0],migration[1]) for migration in migrations])
        outline+="{}{}".format(outfile_delim,outfile_delim.join([str(introgressed_bases),str(introgressed_bases/float(args.window_size))]))
      fout.write("{}\n".format(outline))
      #Restart the variables
      number_of_sites=0
      if args.use_haplotype:
        site_dict=defaultdict(int)
      else:
        stats_components=defaultdict(int)
      window_start += args.window_size
      #Window was tested
      windows_tested+=1
    #Is site biallelic
    if not isBiallelic(spline[3],spline[4]):
      continue
    #Get the alleles for each population
    data_dict={}
    for pop in populations:
      data_dict[pop]=[int(allele) for genotype in spline[pop_index_dict[pop]] for allele in genotype.split(deliminator)]
    #Need alleles in all populations
    #if any([np.shape(data_dict[pop])[0]<1 for pop in populations]):
    #Need variants
    if not (np.unique(sum(data_dict.values(),[])).shape[0]==2):
      continue
    #If using haplotype:
    if args.use_haplotype:
      #Get allele site pattern
      site_key=allele_site_pattern(pop1_allele=spline[pop_index_dict["p1"]],pop2_allele=spline[pop_index_dict["p2"]]
                             ,pop3_allele=spline[pop_index_dict["p3"]],pop4_allele=spline[pop_index_dict["p4"]])
      if not site_key:
        continue
      site_dict[site_key]+=1
    else:
      #Get the derived allele frequency of each population
      pop_derived_freq=get_derived_freq_per_population(data_dict)
      #Update the statistics numerators and denominators
      intro_dict=update_introgression_stats_sums(p_one=pop_derived_freq["p1"],p_two=pop_derived_freq["p2"]
                                           ,p_three=pop_derived_freq["p3"],p_four=pop_derived_freq["p4"]
                                           ,list=statistics_list)
      stats_components={k: stats_components.get(k,0)+intro_dict.get(k,0) for k in intro_dict.keys()}
    #Report on number of windows tested
    if windows_tested == nextReport and args.verbose:
      print(windows_tested,"windows done.") #,time.asctime(time.localtime(time.time())))
      nextReport += args.report

#Calculate last window & write results
#Write new results to outfile
outline="{}".format(outfile_delim.join([str(args.f),str(args.seed),str(window_start),str(window_start + args.window_size - 1),str(args.window_size),str(number_of_sites)]))
#Calculate last window
if args.use_haplotype:
  stats_dict=calculate_introgression_stats(site_dict,statistics_list)
  outline+="{}{}".format(outfile_delim,outfile_delim.join([str(site_dict[key]) for key in sites]))
else:
  intro_stats=calculate_introgression_stats(stats_components,statistics_list)
  outline+="{}{}".format(outfile_delim,outfile_delim.join([str(intro_stats[key]) for key in statistics_list]))
#Count the number of introgressed bases if applicable
if args.count_introgressed:
  introgressed_bases=np.sum([number_intro_bases(window_start,(window_start+args.window_size-1),migration[0],migration[1]) for migration in migrations])
  outline+="{}{}".format(outfile_delim,outfile_delim.join([str(introgressed_bases),str(introgressed_bases/float(args.window_size))]))
fout.write("{}\n".format(outline))

#Close outfile
fout.close()

