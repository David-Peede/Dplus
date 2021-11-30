#!/usr/bin/env python3
import numpy as np
from collections import defaultdict
import gzip
import sys
import statistics
import argparse
import re

#Functions
def get_population_dict(string,header,populations=["pop1","pop2","pop3","outgroup"]):
  dict={}
  for pop in populations:
    search_string=r'{}\[(.*?)\]'.format(pop)
    population_names=re.findall(search_string,string)[0].split(',')
    dict[pop]=[header.index(ind) for ind in population_names]
  return(dict)

def haplo(calls):
  """Convert the one character calls ('K','M','R','S','W','Y',or 'N') to nucleotides ('A','T','C','G') or 'N'. Input the list or array of
     calls for each diploid individual. Get list of nucleotides for each haploid in the input.
     Assume 'N' calls are filtered out.
  """
  output = []
  for call in calls:
    if call in "ACGT":
      output.append(call)
      output.append(call)
    elif call == "K":
      output.append("G")
      output.append("T")
    elif call == "M":
      output.append("A")
      output.append("C")
    elif call == "R":
      output.append("A")
      output.append("G")
    elif call == "S":
      output.append("C")
      output.append("G")
    elif call == "W":
      output.append("A")
      output.append("T")
    elif call == "Y":
      output.append("C")
      output.append("T")
    elif call=="N":
      output.append("N")
      output.append("N")
    else:
      print ("WARNING", call, "is not recognised as a valid base or ambiguous base")
      output.append("N")
      output.append("N")
  return output

def isBiallelic(alleles):
  """Is this allele biallelic? Input: a list of array  with all of the alleles for the four populations.
     Output a boolean True or False. True if allele biallelic. False if allele not biallelic.
  """
  return(np.unique(alleles).shape[0]==2)

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
  assert(alleles.shape[0]==2),"Alleles used are not biallelic!"
  #ancestral allele: fixed allele in outgroup or major allele in all four populations
  if (np.unique(dict["outgroup"]).shape[0]==1):
    ancestral_allele=np.unique(dict["outgroup"])[0]
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

def nucleotide_diversity_PerSite(population_alleles):
  """ Calculate the nucleotide diversity at a site for a given population.
      Input the list or array of alleles ["A","T","C","G"] for this biallelic site for the population.
      Output a float of the nucleotide diversity at the given site
  """
  #Get minor allele for ppulation
  alleles=np.unique(population_alleles)
  minor_allele=alleles[alleles!=get_major_allele(population_alleles)]
  #Count the number of minor alleles
  mA = population_alleles.count(minor_allele)
  #Count the total number of sequences
  seqTotal=float(np.shape(population_alleles)[0])
  #Calculate nucleotide diversity at this position
  nd = (2*mA*(seqTotal-mA))/(seqTotal*(seqTotal-1))
  return nd

def calculate_nd_stats(dict,nd_populations):
  results_dict=defaultdict(int)
  if "pop1" in nd_populations:
    results_dict["pop1_pi"]=nucleotide_diversity_PerSite(dict["pop1"])
  if "pop2" in nd_populations:
    results_dict["pop2_pi"]=nucleotide_diversity_PerSite(dict["pop2"])
  if "pop3" in nd_populations:
    results_dict["pop3_pi"]=nucleotide_diversity_PerSite(dict["pop3"])
  if "outgroup" in nd_populations:
    results_dict["outgroup_pi"]=nucleotide_diversity_PerSite(dict["outgroup"])
  if "allPops" in nd_populations:
    results_dict["allPops_pi"]=nucleotide_diversity_PerSite(sum(dict.values(),[]))
  return results_dict


##Parameters
parser = argparse.ArgumentParser()
#Specify file names
parser.add_argument("-o","--outfile",type=str,help="Path to outfile.")
parser.add_argument("-i","--infile",type=str,help="Path to data file.")

parser.add_argument("--d_statistic",action="store_true",help="To run the D statistic add option.")
parser.add_argument("--dplus_statistic",action="store_true",help="To run the D+ statistic add option.")
parser.add_argument("--fd_statistic",action="store_true",help="To run the fD statistic add option.")
parser.add_argument("--fdm_statistic",action="store_true",help="To run the fDM statistic add option.")
parser.add_argument("--df_statistic",action="store_true",help="To run the df statistic add option.")

parser.add_argument("--report",nargs='?',const=100,type=int,help="When to report progress on windows.")

parser.add_argument("-ws","--window_size",type=int,help="Size of the sliding window.")

parser.add_argument("-p","--populations",type=str,dest="population_string",help="Include the name for each population as following: 'pop1[name,name,name];pop2[name,name];pop3[name,name,name];outgroup[name,name]'.")

parser.add_argument("--verbose",action="store_true",help="Verbose mode.")

args=parser.parse_args()

#args.window_size=5000;args.report=100
#args.population_string="pop1[ag108,ag572,ag112,ag569];pop2[am216,am160,am48,am293];pop3[tiP86,tiP313,tiP84,tiP57];outgroup[hec273,eth67,ser202,par371]"
#args.test=True;args.d_statistic=True

statistics_list=[]
#Do D
if args.d_statistic:
  statistics_list.append("D")
#Do D+
if args.dplus_statistic:
  statistics_list.append("D+")
#Do fD
if args.fd_statistic:
  statistics_list.append("fD")
#Do fDM
if args.fdm_statistic:
  statistics_list.append("fDM")
#Do df
if args.df_statistic:
  statistics_list.append("df")

nd_populations=["pop2"]

nextReport=args.report
windows_tested=0

#Minimum number of sites per window
min_number_sites=3000

populations=["pop1","pop2","pop3","outgroup"]

#Open outfile
fout=open(args.outfile,"w")
outfile_header="{}\t".format("\t".join(["scaffold","start","stop","window_size","number_of_sites"]))
outfile_header+="{}\t".format("\t".join(statistics_list))
outfile_header+="\t".join(["{}_pi".format(pop) for pop in nd_populations])
fout.write("{}\n".format(outfile_header))

#Read in header first
header=True

with gzip.open(args.infile,"rt") as file:
  for line in file:
    if header:
      names=line.split()
      pop_index_dict=get_population_dict(args.population_string,names)
      header = False
      continue
    #Read in the file
    spline=np.asarray(line.split())
    #Start first window
    if not "scaffold" in locals():
      scaffold=spline[0]
      window_start=1
      stats_components=defaultdict(int)
      nd_stats=defaultdict(int)
      number_of_sites=0
    #Get information from site for current window
    number_of_sites+=1
    #Is a new window needed?
    if spline[0] != scaffold or int(spline[1]) > (window_start + args.window_size - 1):
      #Are there enough sites in the window
      if number_of_sites >= min_number_sites:
        #Calculate the stats for the window
        intro_stats=calculate_introgression_stats(stats_components,statistics_list)
        #Write new results to outfile
        outline="{}\t".format("\t".join([scaffold,str(window_start),str(window_start + args.window_size - 1),str(args.window_size),str(number_of_sites)]))
        outline+="{}\t".format("\t".join([str(intro_stats[key]) for key in statistics_list]))
        outline+="\t".join([str(nd_stats["{}_pi".format(key)]) for key in nd_populations])
        fout.write("{}\n".format(outline))
      #Restart the variables
      number_of_sites=0
      stats_components=defaultdict(int)
      nd_stats=defaultdict(int)
      #Same scaffold or start a new window?
      if spline[0] == scaffold:
        window_start += args.window_size
      else:
        scaffold=spline[0]
        window_start = 1
      #Window was tested
      windows_tested+=1
    #Is site biallelic
    if not isBiallelic(haplo(spline[sum(pop_index_dict.values(),[])][spline[sum(pop_index_dict.values(),[])]!="N"])):
      continue
    #Get the alleles for each population
    data_dict={}
    for pop in populations:
      calls_array=spline[pop_index_dict[pop]]
      data_dict[pop]=haplo(calls_array[calls_array!="N"])
    #Need alleles in all populations
    if any([np.shape(data_dict[pop])[0]<1 for pop in populations]):
      continue
    #Get the derived allele frequency of each population
    pop_derived_freq=get_derived_freq_per_population(data_dict)
    #Update the statistics numerators and denominators
    intro_dict=update_introgression_stats_sums(p_one=pop_derived_freq["pop1"],p_two=pop_derived_freq["pop2"]
                                           ,p_three=pop_derived_freq["pop3"],p_four=pop_derived_freq["outgroup"]
                                           ,list=statistics_list)
    stats_components={k: stats_components.get(k,0)+intro_dict.get(k,0) for k in intro_dict.keys()}
    #Update nucleotide diversity
    nd_site_dict=calculate_nd_stats(data_dict,nd_populations)
    nd_stats={k: nd_stats.get(k,0) + nd_site_dict.get(k,0) for k in nd_site_dict.keys()}
    #Report on number of windows tested
    if windows_tested == nextReport and args.verbose:
      print(windows_tested,"windows done.") #,time.asctime(time.localtime(time.time())))
      nextReport += args.report


#Close outfile

fout.close()
