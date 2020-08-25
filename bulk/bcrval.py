#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Compare TRUST4 and MiXCR on the RNA-seq with BCR-seq validations from iRep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import numpy as np
import scipy as sp
import seaborn as sns
import pandas as pd
from itertools import combinations 
from scipy.cluster import hierarchy
#import matplotlib.style
import string


# In[3]:


irepPath = "/Users/lsong/Song/MyProjects/bcr/bcrval/irep"
trust4Path = "/Users/lsong/Song/MyProjects/bcr/bcrval/TRUST4"
trust4FqPath = "/Users/lsong/Song/MyProjects/bcr/bcrval/TRUST4_FQ"
mixcrPath = "/Users/lsong/Song/MyProjects/bcr/bcrval/mixcr"
mixcrNoFilterPath = "/Users/lsong/Song/MyProjects/bcr/bcrval/mixcr_NoFilter"


# In[4]:


def GetChainType( v, j, c ):
    s = ""
    if (v == '-' or v == ""):
        v = "*"
    if (j == '-' or j == ""):
        j = "*"
    if (c == '-' or c == ""):
        c = "*"
        
    if ( v != "*" ):
        s = v 
    elif ( j != "*"):
        s = j
    elif ( c != "*" ):
        s = c
    else:
        return -1

    if ( s[0:3] == "IGH" ):
        return 0
    elif ( s[0:3] == "IGK" ):
        return 1
    elif ( s[0:3] == "IGL" ):
        return 2 
    elif ( s[0:3] == "TRA" ):
        return 3
    elif ( s[0:3] == "TRB" ):
        return 4
    elif ( s[0:3] == "TRG" ):
        return 5
    elif ( s[0:3] == "TRD" ):
        return 6
    
    return -1 

chainTypeToName = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD" )

def GetMainGeneName( g ): # remove the allele type from the gene name
    fields = g.split( "*", 1 )
    return fields[0]


# In[5]:


def ParseLine(line, tool):
    ret = ["-", "-", "-", "-", "", 0.0] # V, D, J, C, cdr3nt, count
    if (tool == "trust4"): # parse TRUST4's report
        if (line[0] == "#"):
            return ret
        cols = line.split("\t") ;
        if (cols[3] == "partial"):
            return ret
        if (GetChainType(cols[4], cols[6], cols[7]) != 0):
            return ret
        for g in [4, 5, 6, 7]:
            if (cols[g] != "*"):
                ret[g - 4] = GetMainGeneName(cols[g])
        ret[4] = cols[2][3:-3]
        ret[5] = float(cols[0])
        
    elif (tool == "mixcr"): # parse mixcr's result
        if (line[0] == "#"):
            return ret
        cols = line.split("\t") ;
        if (GetChainType(cols[5], cols[7], cols[8]) != 0):
            return ret
        for g in [5, 6, 7, 8]:
            if (cols[g] != ""):
                ret[g - 5] = GetMainGeneName(cols[g])
        ret[4] = cols[23][3:-3]
        ret[5] = float(cols[1])
        
    elif (tool == "irep"): # Parse irep's result
        cols = line.split(",")
        for g in [1, 2, 3, 4]:
            if (cols[g] != "-"):
                ret[g - 1] = GetMainGeneName(cols[g][1:]) 
        ret[4] = cols[5].upper()
        ret[5] = float(cols[6])
    return ret


# In[6]:


# Return the list of qualified lines from each tool
def ReadFile(file, tool, onlyCDR3):
    retList = []
    fp = open(file, "r")
    line = fp.readline() # all the tools have headers for now.
    for line in fp:
        line = line.rstrip()
        l = ParseLine(line, tool)
        if (l[3] == ""): # empty CDR3
            continue
        l = [l[0], l[2], l[3], l[4], l[5]] # remove D gene
        if (onlyCDR3 == False and tool == "irep"):
            for g in [0, 1, 2]: # Only use irep report with all the component 
                if (l[g] == "-"):
                    continue
        # For isotype ignore the numerical numbers, since irep put all IGHA in IGHA2
        if (l[2][-1].isnumeric()):
            l[2] = l[2][:-1]
        
        if (onlyCDR3):
            l[0] = l[1] = l[2] = "-"
        retList.append(l) ;
    
    # Coalesce same entries
    tmpList = sorted(retList)
    retList = []
    i = 0 
    while ( i < len(tmpList) ):
        j = i + 1
        while (j < len(tmpList)):
            if (tmpList[j][:4] != tmpList[i][:4]):
                break
            j += 1
        counts = [] 
        for k in range(i, j):
            counts.append(tmpList[k][4])
        l = tmpList[i][:]
        l[4] = max(counts) # Pick the largest one as representative.
        retList.append(l)
        i = j
    
    # Compute the fraction
    readSum = 0
    for b in retList:
        readSum += b[4]
    for i in range(len(retList)):
        retList[i].append( retList[i][4] / readSum ) # index 5 is for fraction
    
    return retList


# In[7]:


def GetSensitivity(truthList, positiveList):
    trueSet = set({})
    positiveSet = set({})
    for b in truthList:
        trueSet.add((b[0], b[1], b[2], b[3])) 

    if (len(trueSet) == 0):
        return 0
    for b in positiveList:
        #if (b[5] < minFraction):
        #    continue
        positiveSet.add((b[0], b[1], b[2], b[3])) 
    return len(trueSet.intersection(positiveSet)) / len(trueSet) * 100.0, len(trueSet.intersection(positiveSet)),            len(positiveSet)

def GetPrecision(truthList, positiveList, onlyCDR3):
    trueSet = set({})
    positiveSet = set({})
    for b in truthList:
        trueSet.add((b[0], b[1], b[2], b[3])) 

    if (len(trueSet) == 0):
        return 0
    for b in positiveList:
        #if (b[5] <= minFraction):
        #    continue
        if (onlyCDR3 == False and (b[0] == "-" or b[1] == "-" or b[2] == "-")):
            continue
        positiveSet.add((b[0], b[1], b[2], b[3])) 
    #print(len(trueSet.intersection(positiveSet)), len(positiveSet))
    return len(trueSet.intersection(positiveSet)) / len(positiveSet) * 100.0, len(trueSet.intersection(positiveSet))

def GetCorrelation(truthList, positiveList, onlyCDR3):
    trueSet = set({})
    positiveSet = set({})
    for b in truthList:
        trueSet.add((b[0], b[1], b[2], b[3])) 

    if (len(trueSet) == 0):
        return 0
    for b in positiveList:
        if (onlyCDR3 == False and (b[0] == "-" or b[1] == "-" or b[2] == "-")):
            continue
        positiveSet.add((b[0], b[1], b[2], b[3])) 
    
    intersectionSet = trueSet.intersection(positiveSet)
    trueDict = {}
    positiveDict = {}
    for b in truthList:
        if ((b[0], b[1], b[2], b[3]) in intersectionSet):
            trueDict[(b[0], b[1], b[2], b[3])] = b[4]
    for b in positiveList:
        if ((b[0], b[1], b[2], b[3]) in intersectionSet):
            positiveDict[(b[0], b[1], b[2], b[3])] = b[4]
    t = []
    p = []
    for key in trueDict.keys():
        t.append(trueDict[key])
        p.append(positiveDict[key])
    return sp.stats.pearsonr(t, p)[0]


# In[13]:


irepListAll = {}
trust4ListAll = {}
mixcrListAll = {}
mixcrNoFilterListAll = {}
trust4FqListAll = {}

for rdlen in ["150"]:
    irepListAll[rdlen] = {}
    trust4ListAll[rdlen] = {}
    mixcrListAll[rdlen] = {}
    for sample in ["20", "83", "94", "97", "116", "122"]:
        irepListAll[rdlen][sample] = ReadFile(irepPath+"/FZ-"+sample+".csv", "irep", False)
        trust4ListAll[rdlen][sample] = ReadFile(trust4Path+"/FZ-"+sample+"_"+rdlen+"bp_report.tsv", "trust4", False)
        mixcrListAll[rdlen][sample] = ReadFile(mixcrPath+"/FZ-"+sample+"_"+rdlen+"bp_full_clones.txt", "mixcr", False)

trust4FqListAll["150"] = {}
mixcrNoFilterListAll["150"] = {}
for sample in ["20", "83", "94", "97", "116", "122"]:
    trust4FqListAll["150"][sample] = ReadFile(trust4FqPath+"/FZ-"+sample+"_"+rdlen+"bp_report.tsv", "trust4", False)
    mixcrNoFilterListAll["150"][sample] = ReadFile(mixcrNoFilterPath+"/FZ-"+sample+"_full_clones.txt", "mixcr", False)





# evaluation on the all-component (V, J, C, CDR3nt) performance
sns.set( font_scale=2, style="white")
for rdlen in ["150"]:
    fig, axes = plt.subplots( 1, 6, figsize=(35, 10))
                             #gridspec_kw={'width_ratios': widthRatio} )
    #fig.tight_layout( w_pad=2 )
    pr = pd.DataFrame(columns=["Sensitivity", "recalls", "Precision", "topN", "Method", "Sample", "Correlation"]) # precision-recall 
    axIdx = 0 
    for sample in ["20", "83", "94", "97", "116", "122"]:
    #for sample in ["116"]:
        #irepList = ReadFile(irepPath+"/FZ-"+sample+".csv", "irep", False)
        #trust4List = ReadFile(trust4Path+"/FZ-"+sample+"_"+rdlen+"bp_report.tsv", "trust4", False)
        #mixcrList = ReadFile(mixcrPath+"/FZ-"+sample+"_"+rdlen+"bp_full_clones.txt", "mixcr", False)
        irepList = irepListAll[rdlen][sample]
        trust4List = trust4ListAll[rdlen][sample]
        trust4FqList = trust4ListAll[rdlen][sample]
        mixcrList = mixcrListAll[rdlen][sample]
        mixcrNoFilterList = mixcrNoFilterListAll[rdlen][sample]
        
        irepList.sort(key = lambda b:b[4], reverse=True) # sort by abundance
        trust4List.sort(key = lambda b:b[4], reverse=True)
        trust4FqList.sort(key = lambda b:b[4], reverse=True)
        mixcrList.sort(key = lambda b:b[4], reverse=True)
        mixcrNoFilterList.sort(key = lambda b:b[4], reverse=True)
        
        samplePr = pd.DataFrame()
        topN = [100]
        maxLen = max([len(trust4List), len(trust4FqList), len(mixcrList), len(mixcrNoFilterList)])
        for i in range(500, maxLen, 500):
            topN.append(i)
        topN.append(-1)
        trust4Recall = [GetSensitivity(irepList, trust4List[:n]) for n in topN]
        trust4FqRecall = [GetSensitivity(irepList, trust4FqList[:n]) for n in topN]
        mixcrRecall = [GetSensitivity(irepList, mixcrList[:n]) for n in topN]  
        mixcrNoFilterRecall = [GetSensitivity(irepList, mixcrNoFilterList[:n]) for n in topN]  
        trust4Precision = [GetPrecision(irepList, trust4List[:n], False)[0] for n in topN] 
        trust4FqPrecision = [GetPrecision(irepList, trust4FqList[:n], False)[0] for n in topN] 
        mixcrPrecision = [GetPrecision(irepList, mixcrList[:n], False)[0] for n in topN] 
        mixcrNoFilterPrecision = [GetPrecision(irepList, mixcrNoFilterList[:n], False)[0] for n in topN] 
        
        trust4Correlation = GetCorrelation(irepList, trust4List, False) 
        trust4FqCorrelation = GetCorrelation(irepList, trust4FqList, False) 
        mixcrCorrelation = GetCorrelation(irepList, mixcrList, False)
        mixcrNoFilterCorrelation = GetCorrelation(irepList, mixcrNoFilterList, False)
        
        samplePr["topN"] = topN + topN + topN + topN
        samplePr["Sensitivity"] = [x[0] for x in trust4Recall] + [x[0] for x in trust4FqRecall] + [x[0] for x in mixcrRecall] + [x[0] for x in mixcrNoFilterRecall]
        samplePr["# of recalls"] = [x[1] for x in trust4Recall] + [x[1] for x in trust4FqRecall] + [x[1] for x in mixcrRecall] + [x[1] for x in mixcrNoFilterRecall]
        samplePr["Precision"] = trust4Precision + trust4FqPrecision + mixcrPrecision + mixcrNoFilterPrecision
        samplePr["Method"] = ["T4_BAM(%.2f)"%trust4Correlation] * len(topN) +                             ["T4_FQ(%.2f)"%trust4FqCorrelation] * len(topN) +                             ["MiXCR(%.2f)"%mixcrCorrelation] * len(topN) +                             ["MiXCR_0(%.2f)"%mixcrNoFilterCorrelation] * len(topN)
        #samplePr["Method"] = ["TRUST4"] * len(topN) + ["MiXCR"] * len(topN)
        samplePr["Sample"] = ["FZ-"+sample] * 4 * len(topN)
        samplePr["Correlation"] = [trust4Correlation] * len(topN) + [trust4FqCorrelation] * len(topN) + [mixcrCorrelation] * len(topN) + [mixcrNoFilterCorrelation] * len(topN)
        samplePr["Style"] = [0] * len(topN) + [1] * len(topN) + [0] * 2 * len(topN)
        ax = axes[axIdx]
        #hasLegend = False
        #if (axIdx == 5):
        #    hasLegend = "brief"
        snsFig = sns.lineplot(x="# of recalls", y="Precision", hue="Method", style="Method",                      data=samplePr, palette="bright", ax=ax, sizes=(3,3), size="Method" )
        ax.set_title( "FZ-"+str(sample), fontdict={'fontsize':"medium", "fontweight":"bold"} )
        ax.set(ylim=(0,100))
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles[1:], labels=labels[1:])
        if (axIdx > 0):
            snsFig.set(ylabel="") 
        pr = pr.append(samplePr)
        axIdx += 1
        
    #snsFig = sns.relplot(x="Sensitivity", y="Precision", hue="Sample", style="Method", \
    #                     kind="line", data=pr, palette="muted", height=10)
    #plt.suptitle( rdlen+"bp", fontsize="xx-large", fontweight="bold" )
    plt.savefig( "bcrval_PRcurve_wMiXCR0" + ".pdf", bbox_inches="tight", format="pdf" )


# In[25]:


display(trust4Recall)
display(trust4FqRecall)


# In[113]:


sns.set( font_scale=2)
fig, axes = plt.subplots( 1, 6, figsize=(30, 10))
axIdx = 0
for sample in ["20", "83", "94", "97", "116", "122"]:
    data = pr.loc[pr["Sample"]=="FZ-"+sample]
    ax = axes[axIdx]
    hasLegend = False
    if (axIdx == 5):
        hasLegend = "brief"
    snsFig = sns.lineplot(x="# of recall", y="Precision", hue="Method",                     data=data, palette="bright", ax=ax, sizes=(3,3), size="Method", legend=hasLegend )
    ax.set(ylim=(0,100))
    if (axIdx > 0):
        snsFig.set(ylabel="")
    axIdx += 1
    
#ax.set_xscale("log")


# Put six data sets in the same figure, 
for rdlen in ["150"]:
    data = pd.DataFrame()
    if (1):
        trust4Recall = []
        mixcrRecall = []
        trust4Precision = []
        mixcrPrecision = []
        for sample in ["20", "83", "94", "97", "116", "122"]: 
            irepList = irepListAll[rdlen][sample]
            trust4List = trust4ListAll[rdlen][sample]
            mixcrList = mixcrListAll[rdlen][sample]

            trust4Recall.append(GetSensitivity(irepList, trust4List)[1])
            mixcrRecall.append(GetSensitivity(irepList, mixcrList)[1])
            trust4Precision.append(GetPrecision(irepList, trust4List, False)[0])
            mixcrPrecision.append(GetPrecision(irepList, mixcrList, False)[0])
    data["Recall"] = trust4Recall + mixcrRecall
    data["Precision"] = trust4Precision + mixcrPrecision
    data["Method"] = ["TRUST4"] * len(trust4Recall) + ["MiXCR"] * len(mixcrRecall)
    data["Sample"] = ["FZ-20", "FZ-83", "FZ-94", "FZ-97", "FZ-116", "FZ-122"] * 2
    
    sns.set(font_scale=1.75)
    sns.set_style("white")
    fig = plt.figure(figsize=(7,7))
    plt.plot(trust4Recall, trust4Precision, "ro", markersize=9)
    plt.plot(mixcrRecall, mixcrPrecision, "bo", markersize=9)
    samples = ["FZ-20", "FZ-83", "FZ-94", "FZ-97", "FZ-116", "FZ-122"]
    for i in range(len(trust4Recall)):
        plt.plot([trust4Recall[i], mixcrRecall[i]], [trust4Precision[i], mixcrPrecision[i]], ls="--", c="g", )
        x = trust4Recall[i] * 0.9
        y = trust4Precision[i] + 0.5
        plt.text(x, y, s=samples[i], fontdict={"size":15})
    plt.xlim(0,8000)
    plt.ylim(20, 50)
    plt.xlabel("# of recalls")
    plt.ylabel("Precision(%)")
    
    trust4Patch = plt.plot([], [], marker="o", color="r", ms=9, ls="", label="TRUST4")[0]
    mixcrPatch = plt.plot([], [], marker="o", color="b", ms=9, ls="", label="MiXCR")[0]
    plt.legend(handles=[trust4Patch, mixcrPatch], loc="lower right", prop={"size":20})
    
    plt.savefig( "bcrval_"+rdlen+"bp_pr6samples" + ".pdf", bbox_inches="tight", format="pdf" )
