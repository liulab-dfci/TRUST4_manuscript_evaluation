#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import scipy as sp

from statsmodels.formula.api import ols

def GetChainType( v, j, c ):
    s = ""
    if ( j != "*" ):
        s = j 
    elif ( c != "*" ):
        s = c
    elif ( v != "*" ):
        s = v
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

def GetGeneType( s ):
    if ( s[3] == "V" ):
        return 0 
    elif ( s[3] == "J" ):
        return 2
    elif (s[3] == "D"):
        if (len(s) == 4):
            return 3
        else:
            return 1
    elif ( s[3] == "-" ):
        return -1
    else:
        return 3
    
chainTypeToName = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD" )
geneTypeToName = ("V", "J", "C" )

def GetMainGeneName( g ): # remove the allele type from the gene name
    fields = g.split( "*", 1 )
    return fields[0]

# Read in the gene coordination
fp = open( "bcrtcr.fa" )
geneCoord = {}
chainStrand = [1,1,1,1,1,1,1]
for line in fp:
    if ( line[0] != '>' ):
        continue
    line = line.rstrip()
    cols = line.split()
    #if ( cols[4] == '+' ):
    #    geneCoord[ cols[0][1:] ] = int( cols[2] ) 
    #else:
    #    geneCoord[ cols[0][1:] ] = -int( cols[2] ) 
    geneCoord[ cols[0][1:] ] = int( cols[2] ) 
    if ( cols[0][4] != 'V' and cols[0][4] != 'J' ):# use constant gene to determine chain's strand  
        chainType = GetChainType( cols[0][1:], cols[0][1:], cols[0][1:] )
        if ( cols[4] == "-" ):
            chainStrand[ chainType ] = -1
fp.close()

# Read in the sample information.
# b4055d85-4192-4a3e-87bc-cf7bf7fa13db	d0bd5182-d56f-4593-bdfa-bd791650d8f8_gdc_realn_rehead.bam	67ceb2f7d674453cbfddf7a196a5b520	8084198100	released	 TCGA-13-0804-01A-01R-1564-13	 TCGA-OV
fp = open( "gdc_manifest.kraken_info.txt" )
idToCancerType = {}
idStringToNum = {}
idNumToString = {}
cancerSamples = []
cancerSampleSet = set({})
normalSamples = []
normalSampleSet = set({})
sampleSurvival = []
k = 0 
for line in fp:
    line = line.rstrip("\n")
    cols = line.split("\t")
    for i in range(len(cols)):
        cols[i] = cols[i].strip()
    idToCancerType[ k ] = cols[6].split("-")[1] 
    idStringToNum[ cols[5] ] = k
    idNumToString[k] = cols[5]
    
    sampleType = int( cols[5].split( "-" )[3][0:2] )
    if ( sampleType < 10 ):
        cancerSamples.append( k )  
        cancerSampleSet.add(k)
    else:
        normalSamples.append(k)
        normalSampleSet.add(k)
    daysToBirth = -1
    daysToDeath = -1
    daysToFollowUp = -1
    if (len(cols) > 10):
        if (cols[7] != ""):
            daysToBirth = -int(cols[7])
        if (cols[8] != "" and cols[8] != "null"):
            daysToFollowUp = int(cols[8])
        if (cols[9] != ""):
            daysToDeath = int(cols[9])
        sampleSurvival.append([daysToBirth, daysToFollowUp, daysToDeath, cols[10]])
    else:
        sampleSurvival.append([-1, -1, -1, "NotReported"])
       
    k += 1
fp.close()

# Load the cancer subtype information
fp = open("annotation_of_TCGA_samples.txt")
idToSubCancerType = {}
header = fp.readline()
for line in fp:
    line = line.rstrip()
    cols = line.split("\t")
    if (cols[0] not in idStringToNum):
        continue
    idToSubCancerType[ idStringToNum[cols[0]] ] = cols[3]
fp.close()

cancerTypeToSamples = {}
subCancerTypeToSamples = {}
for sampleId in cancerSamples:
    cancerType = idToCancerType[ sampleId ] 
    if ( cancerType not in cancerTypeToSamples ):
        cancerTypeToSamples[ cancerType ] = [ sampleId ]
    else:
        cancerTypeToSamples[ cancerType ].append( sampleId )
    
    if (sampleId not in idToSubCancerType):
        continue
    subType = idToSubCancerType[sampleId] 
    if ( subType not in subCancerTypeToSamples ):
        subCancerTypeToSamples[ subType ] = [ sampleId ]
    else:
        subCancerTypeToSamples[ subType ].append( sampleId )    
        
# Get the tumor/normal pairs.
tmp = {}
for b in idStringToNum.keys():
    sampleType = int( b.split( "-" )[3][0:2] )
    if ( sampleType >= 10 ):
        patientBarcode = (b.split("-")[2], b.split("-")[-1])
        tmp[patientBarcode] = [-1, -1]
for b in idStringToNum.keys():
    sampleType = int( b.split( "-" )[3][0:2] )
    patientBarcode = (b.split("-")[2], b.split("-")[-1])
    if (patientBarcode not in tmp):
        continue
    if (sampleType < 10):
        tmp[patientBarcode][0] = idStringToNum[b]
    else:
        tmp[patientBarcode][1] = idStringToNum[b]
tumorNormalPair = {}
for b in tmp.keys():
    if (tmp[b][0] >= 0 and tmp[b][1] >= 0):
        tumorNormalPair[b] = tmp[b][:]

# Read in the tumor purity information
fp = open("tumor_purity.txt")
tmp = {}
sampleTumorPurity = []
header = fp.readline()
for line in fp:
    line = line.rstrip()
    cols = line.split("\t")
    tmp[cols[0]] = float(cols[1])
fp.close()

for i in range(len(idNumToString)):
    b = idNumToString[i]
    prefix = "-".join(b.split("-")[0:4])
    if (prefix not in tmp):
        sampleTumorPurity.append(-1)
    else:
        sampleTumorPurity.append(tmp[prefix])
tmp = {}

# Read in the deconvolution information.
sampleDeconv = {}
for cancer in cancerTypeToSamples:
    if (cancer == "LAML"):
        continue
    fp = open("ImmuneDeconv/TCGA-"+cancer)
    header = fp.readline()
    headerCols = header.rstrip().split()
    barcodePrefixToFull = {}
    
    for i in range(1, len(headerCols)):
        barcodePrefixToFull[headerCols[i]] = ""
    for sampleId in cancerTypeToSamples[cancer]:
        b = idNumToString[sampleId]
        #prefix = "-".join(b.split("-")[0:4])
        prefix = b[0:15]
        #print(prefix)
        if (prefix in barcodePrefixToFull):
            barcodePrefixToFull[prefix] = b
    
    for line in fp:
        cols = line.rstrip().split("\t")
        cellType = cols[0]
        for i in range(1,len(cols)):
            if (barcodePrefixToFull[headerCols[i]] != ""):
                sampleId = idStringToNum[ barcodePrefixToFull[headerCols[i]] ]
                if (cellType not in sampleDeconv):
                    sampleDeconv[cellType] = {}
                sampleDeconv[cellType][sampleId] = float(cols[i])
    fp.close()
    
# Read in the sequence depth.
fp = open("seqDepth/seqDepth.txt")
sampleSeqDepth = {}
for line in fp:
    cols = line.rstrip().split("\t")
    sampleSeqDepth[ idStringToNum[cols[0]] ] = float(cols[1])
fp.close()

   

# Plot the CDR3 similarity betwee/within samples
k = 0
withins = [0] * 21
for cancerType in ["COAD"]: #cancerTypeToSamples:
    #print(cancerType)
    fp = open("./within/TCGA-%s_similarity_within_histo.out"%cancerType)
    for line in fp:
        cols = line.rstrip().split()
        withins[int(cols[0])] += int(cols[2])
    fp.close()
    
withins = np.array(withins) / sum(withins)   

xv = []
yv = []
typev = []

for i in range(21):
    xv.append(i / 20)
    yv.append(withins[i] * 100)
    typev.append("Intra")

fp = open("./between/TCGA-COAD_similarity_between_hitso.out")
for line in fp:
        cols = line.rstrip().split()
        xv.append(float(cols[1]))
        yv.append(float(cols[3]))
        typev.append("Inter")
fp.close()

sns.set(font_scale=2, style="white")
fig = plt.subplots(figsize=(15, 10) )
df = pd.DataFrame()
df["Similarity"] = xv
df["Fraction of CDR3 pairs (%)"] = yv
df["Type"] = typev
snsFig = sns.barplot(x="Similarity", y="Fraction of CDR3 pairs (%)", data=df, hue="Type")
snsFig.set(title= "TCGA-COAD")
snsFig.legend()
for tick in snsFig.get_xticklabels():
    tick.set_rotation(45) 
plt.show()    


# In[5]:


import lifelines as ll
from lifelines.plotting import *

def GetSurvivalEvent(sampleId, limit = -1):
    if (sampleSurvival[sampleId][3] == "NotReported"):
        return -1, -1, 0
    event = 0 if (sampleSurvival[sampleId][3] == "Alive") else 1
    time = sampleSurvival[sampleId][1] if (sampleSurvival[sampleId][3] == "Alive")                 else sampleSurvival[sampleId][2]
    if (limit != -1 and time > limit):
        event = 0
        time = limit
    age = sampleSurvival[sampleId][0]
    return time, event, age

def ComputeEntropy(a):
    s = sum(a)
    if (len(a) == 1):
        return (-1, -1)
    fraction = np.array(a) / s
    return sp.stats.entropy(fraction, base=2), sp.stats.entropy(fraction, base=2) / math.log(len(a), 2)


# In[111]:


# Compute the entropy for each isotype
sampleIsotypeEntropy = {}

sampleIsotypeReadCnt = {}
fp = open("tcga_cluster_0.8.out")
prevSampleId = -1 
CDR3s = {}
clusters = {}
isotypeList = ["IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "*"]
for i in isotypeList:
    CDR3s[i] = {}
    clusters[i] = {}
    
for line in fp:
    cols = line.rstrip().split()
    sampleId = idStringToNum[cols[0]]
    if (sampleId != prevSampleId and prevSampleId != -1):
        #if (prevSampleId == 2394):
        #    break
        sampleIsotypeEntropy[prevSampleId] = {}
        for i in isotypeList:
            if (prevSampleId not in sampleIsotypeReadCnt
                or i not in sampleIsotypeReadCnt[prevSampleId]):
                continue
            #print(CDR3s)
            sampleIsotypeEntropy[prevSampleId][i] =             {
                "ByCDR3nt": ComputeEntropy(list(CDR3s[i].values())),
                "ByCluster": ComputeEntropy(list(clusters[i].values()))
            }
        CDR3s = {}
        clusters = {}  
        for i in isotypeList:
            CDR3s[i] = {}
            clusters[i] = {}
    
    chainType = GetChainType(cols[3], cols[5], cols[6])
    if (chainType != 0 ):
        prevSampleId = sampleId
        continue
    isotype = "*"
    if (cols[6] != "*"):
        isotype = GetMainGeneName(cols[6])
    if (sampleId not in sampleIsotypeReadCnt):
        sampleIsotypeReadCnt[sampleId] = {}
    if (isotype not in sampleIsotypeReadCnt[sampleId]):
        sampleIsotypeReadCnt[sampleId][isotype] = 0
    sampleIsotypeReadCnt[sampleId][isotype] += float(cols[11])
    if (cols[9] not in CDR3s[isotype]):
        CDR3s[isotype][cols[9]] = 0
    CDR3s[isotype][cols[9]] += float(cols[11])
    if (cols[1] not in clusters[isotype]):
        clusters[isotype][cols[1]] = 0
    clusters[isotype][cols[1]] += float(cols[11])
    prevSampleId = sampleId
fp.close()

sampleIsotypeEntropy[prevSampleId] = {}
for i in isotypeList:
    if (i not in sampleIsotypeReadCnt[prevSampleId]):
        continue
    sampleIsotypeEntropy[prevSampleId][i] =     {
        "ByCDR3nt": ComputeEntropy(list(CDR3s[i].values())),
        "ByCluster": ComputeEntropy(list(clusters[i].values()))
    }


# In[7]:


# Sample isotype read count
fp = open("tcga_simpleRepNoPartial.out")
sampleIsotypeReadCnt = {}
for line in fp: 
    line = line.rstrip()
    cols = line.split()
    chainType = GetChainType(cols[5], cols[7], cols[8])
    if (chainType != 0):
        continue
    count = float( cols[1] ) 
    sampleId = idStringToNum[ cols[0] ]
    isotype = "*"
    if (cols[8] != "*"):
        isotype = GetMainGeneName(cols[8])
    if (sampleId not in sampleIsotypeReadCnt):
        sampleIsotypeReadCnt[sampleId] = {}
    if (isotype not in sampleIsotypeReadCnt[sampleId]):
        sampleIsotypeReadCnt[sampleId][isotype] = 0
    sampleIsotypeReadCnt[sampleId][isotype] += count
fp.close()


# In[285]:


# Obtain the gene expression for each sample
sampleTPM = {}
fp = open("TPM_of_TCGA.trbighc.txt")
header = fp.readline()
colIdxToSampleId = []
cols = header.rstrip().split()
for c in cols:
    if (c not in idStringToNum):
        colIdxToSampleId.append(-1)
        sampleTPM[-1] = []
    else:
        colIdxToSampleId.append( idStringToNum[c] )
        sampleTPM[idStringToNum[c]] = {}
        
geneNameToId = {}
geneIdToName = []
for line in fp:
    cols = line.rstrip().split()
    geneName = cols[0]
    for i in range(1, len(cols)):
        sampleId = colIdxToSampleId[i - 1]
        if (sampleId == -1):
            continue
        #if (geneName[0] == "T"):
        #    sampleTPM[sampleId][0] += float(cols[i])
        #else:
        sampleTPM[sampleId][geneName] = float(cols[i])
    geneNameToId[geneName] = len(geneIdToName)
    geneIdToName.append(geneName)
fp.close()


# In[286]:


sampleTPM[0]


# In[314]:


# Plot isotype fraction and 
entropyMetric = "ByCluster"
entropyType = 1
cancerType = "COAD"

isotypeList = ["IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "*"]
isotypeToId = {}
chosenIsotypes = ["IGHG1", "IGHA1"]
for i in range(len(isotypeList)):
    isotypeToId[isotypeList[i]] = i

chosenSamples = {}
for metric in ["Fraction", "Clonality"]:
    for isotype in chosenIsotypes:
        i = 0
        infos = []
        if (metric == "Clonality"):
            for sampleId in cancerTypeToSamples[cancerType]:
                time, event, age = GetSurvivalEvent(sampleId)
                if (time == -1 or age == -1): #or sampleTumorPurity[sampleId] == -1):
                    continue
                if (sampleId not in sampleIsotypeEntropy):
                    continue
                #if (sampleReadCnt[sampleId][chainType] < 10):
                #    continue
                if (sampleId not in sampleIsotypeEntropy
                   or isotype not in sampleIsotypeEntropy[sampleId]
                   or entropyMetric not in sampleIsotypeEntropy[sampleId][isotype]
                   or sampleIsotypeEntropy[sampleId][isotype][entropyMetric][entropyType] == -1):
                    continue
                    
                if (sampleId not in chosenSamples):
                    chosenSamples[sampleId] = 0 
                #chosenSamples[sampleId] += 1
        elif (metric == "Fraction"):
             for sampleId in cancerTypeToSamples[cancerType]:
                isotypeCount = [0] * len(isotypeList)
                time, event, age = GetSurvivalEvent(sampleId)
                if (time == -1 or age == -1): #or sampleTumorPurity[sampleId] == -1):
                    continue
                #if (sampleId not in sampleIsotypeReadCnt):
                #    continue
                #for i in sampleIsotypeReadCnt[sampleId]:
                #    isotypeCount[ isotypeToId[i] ] += sampleIsotypeReadCnt[sampleId][i]
                #if (sum(isotypeCount) - isotypeCount[9] == 0):
                #    continue
                if (sampleId not in sampleTPM or sum(sampleTPM[sampleId].values()) == 0):
                    continue
                #print(isotypeCount)
                #tmp = isotypeCount[isotypeToId[chosenIsotype]] / sum(isotypeCount)
                #infos.append([isotypeCount[isotypeToId[chosenIsotype]] / sum(isotypeCount), 
                #              time, event, age])
                if (sampleId not in chosenSamples):
                    chosenSamples[sampleId] = 0 
                chosenSamples[sampleId] += 1
                
                    
                
sns.set(font_scale=1.25, style="white")
axIdx = 0 
figWidth = 15
figHeight = 5 * len(chosenIsotypes)
fig, axes = plt.subplots(len(chosenIsotypes), 3, figsize=(figWidth, figHeight ) ) 
for metric in ["Abundance", "Fraction", "Clonality"]:
    for isotype in chosenIsotypes:
        i = 0
        infos = []
        if (metric == "Clonality"):
            for sampleId in cancerTypeToSamples[cancerType]:
                time, event, age = GetSurvivalEvent(sampleId)
                if (time == -1 or age == -1): #or sampleTumorPurity[sampleId] == -1):
                    continue
                if (sampleId not in sampleIsotypeEntropy):
                    continue
                #if (sampleReadCnt[sampleId][chainType] < 10):
                #    continue
                entropy = 1
                if (sampleId not in sampleIsotypeEntropy
                   or isotype not in sampleIsotypeEntropy[sampleId]
                   or entropyMetric not in sampleIsotypeEntropy[sampleId][isotype]
                   or sampleIsotypeEntropy[sampleId][isotype][entropyMetric][entropyType] == -1):
                    entropy = 2 # 1-entropy=-1
                else:
                    entropy = sampleIsotypeEntropy[sampleId][isotype][entropyMetric][entropyType]
                if (sampleId not in chosenSamples or chosenSamples[sampleId] != len(chosenIsotypes)):
                    continue
                #print(sampleId, sampleEntropy[sampleId][chainType] )
                infos.append([1 - entropy, 
                              time, event, age])#, sampleTumorPurity[sampleId]])
        elif (metric == "Fraction"):
            for sampleId in cancerTypeToSamples[cancerType]:
                #isotypeCount = [0] * len(isotypeList)
                time, event, age = GetSurvivalEvent(sampleId)
                if (time == -1 or age == -1): #or sampleTumorPurity[sampleId] == -1):
                     continue
                if (sampleId not in sampleIsotypeReadCnt):
                    continue
                #for i in sampleIsotypeReadCnt[sampleId]:
                #    isotypeCount[ isotypeToId[i] ] += sampleIsotypeReadCnt[sampleId][i]
                #if (sum(isotypeCount) - isotypeCount[9] == 0):
                #    continue
                #print(isotypeCount)
                #tmp = isotypeCount[isotypeToId[chosenIsotype]] / sum(isotypeCount)
                #infos.append([isotypeCount[isotypeToId[chosenIsotype]] / sum(isotypeCount), 
                #              time, event, age])
                if (sampleId not in chosenSamples or chosenSamples[sampleId] != len(chosenIsotypes)):
                    continue
                tmp = sampleTPM[sampleId][isotype] / (sum(sampleTPM[sampleId].values()) - sampleTPM[sampleId]["TRBC2"])
                #tmp = (isotypeCount[isotypeToId[isotype]] ) / \
                #    (sum(isotypeCount) - isotypeCount[9])
                infos.append([tmp, time, event, age])  
        elif (metric == "Abundance"):
            for sampleId in cancerTypeToSamples[cancerType]:
                time, event, age = GetSurvivalEvent(sampleId)
                if (time == -1 or age == -1): #or sampleTumorPurity[sampleId] == -1):
                     continue
                #if (sampleId not in sampleIsotypeReadCnt):
                #    continue
                if (sampleId not in chosenSamples or chosenSamples[sampleId] != len(chosenIsotypes)):
                    continue
                tmp = sampleTPM[sampleId][isotype]
                infos.append([tmp, time, event, age])      
        if (len(infos) == 0):
            continue
            
        if (metric == "Clonality"): # used median to represent the p-value
            x = []
            for i in infos:
                if (i[0] != -1):
                    x.append(i[0])
            m = np.median(x)
            for i in range(len(infos)):
                if (infos[i][0] == -1):
                    infos[i][0] = m
        df = pd.DataFrame(infos, columns=[metric, "time", "event", "age"])#, "purity"])
        lowCutoff = np.quantile(df[metric], 0.5)
        highCutoff = np.quantile(df[metric], 0.5)
        highDf = df[ df[metric] >= highCutoff]
        lowDf = df[ df[metric] < lowCutoff]
        
        ax = axes[axIdx%len(chosenIsotypes)][axIdx//len(chosenIsotypes)]
        kmf = ll.KaplanMeierFitter()
        label = isotype + " " + metric
        f1 = kmf.fit(highDf["time"], highDf["event"], label="High " + label)
        ax = kmf.plot(ax=ax, ci_show=False)
        kmf2 = ll.KaplanMeierFitter()
        f2 = kmf2.fit(lowDf["time"], lowDf["event"], label="Low " + label)
        ax = kmf2.plot(ax = ax, ci_show=False)
        #add_at_risk_counts(f1, f2, ax = ax)
        print(len(lowDf), len(highDf))
        cph = ll.CoxPHFitter()
        cph.fit(df, "time", "event")
        print(cph.summary)
        p = cph.summary.loc[metric]["p"]
        hr = cph.summary.loc[metric]["exp(coef)"]
        if (p > 0.01):
            ax.text(0, 0, "cox (hr=%.3f, p=%.3f)"%(hr, p))
        else:
            ax.text(0, 0, "cox (hr=%.3f, p=%.2e)"%(hr, p))
        ax.set(ylim=(-0.05, 1.05))
        axIdx += 1
plt.tight_layout()
plt.show()
