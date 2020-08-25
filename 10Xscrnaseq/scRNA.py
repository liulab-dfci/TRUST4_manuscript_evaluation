#!/usr/bin/env python
# coding: utf-8

# In[159]:


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
#import matplotlib.style as mpl
import string
import re

#mplDefaultParams = mpl.rcParams


# In[160]:


import sys
sys.executable


# In[161]:


def GetChainType( v, j, c ):
    s = ""
    if (v == '-'):
        v = "*"
    if (j == '-'):
        j = "*"
    if (c == '-'):
        c = "*"
    if ( j != "*" ):
        s = j 
    elif ( c != "*" and c != "-"):
        s = c
    elif ( v != "*" and v != "-"):
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

chainTypeToName = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD" )

def GetMainGeneName( g ): # remove the allele type from the gene name
    fields = g.split( "*", 1 )
    ret = fields[0]
    for c in ["IGLC", "IGKC", "TRAC", "TRBC", "TRGC", "TRDC"]:
        if (c in ret):
            ret = ret[0:4]
            return ret
    return ret


# In[162]:


#dataset = "10X_BMMC_PBMC"
#dataset = "10X_PBMC_8K"
dataset = "10X_PBMC5K_protein"
#dataset = "10X_BMMC_CLL"
#dataset = "10X_vdj_v1_hs_pbmc3"
#dataset = "10X_vdj_nextgem_hs_pbmc3"
#dataset = "10X_vdj_v1_hs_nsclc"
#dataset = "10X_vdj_nextgem_hs_pbmc3_umi"

SeuratFile = ""
umapFile = ""
trustBarcodeReportFile = ""
trustReportFile = ""

path="/Users/lsong/Song/MyProjects/TRUST4/scRNA/" + dataset + "/"
elif (dataset == "10X_PBMC5K_protein"):
    SeuratFile = path+"10X_PBMC5K_protein_Seurat_metadata.tsv"
    umapFile = path+"10X_PBMC5K_protein_Seurat_umap.tsv"
    trustBarcodeReportFile = path+"TRUST_5k_pbmc_protein_v3_nextgem_possorted_genome_bam_barcode_report.tsv"
    trustReportFile = path+"TRUST_5k_pbmc_protein_v3_nextgem_possorted_genome_bam_report.tsv"
elif (dataset == "10X_BMMC_CLL"):
    SeuratFile = path+"10X_BMMC_CLL_Seurat_metadata.tsv"
    umapFile = path+"10X_BMMC_CLL_Seurat_umap.tsv"
    trustBarcodeReportFile = path+"TRUST_scRNA_BMMC_CLL_donor_mRNA_barcode_report.tsv"
    trustReportFile = path+"TRUST_scRNA_BMMC_CLL_donor_mRNA_report.tsv"
elif (dataset == "10X_vdj_v1_hs_pbmc3"):
    SeuratFile = path+"10X_vdj_v1_pbmc3_Seurat_metadata.tsv"
    umapFile = path+"10X_vdj_v1_pbmc3_Seurat_umap.tsv"
    trustBarcodeReportFile = path+"TRUST_vdj_v1_hs_pbmc3_possorted_genome_bam_barcode_report.tsv"
    #trustBarcodeReportFile = path+"tmp.tsv"
    trustReportFile = path+"TRUST_vdj_v1_hs_pbmc3_possorted_genome_bam_report.tsv"
elif (dataset == "10X_vdj_nextgem_hs_pbmc3"):
    SeuratFile = path+"10X_vdj_nextgem_hs_pbmc3_Seurat_metadata.tsv"
    umapFile = path+"10X_vdj_nextgem_hs_pbmc3_Seurat_umap.tsv"
    trustBarcodeReportFile = path+"TRUST_vdj_nextgem_hs_pbmc3_possorted_genome_bam_barcode_report.tsv"
    trustReportFile = path+"TRUST_vdj_nextgem_hs_pbmc3_possorted_genome_bam_report.tsv"
elif (dataset == "10X_vdj_v1_hs_nsclc"):
    SeuratFile = path+"10X_vdj_v1_hs_nsclc_Seurat_metadata.tsv"
    umapFile = path+"10X_vdj_v1_hs_nsclc_Seurat_umap.tsv"
    trustBarcodeReportFile = path+"10X_vdj_v1_hs_nsclc_5gex_barcode_report.tsv"
    trustReportFile = path+"10X_vdj_v1_hs_nsclc_5gex_report.tsv"  


# In[ ]:





# In[163]:


def IsProductiveAa(cdr3aa):
    if ("_" in cdr3aa or cdr3aa == "*" or cdr3aa == "None"):
        return False
    return True

def IsProductiveCell(l):
    for i in [2, 3]:
        if (l[i] == "*"):
            continue 
        if (IsProductiveAa(l[i].split(",")[5])):
            return 1
    return 0
        

def GetChainState(barcodeInfo, considerProductive):
    chainState = {}
    for b in barcodeInfo.keys():
        clusterName = barcodeInfo[b][1]
        chain1 = 1 
        if ( barcodeInfo[b][2] == "*" ):
            chain1 = 0 
        chain2 = 1
        if ( barcodeInfo[b][3] == "*" ):
            chain2 = 0 
            
        productive = IsProductiveCell(barcodeInfo[b])
        
        #hasChainStatus = productive | (chain1 << 1) | (chain2 << 2) 
        hasChainStatus = chain1 | (chain2 << 1)
        if (considerProductive):
            hasChainStatus = productive | (chain1 << 1) | (chain2 << 2) 
            
        if (clusterName in chainState):
            chainState[clusterName][0] += 1
            chainState[clusterName][ hasChainStatus ] += 1
        else:
            chainState[clusterName] = [0, 0, 0, 0]
            if (considerProductive):
                chainState[clusterName] = [0, 0, 0, 0, 0, 0, 0, 0]
            chainState[clusterName][0] = 1
            chainState[clusterName][ hasChainStatus ] = 1
    data = pd.DataFrame.from_dict(chainState, orient="index")
    data = data.reindex(sorted(data.columns), axis=1)
    if (considerProductive):
        data = data.drop(columns=[1])
    return data

def CellTypeRank(t):
    if ("T" in t):
        return (0, t)
    elif ("B" in t):
        return (1, t)
    elif ("ther" in t):
        return (3, t)
    elif (t == "Sum"):
        return (4, t)
    elif (t == "Missing"):
        return (5, t)
    else:
        return (2, t)

def GetCellTypeComposition(barcodeSeurat, barcodeInfo):
    # Get the total number of Seurat cells for each cell types.
    cellTypes = {}
    for b in barcodeSeurat.keys():
        clusterName = barcodeSeurat[b][2]
        if (clusterName in cellTypes):
            cellTypes[clusterName] += 1
        else:
            cellTypes[clusterName] = 1
            
    infoCellTypes = {"B":0, "abT":0, "gdT":0}
    for b in barcodeInfo.keys():
        infoCellTypes[barcodeInfo[b][1]] += 1
        
    cellTypeComposition = {}
    cellTypeComposition["Sum"] = {'#_of_cells': len(barcodeSeurat),'B':0, 'abT':0, 'gdT':0 }
    cellTypeComposition["Missing"] = {'#_of_cells': 0,'B':0, 'abT':0, 'gdT':0 }
    
    for b in barcodeSeurat.keys():
        clusterName = barcodeSeurat[b][2]
        if (b not in barcodeInfo):
            continue
        #if (clusterName == "CD8 T" and barcodeTrust[b][1] == 'B' ):
        #    print(b)
        if ( clusterName not in cellTypeComposition ):
            cellTypeComposition[clusterName] = {'#_of_cells': cellTypes[clusterName],'B':0, 'abT':0, 'gdT':0 }
        cellTypeComposition[clusterName][ barcodeInfo[b][1] ] += 1
        cellTypeComposition["Sum"][ barcodeInfo[b][1] ] += 1
    #print(cellTypeComposition)
    for t in ["B", "abT", "gdT"]:
        cellTypeComposition["Missing"][t] = infoCellTypes[t] - cellTypeComposition["Sum"][t] 
    cnt = 0
    for b in barcodeInfo:
        if (b not in barcodeSeurat):
            cnt += 1
    data = pd.DataFrame.from_dict(cellTypeComposition)
    data = data.reindex(sorted(data.columns, key=CellTypeRank), axis=1)
    return data, cnt
        


# In[164]:


barcodeSeurat = {}
fp = open( SeuratFile, "r" )
fp.readline()
for line in fp:
    line = line.rstrip()
    cols = line.split("\t")
    if (len(cols) > 7 and cols[7] != "RNA"): # Seurat's cocluster
        continue
    barcodeSeurat[ cols[0].split("-")[0] ] = [int(cols[2]), int(cols[3]), cols[6]]

fp = open(umapFile, "r")
fp.readline()
for line in fp:
    line = line.rstrip()
    cols = line.split("\t")
    barcodeSeurat[cols[0].split("-")[0]] += [float(cols[1]), float(cols[2])]
fp.close()


# In[165]:


barcodeTrust = {}
fp = open( trustBarcodeReportFile, "r" )
fp.readline()
for line in fp:
    line = line.rstrip()
    cols = line.split("\t")
    barcodeTrust[ cols[0].split("-")[0] ] = cols
fp.close()
    


# In[166]:


# Summary on the number of cell, and the number of cells in different clusters.
print( len(barcodeSeurat) )
cellTypes = {}
for b in barcodeSeurat.keys():
    clusterName = barcodeSeurat[b][2]
    if (clusterName in cellTypes):
        cellTypes[clusterName] += 1
    else:
        cellTypes[clusterName] = 1
#print(cellTypes)
#pd.DataFrame.from_dict(cellTypes, orient="index")


# In[167]:


# Summary on the number of cell, and the number of cells in different clusters.
print( len(barcodeTrust) )
display( GetChainState(barcodeTrust, False))


# In[ ]:





# In[168]:


# Look at the reads/cell v.s. the chainState
fig, axes = plt.subplots(1, 2)
axIdx = 0 
for cellType in ["B", "abT"]:
    barcodeChainState = {}
    for b in barcodeTrust.keys():
        if (b not in barcodeSeurat):
            continue
        if (barcodeTrust[b][1] != cellType):
            continue
        chain1 = 1 
        if ( barcodeTrust[b][2] == "*" ):
            chain1 = 0 
        chain2 = 1
        if ( barcodeTrust[b][3] == "*" ):
            chain2 = 0 

        hasChainStatus = (chain1 | (chain2 << 1))  
        barcodeChainState[b] = [2 if hasChainStatus == 3 else 1, barcodeSeurat[b][0]]
    data = pd.DataFrame()
    if (len(barcodeChainState)>0):
        data = pd.DataFrame.from_dict(barcodeChainState, orient="index")
    data.columns = ["Chains", "Reads"]
    #print(data.loc[(data["Chains"]==2)])
    ax = axes[axIdx]
    snsFig = sns.boxplot(x="Chains", y="Reads", data = data, ax=ax, palette="muted")
    ax.set(ylim=(0, 15000))
    if (axIdx > 0):
        snsFig.set(ylabel="")
        ax.get_yaxis().set_visible(False)
    else:
        snsFig.set(ylabel="UMIs")
    snsFig.set(xlabel="Chains (%s)"%cellType )
            
    axIdx += 1


# In[ ]:





# In[169]:


cellTypeCompositionDf, notInSeurat = GetCellTypeComposition(barcodeSeurat, barcodeTrust)
print(notInSeurat)
display(cellTypeCompositionDf)

for chainCount in [1, 2]:
    barcodeTest = {}
    for b in barcodeTrust:    
        cnt = 0
        if (barcodeTrust[b][2] != "*"):
            cnt += 1
        if (barcodeTrust[b][3] != "*"):
            cnt += 1
        if (cnt != chainCount):
            continue 
        #if (cnt == 2):
        #    if (b not in barcodeSeurat):
        #        print(b)
        barcodeTest[b] = barcodeTrust[b][:]
    cellTypeCompositionDf, cnt = GetCellTypeComposition(barcodeSeurat, barcodeTest)
    print(cnt)
    display(cellTypeCompositionDf)
    
    


# In[170]:


# Calculate general repertoire stats for each cell type
#for cellType in ["CD8 T", "CD4 Memory T", "CD4 Naive T"]:
#    for 


# In[171]:


# umap the results
Seurat = pd.DataFrame.from_dict(barcodeSeurat, orient="index")
Seurat = Seurat.reset_index()
Seurat.columns = ["barcode", "nCount", "nFeature", "cluster", "UMAP_1", "UMAP_2"]

trust4Cell = []
for b in barcodeSeurat.keys():
    if (b in barcodeTrust):
        clusterName = barcodeTrust[b][1]  
    else:
        clusterName = "NA"
    trust4Cell.append(clusterName)

Seurat["TRUST4"] = trust4Cell

trust4Cell = []
for b in barcodeSeurat.keys():
    cGene = "*"
    if (b in barcodeTrust and barcodeTrust[b][2] != "*"):
        cGene = barcodeTrust[b][2].split(",")[3]
    trust4Cell.append(cGene)
Seurat["C gene"] = trust4Cell

#snsFig = sns.relplot(x="UMP_1", y="UMP_2", hue="TRUST4",
#                     data=Seurat.loc[(Seurat["cluster"]=="CD8 T")], alpha=0.75, palette="bright" )

snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue="cluster", data = Seurat, alpha=0.75, palette="bright" )


# In[172]:


name=""
if (dataset == "10X_BMMC_PBMC"):
    name="CD8 T"
else:
    name = "CD8Tcells"
snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue="TRUST4",
                    data=Seurat.loc[(Seurat["cluster"]==name) & (Seurat["TRUST4"]!="NA")], 
                     alpha=0.75, palette="bright", hue_order=["B", "abT", "gdT"] )


# In[173]:


# my own way to visualize both Seurat and TRUST4
sns.set(font_scale=1.5, style="white")
fig = plt.figure(figsize=(10, 10))
plt.grid(False)
#plt.set_facecolor([1, 1, 1])
cellTypeToScale = {}
k = 0
for i in Seurat["cluster"].unique():
    cellTypeToScale[i] = 0.4 + k / (len(Seurat["cluster"].unique()) - 1) * 0.45
    k += 1

for cellType in cellTypeToScale:
    SeuratCellType = Seurat.loc[ Seurat["cluster"] == cellType ]
    plt.plot(SeuratCellType["UMAP_1"], SeuratCellType["UMAP_2"], "o", marker="o", ms=1,
            color=(cellTypeToScale[cellType], cellTypeToScale[cellType], cellTypeToScale[cellType]))


# Plot TRUST4's result
trustCellTypeToColor = {"B":"r", "abT":"cyan", "gdT":"g"}
for trustCellType in trustCellTypeToColor:
    SeuratTrustCellType = Seurat.loc[ Seurat["TRUST4"] == trustCellType]
    size = 1
    if (trustCellType == "gdT"):
        size=3
    if ("vdj" not in dataset and trustCellType == "abT"):
        size = 3
    plt.plot(SeuratTrustCellType["UMAP_1"], SeuratTrustCellType["UMAP_2"], "o", marker="o", ms=size,
            color = trustCellTypeToColor[trustCellType])

# Add the annotation
for cellType in cellTypeToScale:    
    SeuratCellType = Seurat.loc[ Seurat["cluster"] == cellType ]
    x = np.median(SeuratCellType["UMAP_1"]) - 1
    y = np.median(SeuratCellType["UMAP_2"])
    plt.text(x, y, cellType, fontdict={"size":20, "weight":"bold"})
plt.xlabel("UMAP_1")
plt.ylabel("UMAP_2")

# Add legend
patches = []
for trustCellType in trustCellTypeToColor:
    size = 5
    if (trustCellType == "gdT"):
        size = 10
    if ("vdj" not in dataset and trustCellType == "abT"):
        size = 10
    patches.append(plt.plot([], [], marker="o", color=trustCellTypeToColor[trustCellType], 
                            ms=size, ls="", label=trustCellType)[0])
plt.legend(handles=patches, loc="lower left")

#plt.savefig( dataset+"_"+"Seurat_trust4_umap" + ".pdf", bbox_inches="tight", format="pdf" )


# In[174]:


SeuratTrustCellType = Seurat.loc[ (Seurat["TRUST4"] == "gdT") ]
SeuratTrustCellType["barcode"]


# In[175]:


# Color the isotypes in B cell cluster.
snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue="C gene", style = "cluster",
                    data=Seurat.loc[(Seurat["cluster"].str.contains("B") )
                                      & (Seurat["C gene"].str.contains("IGH"))],
                    alpha=0.75, palette="bright")


# In[176]:


a=Seurat.loc[Seurat["cluster"].str.contains("B")]


# In[177]:


# Color the isotypes in B cell cluster.
snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue="C gene",
                    data=Seurat.loc[Seurat["cluster"].str.contains("Plasma") 
                                      & (Seurat["C gene"].str.contains("IGH"))],
                    alpha=0.75, palette="bright")


# In[178]:


# Explore the distribution of the cells with same CDR3.
cdr3Cnt = [{}, {}, {}, {}, {}, {}, {}]
for b in barcodeTrust.keys():
    for i in [2, 3]:
        if (barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue
        if (barcodeTrust[b][i] == "*"):
            continue    
        cols = barcodeTrust[b][i].split(",")
        chainType = GetChainType(cols[0], cols[2], cols[3])
        if (cols[4] not in cdr3Cnt[chainType]):
            cdr3Cnt[chainType][cols[4]] = 0  
        cdr3Cnt[chainType][cols[4]] += 1

for chainType in range(1):
    maxCdr3 = ""
    for cdr3 in cdr3Cnt[chainType]:
        if (maxCdr3 == "" or cdr3Cnt[chainType][cdr3] > cdr3Cnt[chainType][maxCdr3]):
            maxCdr3 = cdr3 
    print(cdr3Cnt[chainType][maxCdr3])
    data = Seurat.copy()
    chosenBarcode = set({})
    for b in barcodeTrust.keys():
        idx = 3 
        if (chainType == 0 or chainType == 3):
            idx = 2 
        if (barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue    
        if (barcodeTrust[b][idx] == "*"):
            continue
        cols = barcodeTrust[b][idx].split(",")
        if (cols[4] == maxCdr3):
            chosenBarcode.add(b)
    mark = []
    for b in barcodeSeurat.keys():
        if (b in chosenBarcode):
            mark.append(1) 
        else:
            mark.append(0)
    data[ maxCdr3 ] = mark
    display(data.loc[(data[maxCdr3]==1)])
    snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue=maxCdr3, sizes=maxCdr3,
                    data=data.loc[data["cluster"].str.contains("B")],
                    alpha=0.75, palette="bright")
    #snsFig.set(xlim=(-14,-10))
    #snsFig.set(ylim=(-4, -1))


# In[179]:


for chainType in range(1):
    maxCdr3 = "TGTGCGAGGAAAGGGGGTGATTTCTTTGGTTTTGATATCTGG"
    data = Seurat.copy()
    chosenBarcode = set({})
    for b in barcodeTrust.keys():
        idx = 3 
        if (chainType == 0 or chainType == 3):
            idx = 2 
        if (barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue    
        if (barcodeTrust[b][idx] == "*"):
            continue
        cols = barcodeTrust[b][idx].split(",")
        if (cols[4] == maxCdr3):
            chosenBarcode.add(b)
    mark = []
    for b in barcodeSeurat.keys():
        if (b in chosenBarcode):
            mark.append(1) 
        else:
            mark.append(0)
    data[ maxCdr3 ] = mark
    display(data.loc[(data[maxCdr3]==1)])
    snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue=maxCdr3, sizes=maxCdr3,
                    data=data.loc[data["cluster"].str.contains("B")],
                    alpha=0.75, palette="bright")
    #snsFig.set(xlim=(-14,-10))
    #snsFig.set(ylim=(-4, -1))


# In[180]:


# Explore plasma cells with the same CDR3s.
cdr3Cnt = [{}, {}, {}, {}, {}, {}, {}]
for b in barcodeTrust.keys():
    for i in [2, 3]:
        if (barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue
        if (barcodeTrust[b][i] == "*"):
            continue    
        cols = barcodeTrust[b][i].split(",")
        chainType = GetChainType(cols[0], cols[2], cols[3])
        if (cols[4] not in cdr3Cnt[chainType]):
            cdr3Cnt[chainType][cols[4]] = 0  
        cdr3Cnt[chainType][cols[4]] += 1

for cdr3 in cdr3Cnt[0]:
    data = Seurat.copy()
    chosenBarcode = set({})
    for b in barcodeTrust.keys():
        if (barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue    
        if (barcodeTrust[b][idx] == "*"):
            continue
        cols = barcodeTrust[b][idx].split(",")
        if (cols[4] == cdr3):
            chosenBarcode.add(b)
    mark = []
    for b in barcodeSeurat.keys():
        if (b in chosenBarcode and barcodeSeurat[b][2] == "PlasmaCells"):
            mark.append(1) 
        else:
            mark.append(0)
    data[ cdr3 ] = mark
    if (sum(mark)>1):
        display(data.loc[(data[cdr3]==1)])
#snsFig.set(xlim=(-14,-10))
#snsFig.set(ylim=(-4, -1))


# In[181]:


barcodeSeurat


# In[182]:


for chainType in [1]:
    maxCdr3 = "TGTCAACAGTATGGTAACTCATTGTGGACGTTC"
    data = Seurat.copy()
    chosenBarcode = set({})
    for b in barcodeTrust.keys():
        idx = 3 
        if (chainType == 0 or chainType == 3):
            idx = 2 
        if (barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue 
        if (barcodeTrust[b][idx] == "*"):
            continue
        cols = barcodeTrust[b][idx].split(",")
        if (cols[4] == maxCdr3):
            chosenBarcode.add(b)
    mark = []
    for b in barcodeSeurat.keys():
        if (b in chosenBarcode):
            mark.append(1) 
        else:
            mark.append(0)
    data[ maxCdr3 ] = mark
    display(data.loc[(data[maxCdr3]==1)])
    snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue=maxCdr3, sizes=maxCdr3,
                    data=data.loc[data["cluster"].str.contains("B")],
                    alpha=0.75, palette="bright")


# In[183]:


# Compute the eveness(clonality) for B cells' heavy chain, light chain.
for i in [2, 3]:
    cdr3Cnt = {}
    totalCnt = 0 
    for b in barcodeTrust:
        if (barcodeTrust[b][1] != "B" or barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue
        if (b not in barcodeSeurat or 
                ("Plasma" not in barcodeSeurat[b][2])):
            continue
        #if (b not in barcode10X or barcode10X[b][1] != "B"):
        #    continue
        cols = barcodeTrust[b][i].split(",") 
        if (cols[4] not in cdr3Cnt):
            cdr3Cnt[cols[4]] = 0
        cdr3Cnt[cols[4]] += 1
        totalCnt += 1
    for c in cdr3Cnt:
        cdr3Cnt[c] /= totalCnt 
    entropy = 0
    for c in cdr3Cnt:
        entropy += -cdr3Cnt[c] * math.log(cdr3Cnt[c]) 
    print(entropy, len(cdr3Cnt), totalCnt)
    
for i in [2, 3]:
    break
    cdr3Cnt = {}
    totalCnt = 0 
    for b in barcode10X:
        if (barcode10X[b][1] != "B" or barcode10X[b][2] == "*" or barcode10X[b][3] == "*"):
            continue
        if (b not in barcodeSeurat or 
                ("Plasma" not in barcodeSeurat[b][2])):
            continue
        cols = barcode10X[b][i].split(",") 
        if (cols[4] not in cdr3Cnt):
            cdr3Cnt[cols[4]] = 0
        cdr3Cnt[cols[4]] += 1
        totalCnt += 1
    for c in cdr3Cnt:
        cdr3Cnt[c] /= totalCnt 
    entropy = 0
    for c in cdr3Cnt:
        entropy += -cdr3Cnt[c] * math.log(cdr3Cnt[c]) 
    print(entropy, len(cdr3Cnt), totalCnt)
        


# In[ ]:





# In[ ]:





# In[184]:


# Only analyze the single-cell data with VDJ
barcode10X = {}
barcodeChainAbund = {}
barcodeCellType ={}
barcode10Xbt = {}
if ("vdj" in dataset):
    for i in ["b", "t"]:
        #fp = open(path + dataset[4:] + "_" + i + "_filtered_contig_annotations.csv" )
        fp = open(path + dataset[4:] + "_" + i + "_all_contig_annotations.csv" )
        fp.readline()
        for line in fp:
            cols = line.rstrip().split(",")
            barcode = cols[0].split("-")[0]
            if (barcode not in barcode10X):
                barcode10X[barcode] = [cols[0], "*", "*", "*"]
            v = cols[6] if cols[6] != "None" else "*"
            d = cols[7] if cols[7] != "None" else "*"
            j = cols[8] if cols[8] != "None" else "*"
            c = cols[9] if cols[9] != "None" else "*"
            chainType = GetChainType(v, j, c)
            
            key = (barcode, chainType) 
            productive = int(cols[11]=="True")
            umis = int(cols[15])
            if (key not in barcodeChainAbund):
                barcodeChainAbund[key] = (productive, umis)
            if ((productive, umis) > barcodeChainAbund[key]):
                barcodeChainAbund[key] = (productive, umis)
            if (barcode not in barcode10Xbt):
                barcode10Xbt[barcode] = {}
            barcode10Xbt[barcode][i] = 1
            
    # Determine the cell type for each 10X cell
    for b in barcode10X.keys():
        maxTag = -1
        for i in range(7):
            if ( (b, i) not in barcodeChainAbund ):
                continue 
            if ( maxTag == -1 or barcodeChainAbund[(b, i)] > barcodeChainAbund[(b, maxTag)] ):
                maxTag = i
        
        if (maxTag <= 2):
            barcodeCellType[b] = 0
        elif (maxTag <= 4):
            barcodeCellType[b] = 1
        else:
            barcodeCellType[b] = 2
        
            
    # Finally put CDR3 to 10X cells 
    for i in ["b", "t"]:
        #fp = open(path + dataset[4:] + "_" + i + "_filtered_contig_annotations.csv" )
        fp = open(path + dataset[4:] + "_" + i + "_all_contig_annotations.csv" )
        fp.readline()
        for line in fp:
            cols = line.rstrip().split(",")
            barcode = cols[0].split("-")[0]
            if (barcode not in barcode10X):
                barcode10X[barcode] = [cols[0], "*", "*", "*"]
            v = cols[6] if cols[6] != "None" else "*"
            d = cols[7] if cols[7] != "None" else "*"
            j = cols[8] if cols[8] != "None" else "*"
            c = cols[9] if cols[9] != "None" else "*"
            chainType = GetChainType(v, j, c)
            cellType = 0
            cellTypeName = ""
            if (chainType <= 2):
                cellTypeName = "B"
                cellType = 0
            elif (chainType <= 4):
                cellTypeName = "abT"
                cellType = 1
            else:
                #print(barcode)
                cellTypeName = "gdT"
                cellType = 2
            
            if (barcodeCellType[barcode] != cellType):
                continue
                
            productive = int(cols[11]=="True")
            umis = int(cols[15])
            #if (barcode == "GACCTGGCATCGATGT"):
            #    print((productive, umis), barcodeChainAbund[(barcode, chainType)], cols[11])
                
            idx = 3 
            if (cols[5] == "IGH" or cols[5] == "TRB" or cols[5] == "TRD"):
                idx = 2
            if (cols[13] == "None"):
                continue
                
            if (cols[12] == "None"):
                cols[12] = "*"
            barcode10X[barcode][1] = cellTypeName
            
            if ((productive, umis) < barcodeChainAbund[(barcode, chainType)]
                   or ((productive, umis)==barcodeChainAbund[(barcode, chainType)] 
                           and barcode10X[barcode][idx]!="*")):
                if (productive == 1):
                    barcode10X[barcode].append( [idx, ",".join([v, d, j, c, cols[13], cols[12], str(umis)])] )
            else:
                barcode10X[barcode][idx] = ",".join([v, d, j, c, cols[13], cols[12], str(umis)])
        fp.close()
    
# Clean up the filtered 10X cells
tmp = barcode10X.copy()
barcode10X = {}
for b in tmp.keys():
    if (tmp[b][1] != "*" or tmp[b][1] == "" ):
        barcode10X[b] =tmp[b][:]
        


# In[ ]:





# In[ ]:





# In[185]:


# cell type decomposition test
if ("vdj" in dataset):
    print( len(barcode10X) )
    display(GetChainState(barcode10X, False))
    
    cellTypeCompositionDf, notInSeurat = GetCellTypeComposition(barcodeSeurat, barcode10X)
    print(notInSeurat)
    display(pd.DataFrame.from_dict(cellTypeCompositionDf))
    
    # Use pie chart to visualize the composition.
    #cellTypeCompositionTrust = GetCellTypeComposition(barcodeSeurat, barcodeTrust)
    #cellTypeComposition10X = GetCellTypeComposition(barcodeSeurat, barcode10X) 
    
    
    # umap the results
    tenXCell = []
    for b in barcodeSeurat.keys():
        if (b in barcode10X):
            clusterName = barcode10X[b][1]  
        else:
            clusterName = "NA"
        tenXCell.append(clusterName)
    Seurat["10X"] = tenXCell

    tenXCell = []
    for b in barcodeSeurat.keys():
        cGene = "*"
        if (b in barcode10X and barcode10X[b][2] != "*"):
            cGene = barcode10X[b][2].split(",")[3]
        tenXCell.append(cGene)
    Seurat["C gene"] = tenXCell
    
    snsFig = sns.relplot(x="UMAP_1", y="UMAP_2", hue="C gene", style = "cluster",
                    data=Seurat.loc[Seurat["cluster"].str.contains("B") & (Seurat["C gene"].str.contains("IGH"))
                                    & ~(Seurat["C gene"].str.contains("IGHM"))],
                    alpha=0.75, palette="bright")
    #snsFig = sns.relplot(x="UMP_1", y="UMP_2", hue="TRUST4",
    #                     data=Seurat.loc[(Seurat["cluster"]=="CD8 T")], alpha=0.75, palette="bright" )




# In[187]:


def IsMatch(A, B, onlyCDR3): 
    if (A == "*" or B == "*"):
        return False
    colsA = A.split(",")
    colsB = B.split(",")
    for i in [0, 2, 3]:
        colsA[i] = GetMainGeneName(colsA[i])
        colsB[i] = GetMainGeneName(colsB[i])
    checkList = [0, 2, 3, 4]
    if (onlyCDR3):
        checkList = [4]
    for i in checkList:
        if (colsA[i] != colsB[i]):
            return False
    return True
        


# In[ ]:





# In[188]:


# Compare TRUST and 10X's result
sns.set(font_scale=1.25)
sns.set_style("white")
if ("vdj" in dataset):
    # Find out how many cells share things.
    sharedCnt = 0
    for b in barcodeTrust.keys():
        if (b in barcode10X):
            sharedCnt += 1
    print(len(barcodeTrust), len(barcode10X), sharedCnt) 
    
    cellTypeAssignment = {"*":{"*":0, "B":0, "abT":0}, "B":{}, "abT":{}, "gdT":{"B":0}}
    for b in barcodeTrust.keys():
        trustCellType = barcodeTrust[b][1]
        if (b not in barcode10X):
            tenCellType = "*"
        else:
            tenCellType = barcode10X[b][1]
        if (tenCellType not in cellTypeAssignment[trustCellType]):
            cellTypeAssignment[trustCellType][tenCellType] = 0
        mixCellType = False
        if ( b in barcode10X and "b" in barcode10Xbt[b] and "t" in barcode10Xbt[b]):
            mixCellType = True
        if (False and trustCellType == "abT" and tenCellType == "B"):
            print("abT->B ", b, " ", barcodeSeurat[b][2] if (b in barcodeSeurat) else "NA", " ", mixCellType)
        if (False and trustCellType == "B" and tenCellType == "abT"):
            print("B->abT ", b, " ", barcodeSeurat[b][2] if (b in barcodeSeurat) else "NA", " ", mixCellType)
        cellTypeAssignment[trustCellType][tenCellType] += 1
    for b in barcode10X:
        if (b not in barcodeTrust):
            cellTypeAssignment["*"][barcode10X[b][1]] += 1
    display(pd.DataFrame.from_dict(cellTypeAssignment))
    
    barcodeMatch = {} 
    # match 0: no match. 1: 10X has no such information. 2: match. -1: 10X and T4 both have no information
    for b in barcodeTrust.keys():
        if (b not in barcodeMatch):
            barcodeMatch[b] = [0, 0]
       
        for i in [2, 3]:
            match = ""
            if (barcodeTrust[b][i] != "*"):
                if (b not in barcode10X):
                    #if (i == 3 and barcodeTrust[b][1] == "B" and barcodeTrust[b][2] == "*"):
                    #            print(barcodeTrust[b])
                    match = "T4 only"
                else:
                    if ( IsMatch(barcodeTrust[b][i], barcode10X[b][i], True)):
                        match = "Match"
                    else:
                        for j in range(4, len(barcode10X[b])):
                            if (barcode10X[b][j][0] == i and IsMatch(barcodeTrust[b][i], barcode10X[b][j][1], True) ):
                                match = "Secondary"
                        if (match == "" and barcode10X[b][i] == "*"):
                            match = "T4 only"
                        if (match == ""):  
                            #if (i == 3 and barcodeTrust[b][1] == "abT" and barcodeTrust[b][2] != "*"):
                            #    print(barcodeTrust[b], barcode10X[b])
                            match = "Mismatch"
            else: # TRUST4 has no result.
                if (b not in barcode10X):
                    match = "Missing"
                else:
                    if (barcode10X[b][i] != "*"):
                        match = "10X only" #"T4 missing"
                    else:
                        match = "Missing"
            barcodeMatch[b][i - 2] = match
    #print(barcode10X["AATCCAGAGCTGCGAA"], "\n", barcodeTrust["AATCCAGAGCTGCGAA"], "\n", barcodeSeurat["AATCCAGAGCTGCGAA"])        
    for chainCount in [1]:  
        #fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        axIdx = 0 
        for cellType in ["B", "abT"]:
            fig, axes = plt.subplots(1, 1, figsize=(7, 5))
            data = {"Chain1":{"Missing":0, "10X only":0, "Mismatch":0, "T4 only":0, "Secondary":0, "Match":0}, 
                    "Chain2":{"Missing":0, "10X only":0, "Mismatch":0, "T4 only":0, "Secondary":0, "Match":0}}
            for b in barcodeTrust.keys():
                if (barcodeTrust[b][1] != cellType
                   or (b in barcode10X and barcode10X[b][1] != cellType) ):
                    continue
                if (b not in barcodeSeurat): # only restricted to the cells retained in Seurat
                    continue
                cnt = 0
                if (barcodeTrust[b][2] != "*"):
                    cnt += 1
                if (barcodeTrust[b][3] != "*"):
                    cnt += 1
                if (cnt < chainCount):
                    continue 
                #if (cellType == "B" and barcodeMatch[b][0] == "Mismatch"):
                #    print(b)
                data["Chain1"][barcodeMatch[b][0]] += 1
                data["Chain2"][barcodeMatch[b][1]] += 1
            print("Chain1 precision:", (data["Chain1"]["Match"]+data["Chain1"]["Secondary"]) / 
                (data["Chain1"]["Match"]+data["Chain1"]["Secondary"]+data["Chain1"]["Mismatch"]))
            print("Chain2 precision:", (data["Chain2"]["Match"]+data["Chain2"]["Secondary"]) / 
               (data["Chain2"]["Match"]+data["Chain2"]["Secondary"]+data["Chain2"]["Mismatch"]))
            match = data["Chain1"]["Match"] + data["Chain2"]["Match"]
            secondary = data["Chain1"]["Secondary"] + data["Chain2"]["Secondary"]
            mismatch = data["Chain1"]["Mismatch"] + data["Chain2"]["Mismatch"]
            print(cellType, "precision:", (match + secondary) / (match + secondary + mismatch))
            data = pd.DataFrame.from_dict(data, orient="index")   
            display(data)    
            data = data.stack().reset_index()
            data.columns = ["Chain", "Category", "Count"]
            data = data[data["Count"] != 0]
            #display(data)
    
            if (cellType == "B"):
                data = data.replace(["Chain1", "Chain2"], ["IGH", "IGL/IGK"])
            elif (cellType == "abT"):
                data = data.replace(["Chain1", "Chain2"], ["TRB", "TRA"])
            elif (cellType == "gdT"):
                data = data.replace(["Chain1", "Chain2"], ["TRD", "TRG"])
            #data = data.replace([-1, 0, 1, 2], ["Both missing", "False", "10X missing", "True"])
            display(data)
            snsFig = sns.barplot(x="Chain", y="Count", hue="Category", data=data, palette="bright",
                       hue_order = ["Mismatch", "T4 only", "Secondary", "Match"], ax = axes)
            if (axIdx > 1):
                snsFig.set(ylabel="")
            snsFig.set(xlabel="")
            axIdx += 1
            #plt.legend()
            #plt.savefig(dataset+"_ChainCount_"+cellType+"_CDR3match.pdf", bbox_inches="tight", format="pdf" )
        plt.show()
        
    # Check the imputed chains.
    for b in barcodeTrust:
        for i in [2, 3]:
            if ( "impute" in barcodeTrust[b][i]):
                print("impute: ", b, barcodeMatch[b][i-2])
    # A cell level comparison bewteen TRUST4 and 10X.
    #for cellType in ["B", "abT"]:
    #    for b in barcodeTrust.keys():
    
    # Chekc how many paired-chain cells in TRUST4 in in Seurat
    cnt = 0
    total = 0
    for b in barcodeTrust:
        if (barcodeTrust[b][2] == "*" or barcodeTrust[b][3] == "*"):
            continue
        total += 1
        if (b in barcodeSeurat):
            cnt += 1
    print(cnt, total)
     
    cnt = 0
    total = 0
    for b in barcode10X:
        if (barcode10X[b][2] == "*" or barcode10X[b][3] == "*"):
            continue
        total += 1
        if (b in barcodeSeurat):
            cnt += 1
    print(cnt, total)
        


# In[189]:


barcodeMatch


# In[190]:


# Analyze only considering the barcode from Seurat
# Compare TRUST and 10X's result
sns.set(font_scale=1.25, style="white")
if ("vdj" in dataset):
    barcodeMatch = {} 
    # match 0: no match. 1: 10X has no such information. 2: match. -1: 10X and T4 both have no information
    for b in barcodeSeurat.keys():
        if (b not in barcodeMatch):
            barcodeMatch[b] = [0, 0]
       
        for i in [2, 3]:
            match = ""
            if (b not in barcodeTrust and b not in barcode10X):
                continue
            elif (b not in barcode10X):
                if (barcodeTrust[b][i] != "*"):
                    match = "T4 only"
                else:
                    continue
            elif (b not in barcodeTrust):
                if (barcode10X[b][i] != "*"):
                    match = "10X only"
                else:
                    continue
            else:
                if (barcodeTrust[b][i] == "*" and barcode10X[b][i] == "*"):
                    continue
                elif (barcode10X[b][i] == "*"):
                    match = "T4 only"
                elif (barcodeTrust[b][i] == "*"):
                    match = "10X only"
                else:
                    if (IsMatch(barcodeTrust[b][i], barcode10X[b][i], True)):
                        match = "Match"
                    else:
                        for j in range(4, len(barcode10X[b])):
                            if (barcode10X[b][j][0] == i and IsMatch(barcodeTrust[b][i], barcode10X[b][j][1], True) ):
                                match = "Secondary"
                        if (match == ""):
                            match = "Mismatch"
            #if (b == "AATCCAGAGCTGCGAA"):
            #    print(match)
            barcodeMatch[b][i - 2] = match
            
 
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    axIdx = 0 
    for cellType in ["B", "abT"]:
        data = {"Chain1":{"10X only":0, "Mismatch":0, "T4 only":0, "Secondary":0, "Match":0}, 
                "Chain2":{"10X only":0, "Mismatch":0, "T4 only":0, "Secondary":0, "Match":0}}
        for b in barcodeSeurat.keys():
            if (b in barcodeTrust and barcodeTrust[b][1] != cellType):
                continue
            if (b in barcode10X and barcode10X[b][1] != cellType):
                continue
            if (barcodeMatch[b][0] != 0):
                #if (cellType == "B" and barcodeMatch[b][0] == "Mismatch"):
                #    print(b)
                data["Chain1"][barcodeMatch[b][0]] += 1
            if (barcodeMatch[b][1] != 0):
                data["Chain2"][barcodeMatch[b][1]] += 1
        data = pd.DataFrame.from_dict(data, orient="index")   
        display(data)    
        data = data.stack().reset_index()
        data.columns = ["Chain", "Category", "Count"]
        data = data[data["Count"] != 0]
        #display(data)
        if (cellType == "B"):
            data = data.replace(["Chain1", "Chain2"], ["IGH", "IGL/IGK"])
        elif (cellType == "abT"):
            data = data.replace(["Chain1", "Chain2"], ["TRB", "TRA"])
        elif (cellType == "gdT"):
            data = data.replace(["Chain1", "Chain2"], ["TRD", "TRG"])
        #data = data.replace([-1, 0, 1, 2], ["Both missing", "False", "10X missing", "True"])
        display(data)
        snsFig = sns.barplot(x="Chain", y="Count", hue="Category", data=data, palette="bright",
                   hue_order = ["10X only", "Mismatch", "T4 only", "Secondary", "Match"], ax = axes[axIdx])
        if (axIdx > 0):
            snsFig.set(ylabel="")
        snsFig.set(xlabel="")
        axIdx += 1
        #plt.legend()
        #plt.savefig(dataset+"_ChainCount_"+str(chainCount)+"_CDR3match.pdf", bbox_inches="tight", format="pdf" )
    plt.show()
            
    # A cell level comparison bewteen TRUST4 and 10X of paired chain cells.
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    axIdx = 0 
    for cellType in ["B", "abT"]:
        data = {"10X only":0, "T4 only":0, "Complementary":0, "Mismatch":0, "Match/Secondary":0}
        for b in barcodeSeurat.keys():
            match = ""
            if (b in barcodeTrust and barcodeTrust[b][1] != cellType):
                continue
            if (b in barcode10X and barcode10X[b][1] != cellType):
                continue
            if (b not in barcodeTrust and b not in barcode10X):
                continue
            if ( (barcodeMatch[b][0] == "Match" or barcodeMatch[b][0] == "Secondary")
                and (barcodeMatch[b][1] == "Match" or barcodeMatch[b][1] == "Secondary") ):
                match = "Match/Secondary"
            elif (barcodeMatch[b][0] == "Mismatch" or barcodeMatch[b][1] == "Mismatch"):
                match = "Mismatch"
            elif (barcodeMatch[b][0] == 0 or barcodeMatch[b][1] == 0):
                continue
            elif (barcodeMatch[b][0] == "T4 only" and barcodeMatch[b][1] == "10X only"):
                match = "Complementary"
            elif (barcodeMatch[b][1] == "T4 only" and barcodeMatch[b][0] == "10X only"):
                match = "Complementary"
            elif (barcodeMatch[b][0] == "T4 only" or barcodeMatch[b][1] == "T4 only"):
                match = "T4 only"
            elif (barcodeMatch[b][0] == "10X only" or barcodeMatch[b][1] == "10X only"):
                match = "10X only"
            else:
                break
            
            data[match] += 1
   
        data = pd.DataFrame.from_dict(data, orient="index")
        data = data.reset_index()
        data.columns = ["Category", "Count"]
        data = data[data["Count"] != 0]
        #data = data.replace([-1, 0, 1, 2], ["Both missing", "False", "10X missing", "True"])
        display(data)
        snsFig = sns.barplot(x="Category", y="Count", data=data, palette="bright", 
                   ax = axes[axIdx])
        if (axIdx > 0):
            snsFig.set(ylabel="")
        snsFig.set(xlabel="")
        axes[axIdx].set_title( "BCR" if (axIdx == 0) else "TCR", fontdict={'fontsize':"medium", "fontweight":"bold"} )
        axIdx += 1
        #plt.legend()
        #plt.savefig(dataset+"_ChainCount_"+str(chainCount)+"_CDR3match.pdf", bbox_inches="tight", format="pdf" )
    plt.show()
    
    


# In[191]:


# Check whether the low sensitivity for TCR from the section above is due to the 
# lower expression in TCRs
if (dataset == "10X_vdj_nextgem_hs_pbmc3"):
    data = [[], [], [], []]
    for b in barcodeTrust:
        if b not in barcodeSeurat:
            continue
        tag = 0
        if (barcodeTrust[b][1] == "abT" or barcodeTrust[b][1] == "gdT"):
            tag = 2
        for chain in [2, 3]:
            if (barcodeTrust[b][chain] == "*"):
                continue
            cols = barcodeTrust[b][chain].split(',')
            data[tag + chain - 2].append(float(cols[6]))
    df = pd.DataFrame(columns=["Chain", "Abundance"])
    for c in range(4):
        name = ["IGH", "IGK/IGL", "TRB", "TRA"]
        #tmpName = [name[c]] * len(data[c])
        #tmpData = data[c]
        tmpDf = pd.DataFrame()
        tmpDf["Chain"] = [name[c]] * len(data[c])
        tmpDf["Abundance"] = data[c]
        df = df.append(tmpDf)
        
    print(sp.stats.ranksums(data[0], data[2]))
    print(sp.stats.ranksums(data[1], data[3]))
    
    fig = plt.subplots(figsize=(7, 7))
    snsFig=sns.boxplot(x="Chain", y="Abundance", data=df)
    snsFig.set(ylim=(0,50))
        


# In[192]:


# Use figure to show the annotation of TRUST4, 10X V(D)J and Seurat
if ("vdj" in dataset):
    # Get the total number of Seurat cells for each cell types.
    cellTypes = {}
    for b in barcodeSeurat.keys():
        clusterName = barcodeSeurat[b][2]
        if (clusterName in cellTypes):
            cellTypes[clusterName] += 1
        else:
            cellTypes[clusterName] = 1
    
    cellTypeComposition = {}
    cellTypeComposition["Sum"] = {'#_of_cells': len(barcodeSeurat)}
    cellTypeComposition["Missing"] = {'#_of_cells': 0}
    
    for b in barcodeTrust.keys():
        cnt = 0 
        if (barcodeTrust[b][2] != "*"):
            cnt += 1
        if (barcodeTrust[b][3] != "*"):
            cnt += 1
        if (barcodeTrust[b][1] == "gdT"):
            continue
        key = ""
        if (cnt < 1):
            continue
        else:
            key = "TRUST4_"+barcodeTrust[b][1]#+"_"+str(cnt)
        if (b not in barcodeSeurat):
            if (key not in cellTypeComposition["Missing"]):
                cellTypeComposition["Missing"][key] = 0
            cellTypeComposition["Missing"][key] += 1
            continue
        clusterName = barcodeSeurat[b][2]
            
        if ( clusterName not in cellTypeComposition ):
            cellTypeComposition[clusterName] = {'Seurat': cellTypes[clusterName] }
        if (key not in cellTypeComposition[clusterName]):
            cellTypeComposition[clusterName][key] = 0
        cellTypeComposition[clusterName][key] += 1
        if (key not in cellTypeComposition["Sum"]):
            cellTypeComposition["Sum"][key] = 0
        cellTypeComposition["Sum"][key] += 1
        
    for b in barcode10X.keys():
        cnt = 0 
        if (barcode10X[b][2] != "*"):
            cnt += 1
        if (barcode10X[b][3] != "*"):
            cnt += 1
        if (barcode10X[b][1] == "gdT"):
            continue
        key = ""
        if (cnt < 1):
            continue
        else:
            key = "10X_V(D)J_"+barcode10X[b][1]#+"_"+str(cnt)
        if (b not in barcodeSeurat):
            if (key not in cellTypeComposition["Missing"]):
                cellTypeComposition["Missing"][key] = 0
            cellTypeComposition["Missing"][key] += 1
            continue
        clusterName = barcodeSeurat[b][2]
            
        if ( clusterName not in cellTypeComposition ):
            cellTypeComposition[clusterName] = {'Seurat': cellTypes[clusterName] }
        if (key not in cellTypeComposition[clusterName]):
            cellTypeComposition[clusterName][key] = 0
        cellTypeComposition[clusterName][key] += 1
        if (key not in cellTypeComposition["Sum"]):
            cellTypeComposition["Sum"][key] = 0
        cellTypeComposition["Sum"][key] += 1
    data = pd.DataFrame.from_dict(cellTypeComposition)
    data = data.reindex(sorted(data.columns, key=CellTypeRank), axis=1)
    data = data.fillna(0)
    display(data)
    
    # visualize the result
    df = pd.DataFrame(columns=["CellType", "Assembler", "value"])
    for c in cellTypeComposition:
        if (c == "Sum" or c == "Missing"):
            continue
        for pred in cellTypeComposition[c]:
            df = df.append({"CellType":c, "Method":pred, "#_of_cells":cellTypeComposition[c][pred]}, ignore_index=True)
    fig = plt.subplots(figsize=(15, 7))
    snsFig = sns.barplot(x="CellType", y="#_of_cells", hue="Method", data=df,  
                         hue_order=["Seurat", "TRUST4_abT", "10X_V(D)J_abT", "TRUST4_B", "10X_V(D)J_B"])
    if ("nsclc" in dataset):
        snsFig.set_xticklabels(snsFig.get_xticklabels(), rotation=90)
    snsFig.legend()
    


# In[193]:


# random checkings
if ("vdj" in dataset and "pbmc" in dataset):
    # whether CD4, CD8 T shared TRA or TRB
    for chain in [2, 3]:
        cd8T = {}
        otherT = {}
        for b in barcodeSeurat:
            if ("CD8" in barcodeSeurat[b][2] and b in barcodeTrust
               and barcodeTrust[b][1] == "abT" and barcodeTrust[b][chain] != "*"):
                c = (barcodeTrust[b][chain].split(","))[4]
                if (c=="*"):
                    continue
                if (c not in cd8T):
                    cd8T[c] = 0 
                cd8T[c] += 1
        found = 0
        total = 0
        for b in barcodeSeurat:
            if (("CD4" in barcodeSeurat[b][2] or "Treg" in barcodeSeurat[b][2])
               and b in barcodeTrust and barcodeTrust[b][1] == "abT" and barcodeTrust[b][chain] != "*" ):
                c = (barcodeTrust[b][chain].split(","))[4]
                if (c in cd8T):
                    found += 1
                if (c not in otherT):
                    otherT[c] = 0
                otherT[c] += 1
                total += 1
        repeat = 0        
        for c in otherT:
            if (otherT[c] > 1):
                repeat += 1
        print(found, total, repeat)


# In[194]:


if ("vdj" in dataset and "pbmc" in dataset):
    # whether CD4, CD8 T shared TRA or TRB
    cd8T = {}
    otherT = {}
    barcode = barcode10X
    tmp = 0
    for b in barcodeSeurat:
        if ("CD8" in barcodeSeurat[b][2] and b in barcode
           and barcode[b][1] == "abT" ):
            if (barcode[b][2] == "*"):
                continue
            c1 = (barcode[b][2].split(","))[4]
            if (c1=="*"):
                continue
            if (barcode[b][3] == "*"):
                continue
            c2 = (barcode[b][3].split(","))[4]
            if (c2 == "*"):
                continue
            c = c1 + "\t" + c2
            if (c not in cd8T):
                cd8T[c] = 0 
            cd8T[c] += 1
            tmp += 1
    found = 0
    total = 0
    alreadyFound = {}
    for b in barcodeSeurat:
        if (("CD4" in barcodeSeurat[b][2] or "Treg" in barcodeSeurat[b][2])
           and b in barcode and barcode[b][1] == "abT"):
            if (barcode[b][2] == "*"):
                continue
            c1 = (barcode[b][2].split(","))[4]
            if (c1=="*"):
                continue
            if (barcode[b][3] == "*"):
                continue
            c2 = (barcode[b][3].split(","))[4]
            if (c2 == "*"):
                continue
            c = c1 + "\t" + c2
            if (c in cd8T):
                found += 1
                if (c not in alreadyFound):
                    alreadyFound[c] = 1
            if (c not in otherT):
                otherT[c] = 0
            otherT[c] += 1
            total += 1
    repeat = 0        
    for c in otherT:
        if (otherT[c] > 1):
            repeat += otherT[c]
    print(found, total, repeat, tmp)
    print(len(alreadyFound))


# In[196]:


sns.set(font_scale=1.5, style="white")
if ("nsclc" in dataset):
    fp = open(path + "10X_vdj_v1_hs_nsclc_5gex_full.fa")
    barcodePairChain = {}
    while True:
        header = fp.readline()
        if (not header):
            break
        seq = fp.readline()
        header = header.rstrip()
        cols = header.split(" ") ;
        cols[2] = float(cols[2])
        barcode = (cols[0][1:].split("-"))[0]
        if (barcode not in barcodePairChain):
            barcodePairChain[barcode] = [[], []]
        idx = 1
        chainType = GetChainType(cols[3], cols[5], cols[6])
        if (chainType == 0 or chainType == 3 or chainType == 5):
            idx = 0 
        if (barcode not in barcodeTrust):
            continue
        cdr3 = cols[9].split("=")[1]
        if (barcodeTrust[barcode][idx + 2] == "*"):
            continue
        testCdr3 = barcodeTrust[barcode][idx + 2].split(",")[4]
        if (cdr3 != testCdr3):
            continue ;
        if (len(barcodePairChain[barcode][idx]) == 0
           or (barcodePairChain[barcode][idx][2] < cols[2] )):
            barcodePairChain[barcode][idx] = cols[:]
            
    # Compare the SHM rate bewteen heavy and light chains in plasma cells
    similarity = [[], []]
    for b in barcodeSeurat:
        if (b not in barcodePairChain):
            continue 
        if (barcodeSeurat[b][2] != "PlasmaCells"):
            continue
        if (len(barcodePairChain[b][0]) == 0 or
           len(barcodePairChain[b][1]) == 0):
            continue
        for idx in [0, 1]:
            v = (barcodePairChain[b][idx][3].split(","))[0]
            s = float(v.split(":")[3])
            similarity[idx].append(s)

    data = pd.DataFrame()
    data["Similarity"] = similarity[0] + similarity[1]
    data["Chain"] = ["IGHV"] * len(similarity[0]) + ["IGKV/IGLV"] * len(similarity[1])
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    sns.set(font_scale=1.5)
    sns.boxplot(x = "Chain", y = "Similarity", data=data, ax = ax)
    
    statistics, pvalue = sp.stats.wilcoxon(x=similarity[0], y=similarity[1])
    ax.text(0.3, 82, "p-value=%.2e"%pvalue, fontsize = 16)
    
    print(len(similarity[0]))
    #plt.savefig( "nsclc_PlasmaChainSimilarity.pdf", bbox_inches="tight", format="pdf" )
    
    # Show the correlation of IGH and IGL
    #fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    data = pd.DataFrame()
    data["IGH"] = similarity[0]
    data["IGK/IGL"] = similarity[1]
    sns.set(font_scale=1.5, style="whitegrid")
    snsFig = sns.relplot(x="IGH", y="IGK/IGL", data=data)
    snsFig.set(xlim=[80, 100.1])
    snsFig.set(ylim=[80, 100.1])
    plt.xlabel("IGHV Similarity")
    plt.ylabel("IGKV/IGLV Similarity")
    r = sp.stats.pearsonr(similarity[0], similarity[1])
    print(r)
    plt.text(90, 83, s="cor=%.3f"%r[0] )
    plt.plot([80, 100], [80, 100], ls="--", c="g")
    #plt.savefig( "nsclc_PlasmaChainSimilarityCorrelation.pdf", bbox_inches="tight", format="pdf" )
    plt.show()
    
    # Compare the SHM rate bewteen heavy chain in plasma and memory B cells
    for chain in [0, 1]:
        similarity = [[], []]
        for b in barcodeSeurat:
            if (b not in barcodePairChain):
                continue
            if (len(barcodePairChain[b][chain]) == 0):
                continue
            idx = -1
            if (barcodeSeurat[b][2] == "MemoryBcells"):
                idx = 0 
            elif (barcodeSeurat[b][2] == "PlasmaCells"):
                idx = 1
            if (idx == -1):
                continue
            v = (barcodePairChain[b][chain][3].split(","))[0]
            s = float(v.split(":")[3])
            similarity[idx].append(s) 
        data = pd.DataFrame()
        data["Similarity"] = similarity[0] + similarity[1]
        data["Type"] = ["MemoryB"] * len(similarity[0]) + ["Plasma"] * len(similarity[1])
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        sns.boxplot(x = "Type", y = "Similarity", data=data, ax = ax)
        #sp.stats.wilcoxon(x=similarity[0], y=similarity[1])


# Compare the TCR/BCR-seq part of Cellranger and TRUST4
if (dataset == "10X_vdj_nextgem_hs_pbmc3"):
    barcodeVDJTrust = {}
    fp = open( path + "/IG_barcode_report.tsv", "r" )
    fp.readline()
    for line in fp:
        line = line.rstrip()
        cols = line.split("\t")
        if (cols[1] != "B"):
            continue
        barcodeVDJTrust[ cols[0].split("-")[0] ] = cols
    fp.close()
    
    fp = open(path + dataset[4:] + "_" + "b" + "_all_contig_annotations.csv" )
    fp.readline()
    barcode10XOne = {} # 10X from one file.
    barcodeChainAbund = {}
    for line in fp:
        cols = line.rstrip().split(",")
        barcode = cols[0].split("-")[0]
        if (barcode not in barcode10XOne):
            barcode10XOne[barcode] = [cols[0], "*", "*", "*"]
        v = cols[6] if cols[6] != "None" else "*"
        d = cols[7] if cols[7] != "None" else "*"
        j = cols[8] if cols[8] != "None" else "*"
        c = cols[9] if cols[9] != "None" else "*"
        chainType = GetChainType(v, j, c)

        key = (barcode, chainType) 
        productive = int(cols[11]=="True")
        umis = int(cols[15])
        if (key not in barcodeChainAbund):
            barcodeChainAbund[key] = (productive, umis)
        if ((productive, umis) > barcodeChainAbund[key]):
            barcodeChainAbund[key] = (productive, umis)
    fp.close()
    
    fp = open(path + dataset[4:] + "_" + "b" + "_all_contig_annotations.csv" )
    fp.readline()
    for line in fp:
        cols = line.rstrip().split(",")
        barcode = cols[0].split("-")[0]
        if (barcode not in barcode10XOne):
            barcode10XOne[barcode] = [cols[0], "*", "*", "*"]
        v = cols[6] if cols[6] != "None" else "*"
        d = cols[7] if cols[7] != "None" else "*"
        j = cols[8] if cols[8] != "None" else "*"
        c = cols[9] if cols[9] != "None" else "*"
        chainType = GetChainType(v, j, c)
        cellType = 0
        cellTypeName = ""
        if (chainType <= 2):
            cellTypeName = "B"
            cellType = 0
        elif (chainType <= 4):
            cellTypeName = "abT"
            cellType = 1
        else:
            #print(barcode)
            cellTypeName = "gdT"
            cellType = 2
            

        productive = int(cols[11]=="True")
        umis = int(cols[15])
        #if (barcode == "GACCTGGCATCGATGT"):
        #    print((productive, umis), barcodeChainAbund[(barcode, chainType)], cols[11])

        idx = 3 
        if (cols[5] == "IGH" or cols[5] == "TRB" or cols[5] == "TRD"):
            idx = 2
        if (cols[13] == "None"):
            continue

        if (cols[12] == "None"):
            cols[12] = "*"
        barcode10XOne[barcode][1] = cellTypeName

        if ((productive, umis) < barcodeChainAbund[(barcode, chainType)]
               or ((productive, umis)==barcodeChainAbund[(barcode, chainType)] 
                       and barcode10XOne[barcode][idx]!="*")):
            if (productive == 1):
                barcode10XOne[barcode].append( [idx, ",".join([v, d, j, c, cols[13], cols[12], str(umis)])] )
        else:
            barcode10XOne[barcode][idx] = ",".join([v, d, j, c, cols[13], cols[12], str(umis)])    
    fp.close() 
        

