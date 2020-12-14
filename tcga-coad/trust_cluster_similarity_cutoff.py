# Li Song: the program to re-cluster the CDR3 regions

#!/usr/bin/env python3

import sys
import random

def GetChainType(v, j, c):
	s = ""
	if (v != "*"):
		s = v
	elif (c != "*"):
		s = c 
	elif (j != "*"):
		s = j
	else:
		return -1
	
	if (s[0:3] == "IGH"):
		return 0
	elif (s[0:3] == "IGK"):
		return 1
	elif (s[0:3] == "IGL"):
		return 2
	elif (s[0:3] == "TRA"):
		return 3
	elif (s[0:3] == "TRB"):
		return 4
	elif (s[0:3] == "TRG"):
		return 5
	elif (s[0:3] == "TRD"):
		return 6
	else:
		return -1

def GetMainGeneName(g):
	return g.split("*")[0] 

def GetSimilarity(a, b):
	if (len(a) != len(b)):
		return 0
	diffCnt = 0
	for i in range(len(a)):
		if (a[i] != b[i]):
			diffCnt += 1
	return 1 - diffCnt / len(a)


def PrintSimilarities(cdr3Dict, pairCnt):
	vjCDR3LenList = {}
	# Reorganize the CDR3 list by their V, J and CDR3 length.
	for cdr3 in cdr3Dict:
		key = (GetMainGeneName(cdr3[2]), GetMainGeneName(cdr3[3]), len(cdr3[0]))
		if (key not in vjCDR3LenList):
			vjCDR3LenList[key] = []
		vjCDR3LenList[key].append( cdr3 )
	
	# Collect the total number of possible pairs
	totalCnt = 0
	for l in vjCDR3LenList.values():
		totalCnt += int(len(l) * (len(l) - 1) / 2)

	# Get the sampled down coordinate
	random.seed(17)
	if (pairCnt > totalCnt):
		pairCnt = totalCnt
	samplePoints = random.sample(range(0, totalCnt), pairCnt)
	samplePoints.sort()
	#print(samplePoints, totalCnt, len(vjCDR3LenList))
	
	# Print the similaries
	psum = 0
	tag = 0
	for l in vjCDR3LenList.values():
		cnt = int(len(l) * (len(l) - 1) / 2)
		if (tag >= len(samplePoints)):
			break
		if (samplePoints[tag] >= psum + cnt):
			psum += cnt
			continue
		i = 0
		localPsum = psum
		while (i < len(l) - 1):
			while (tag < len(samplePoints) and localPsum + (len(l) - i - 1) > samplePoints[tag]):
				j = samplePoints[tag] - localPsum
				#print(i, j)
				if (GetSimilarity(l[i][0], l[i + j + 1][0]) == 1):
					print(l)
				print(GetSimilarity(l[i][0], l[i + j + 1][0]))
				tag += 1
			localPsum += (len(l) - i - 1)
			i += 1
		psum += cnt
		#print(psum, localPsum)			

if (__name__ == "__main__"):
	if (len(sys.argv) <= 1):
		print("usage: a.py input_files num_of_pairs [OPTIONS]") 
		exit(1)
	#fpList = open(sys.argv[1])
	#for fname in fpList:
	#	fname = fname.rstrip()
	fp = open(sys.argv[1])
	cdr3s = {}
	for line in fp:
		if (line[0] == "#" or line[0] == "c"):
			continue
		line = line.rstrip("\n")
		cols = line.split("\t")
		skip = False
		if (GetChainType(cols[4], cols[6], cols[7]) != 0):
			continue
			
		for g in [4, 6]: # Must have V, J genes.
			if ( cols[g] == "*"):
				skip = True 
		if (skip or cols[3] == "partial"):
			continue
		key = (cols[2], cols[3], GetMainGeneName(cols[4]), GetMainGeneName(cols[6]))
		cdr3s[key] = 1
	PrintSimilarities(cdr3s, int(sys.argv[2]))
