#!/usr/bin/env python

import glob, os
import argparse
import re
from difflib import SequenceMatcher 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

parser = argparse.ArgumentParser(description="Compare TRUST4 and BASIC", formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage='%(prog)s [options]', add_help=True)
required = parser.add_argument_group("required input parameters")
required.add_argument("-i", default=None, help="Input name", required=True)
args = parser.parse_args()

def rev(string):
    return string[::-1]

def compare(a, b):
  count = 0
  for x, y in zip(a, b):
    if x == y:
      count += 1
    else:
      break
  return count
  
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

inName=args.i.rstrip('/')

trust_dir="PRJNA412649/"+inName+"/trust4/"
basic_dir="PRJNA412649/"+inName+"/basic/"
baldr_dir="PRJNA412649/"+inName+"/baldr/IG-mapped_Unmapped/IgBLAST_quant_sorted_filtered/"

match_count, chain1_count, chain2_count, total_count = 0, 0, 0, 0

file_list = []

out_heavy = open('PRJNA412649/results/IGH-TRUST4-cells.txt', 'a')
out_light = open('PRJNA412649/results/IGKL-TRUST4-cells.txt', 'a')

# Output contains cells seperated by column
# Top row is V-seq, bottom row is J-seq
# 0 = match, 1 = mismatch, x's in both rows = CDR3 mismatch

IGH_V, IGH_J, IGKL_V, IGKL_J = [], [], [], []

# Iterate through all barcode reports in TRUST4 directory
for file in natural_sort(glob.glob(trust_dir+"*_report.tsv")):
    filename = os.path.basename(file)
    file_list.append(filename.split('_report.tsv')[0])
    # Parse TRUST4 barcode report
    with open(trust_dir+filename) as f:
        for line in f:
            if '\tIGH' in line:
                assembly_id_heavy = line.split('\t')[8].rstrip('\n')
                trust4_chain1 = line.split('\t')[2]
                break
    with open(trust_dir+filename) as f:            
        for line in f:
            if '\tIGK' in line or '\tIGL' in line:
                assembly_id_light = line.split('\t')[8].rstrip('\n')
                trust4_chain2 = line.split('\t')[2]
                break

    basic_heavy_chains = 0
    basic_light_chains = 0  

    # Open equivalent BASIC output
    with open(basic_dir+filename.split("_report.tsv")[0]+".fasta") as f:

        # Search for heavy chain
        for name, seq in read_fasta(f):        

            if "heavy_chain" in name:
                trust4_seq = ''
                with open(trust_dir+filename.split("_report.tsv")[0]+"_annot.fa") as f2:
                    for name2, seq2 in read_fasta(f2):
                        if name2.split(' ')[0] == '>'+assembly_id_heavy:
                            trust4_seq = seq2

                            if name2.split(' ')[3] == '*':
                                v_length = 0
                            else:
                                v_first = name2.split(' ')[3].split(',')[0]

                                v_length_ref = int(v_first.split('(')[1].split(')')[0])

                                v_end = int(v_first.split('):(')[1].split('-')[1])
                                cdr3_start = int(name2.split(' ')[9].split('(')[1].split('-')[0])
                                v_ref_end = int(v_first.split('-')[-1].split(')')[0])

                                if v_end >= cdr3_start:
                                    v_length = v_length_ref - (v_end - cdr3_start + 1) 
                                else:
                                    v_length = v_length_ref                                
                                v_length = v_length - (v_length_ref - v_ref_end - 1)

                            if name2.split(' ')[5] == '*':
                                j_length = 0
                            else:
                                j_first = name2.split(' ')[5].split(',')[0]

                                j_length_ref = int(j_first.split('(')[1].split(')')[0])

                                j_start = int(j_first.split('):(')[1].split('-')[0])
                                cdr3_end = int(name2.split(' ')[9].split('-')[1].split(')')[0])
                                j_ref_start = int(j_first.split('-')[-2].split('(')[1])

                                if j_start <= cdr3_end:
                                    j_length = j_length_ref - (cdr3_end - j_start + 1)
                                else:
                                    j_length = j_length_ref
                                j_length = j_length - j_ref_start
                            break

                if trust4_seq != '' and trust4_chain1 in seq:

                    left_side_trust4 = trust4_seq.split(trust4_chain1)[0]
                    right_side_trust4 = trust4_seq.split(trust4_chain1)[1]
                    left_side_basic = seq.split(trust4_chain1)[0]
                    right_side_basic = seq.split(trust4_chain1)[1]

                    left_length = compare(rev(left_side_trust4),rev(left_side_basic))

                    right_length = compare(right_side_trust4,right_side_basic)

                    if (float(left_length) >= float(v_length)):
                        IGH_V.append("0")
                    else:
                        IGH_V.append("1")

                    if (float(right_length) >= float(j_length)):
                        IGH_J.append("0")
                    else:
                        IGH_J.append("1")

                    basic_heavy_chains +=1

        # Search for light chain
            elif "light_chain" in name:
                trust4_seq = ''
                with open(trust_dir+filename.split("_report.tsv")[0]+"_annot.fa") as f2:
                    for name2, seq2 in read_fasta(f2):
                        if name2.split(' ')[0] == '>'+assembly_id_light:
                            trust4_seq = seq2

                            if name2.split(' ')[3] == '*':
                                v_length = 0
                            else:
                                v_first = name2.split(' ')[3].split(',')[0]

                                v_length_ref = int(v_first.split('(')[1].split(')')[0])

                                v_end = int(v_first.split('):(')[1].split('-')[1])
                                cdr3_start = int(name2.split(' ')[9].split('(')[1].split('-')[0])
                                v_ref_end = int(v_first.split('-')[-1].split(')')[0])

                                if v_end >= cdr3_start:
                                    v_length = v_length_ref - (v_end - cdr3_start + 1) 
                                else:
                                    v_length = v_length_ref                                
                                v_length = v_length - (v_length_ref - v_ref_end - 1)

                            if name2.split(' ')[5] == '*':
                                j_length = 0
                            else:
                                j_first = name2.split(' ')[5].split(',')[0]

                                j_length_ref = int(j_first.split('(')[1].split(')')[0])

                                j_start = int(j_first.split('):(')[1].split('-')[0])
                                cdr3_end = int(name2.split(' ')[9].split('-')[1].split(')')[0])
                                j_ref_start = int(j_first.split('-')[-2].split('(')[1])

                                if j_start <= cdr3_end:
                                    j_length = j_length_ref - (cdr3_end - j_start + 1)
                                else:
                                    j_length = j_length_ref
                                j_length = j_length - j_ref_start
                            break

                if trust4_seq != '' and trust4_chain2 in seq:

                    left_side_trust4 = trust4_seq.split(trust4_chain2)[0]
                    right_side_trust4 = trust4_seq.split(trust4_chain2)[1]
                    left_side_basic = seq.split(trust4_chain2)[0]
                    right_side_basic = seq.split(trust4_chain2)[1]

                    left_length = compare(rev(left_side_trust4),rev(left_side_basic))

                    right_length = compare(right_side_trust4,right_side_basic)

                    if (float(left_length) >= float(v_length)):
                        IGKL_V.append("0")
                    else:
                        IGKL_V.append("1")

                    if (float(right_length) >= float(j_length)):
                        IGKL_J.append("0")
                    else:
                        IGKL_J.append("1")    

                    basic_light_chains +=1

    if basic_heavy_chains == 0:
        IGH_V.append("x")
        IGH_J.append("x")          
    if basic_light_chains == 0:  
        IGKL_V.append("x")
        IGKL_J.append("x") 


out_heavy.write('\n'+inName+' BASIC\n') 
out_light.write('\n'+inName+' BASIC\n') 

out_heavy.write('\t'.join(map(str, IGH_V))+'\n') 
out_heavy.write('\t'.join(map(str, IGH_J))+'\n') 
 
out_light.write('\t'.join(map(str, IGKL_V))+'\n')
out_light.write('\t'.join(map(str, IGKL_J))+'\n')
