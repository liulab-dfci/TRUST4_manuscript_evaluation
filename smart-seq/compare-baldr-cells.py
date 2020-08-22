#!/usr/bin/env python

import glob, os
import argparse
import re
from difflib import SequenceMatcher 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

parser = argparse.ArgumentParser(description="Compare BALDR with BASIC/TRUST4", formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage='%(prog)s [options]', add_help=True)
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

file_list = []

# Output contains each cell as a row
# Columns are sorted by V-seq, CDR3-seq, J-seq
# 'b' = BASIC match, 't' = TRUST4 match, '-' = CDR3-seq isn't defined by BALDR, 'x' = No matches

out_file = open('PRJNA412649/results/IGH-baldr-cells.txt', 'a')

trust_dir="PRJNA412649/"+inName+"/trust4/"
basic_dir="PRJNA412649/"+inName+"/basic/"
baldr_dir="PRJNA412649/"+inName+"/baldr/IG-mapped_Unmapped/IgBLAST_quant_sorted_filtered/"


file_count = 0

for file in natural_sort(glob.glob(baldr_dir+"*.IGH.filtered")):

    file_count += 1

    cdr3_basic = False
    v_basic = False
    j_basic = False
    cdr3_trust4 = False
    v_trust4 = False
    j_trust4 = False
    
    with open(file) as f:
        file_name = os.path.basename(file).split('_1.IG-mapped_Unmapped.igblast_tabular.quant.sorted.IGH.filtered')[0]
        file_list.append(file_name)
        lines = f.readlines()
        cdr3_seq_base = lines[0].split('\t')[57]
        if cdr3_seq_base == "-":
            out_file.write("-\t-\t-")
        else:    
            full_seq = lines[0].split('\t')[62].rstrip("\n")
            cdr3_seq = full_seq.split(cdr3_seq_base)[0][-3:] + cdr3_seq_base + full_seq.split(cdr3_seq_base)[1][0:3]
            v_seq = full_seq.split(cdr3_seq)[0]
            j_seq = full_seq.split(cdr3_seq)[1]

            with open(basic_dir+file_name+'.fasta') as f2:
                for basic_name, basic_seq in read_fasta(f2):
                    if cdr3_seq in basic_seq and "heavy_chain" in basic_name:
                        basic_v_seq = basic_seq.split(cdr3_seq)[0]
                        basic_j_seq = basic_seq.split(cdr3_seq)[1]
                        if basic_v_seq.endswith(v_seq) == True:
                            v_basic = True
                        if basic_j_seq.startswith(j_seq) == True: 
                            j_basic = True   
                        cdr3_basic = True

            with open(trust_dir+file_name+'_report.tsv') as f2:
                for line in f2:
                    if '\tIGH' in line:
                        assembly_id = line.split('\t')[8].rstrip('\n')
                        break
                    
                with open(trust_dir+file_name+'_annot.fa') as f3:
                    for trust_name, trust_seq in read_fasta(f3):
                        if trust_name.split(' ')[0] == '>'+assembly_id:
                            trust4_seq = trust_seq
                            break

                    trust4_v_seq = trust4_seq.split(cdr3_seq)[0]
                    trust4_j_seq = trust4_seq.split(cdr3_seq)[1]
                    if trust4_v_seq.endswith(v_seq) == True:
                        v_trust4 = True
                    if trust4_j_seq.startswith(j_seq) == True: 
                        j_trust4 = True   
                    cdr3_trust4 = True

            if cdr3_basic == False and cdr3_trust4 == False:
                out_file.write("x\tx\tx")
            else:
                if v_basic == True:
                   out_file.write("b")
                if v_trust4 == True:
                    out_file.write("t")
                if v_basic == False and v_trust4 == False:
                    out_file.write("x")
                out_file.write("\t")

                if cdr3_basic == True:
                   out_file.write("b")
                if cdr3_trust4 == True:
                    out_file.write("t")
                out_file.write("\t")

                if j_basic == True:
                   out_file.write("b")
                if j_trust4 == True:
                    out_file.write("t")
                if j_basic == False and j_trust4 == False:
                    out_file.write("x")

        out_file.write("\t" + str(file_count) + " " + file_name + "\n")



### LIGHT CHAIN ###


file_list = []

out_file = open('PRJNA412649/results/IGKL-baldr-cells.txt', 'a')

file_count = 0

for file in natural_sort(glob.glob(baldr_dir+"*.IGKL.filtered")):

    file_count += 1

    cdr3_basic = False
    v_basic = False
    j_basic = False
    cdr3_trust4 = False
    v_trust4 = False
    j_trust4 = False
    
    with open(file) as f:
        file_name = os.path.basename(file).split('_1.IG-mapped_Unmapped.igblast_tabular.quant.sorted.IGKL.filtered')[0]
        file_list.append(file_name)
        lines = f.readlines()
        cdr3_seq_base = lines[0].split('\t')[54]
        if cdr3_seq_base == "-":
            out_file.write("-\t-\t-")
        else:    
            full_seq = lines[0].split('\t')[59].rstrip("\n")
            cdr3_seq = full_seq.split(cdr3_seq_base)[0][-3:] + cdr3_seq_base + full_seq.split(cdr3_seq_base)[1][0:3]
            v_seq = full_seq.split(cdr3_seq)[0]
            j_seq = full_seq.split(cdr3_seq)[1]

            with open(basic_dir+file_name+'.fasta') as f2:
                for basic_name, basic_seq in read_fasta(f2):
                    if cdr3_seq in basic_seq and "light_chain" in basic_name:
                        basic_v_seq = basic_seq.split(cdr3_seq)[0]
                        basic_j_seq = basic_seq.split(cdr3_seq)[1]
                        if basic_v_seq.endswith(v_seq) == True:
                            v_basic = True
                        if basic_j_seq.startswith(j_seq) == True: 
                            j_basic = True   
                        cdr3_basic = True

            with open(trust_dir+file_name+'_report.tsv') as f2:
                for line in f2:
                    if '\tIGK' in line or '\tIGL' in line:
                        assembly_id = line.split('\t')[8].rstrip('\n')
                        break
                with open(trust_dir+file_name+'_annot.fa') as f3:
                    for trust_name, trust_seq in read_fasta(f3):
                        if trust_name.split(' ')[0] == '>'+assembly_id:
                            trust4_seq = trust_seq
                            break

                    if cdr3_seq in trust4_seq:

                        trust4_v_seq = trust4_seq.split(cdr3_seq)[0]
                        trust4_j_seq = trust4_seq.split(cdr3_seq)[1]
                        if trust4_v_seq.endswith(v_seq) == True:
                            v_trust4 = True
                        if trust4_j_seq.startswith(j_seq) == True: 
                            j_trust4 = True   
                        cdr3_trust4 = True

            if cdr3_basic == False and cdr3_trust4 == False:
                out_file.write("x\tx\tx")
            else:
                if v_basic == True:
                   out_file.write("b")
                if v_trust4 == True:
                    out_file.write("t")
                if v_basic == False and v_trust4 == False:
                    out_file.write("x")
                out_file.write("\t")

                if cdr3_basic == True:
                   out_file.write("b")
                if cdr3_trust4 == True:
                    out_file.write("t")
                out_file.write("\t")

                if j_basic == True:
                   out_file.write("b")
                if j_trust4 == True:
                    out_file.write("t")
                if j_basic == False and j_trust4 == False:
                    out_file.write("x")

        out_file.write("\t" + str(file_count) + " " + file_name + "\n")
