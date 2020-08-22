#!/usr/bin/env python

import glob, os
import argparse
import re

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

parser = argparse.ArgumentParser(description="Compare TRUST4 single and bulk", formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage='%(prog)s [options]', add_help=True)
required = parser.add_argument_group("required input parameters")
required.add_argument("-i", default=None, help="Input name", required=True)


args = parser.parse_args()

bulkname=args.i.rstrip('/')
dirname="PRJNA412649/"+bulkname+"/trust4/"

sep_IGKL, bulk_IGH, bulk_IGKL = set(), set(), set()

sep_IGH = []

out_file = open('PRJNA412649/results/bulk-TRUST4.txt', 'a')

# Output has cells seperated by row
# V-seq, CDR3-seq, J-seq are seperated by column
# '.' Indicates a mismatch

file_list = []
skip_baldr_count = 0
skip_baldr_basic_count = 0

# Single Cells
for file in natural_sort(glob.glob("PRJNA412649"+bulkname+"baldr/*_1.IG-mapped_Unmapped.igblast_tabular.quant.sorted.IGH.filtered")): 
    with open(file) as f:
        lines = f.readlines()
        cdr3_seq_base = lines[0].split('\t')[57]
        if cdr3_seq_base == "-":
            skip_baldr_count = skip_baldr_count+1
        else:
            full_seq = lines[0].split('\t')[62].rstrip("\n")
            with open('PRJNA412649'+bulkname+'basic/'+os.path.basename(file).split('_1.IG-mapped_Unmapped.igblast_tabular.quant.sorted.IGH.filtered')[0]+'.fasta') as f2:
                if full_seq in f2.read():

                    cdr3_seq = full_seq.split(cdr3_seq_base)[0][-3:] + cdr3_seq_base + full_seq.split(cdr3_seq_base)[1][0:3]
                    v_seq = full_seq.split(cdr3_seq)[0]
                    j_seq = full_seq.split(cdr3_seq)[1]

                    sep_IGH.append((v_seq, cdr3_seq, j_seq))
                    file_list.append(os.path.basename(file).split('_1.IG-mapped_Unmapped.igblast_tabular.quant.sorted.IGH.filtered')[0])
                else:
                    skip_baldr_basic_count = skip_baldr_basic_count+1

print("BALDR Skipped Cells: " + str(skip_baldr_count))
print("BALDR-BASIC Skipped Cells: " + str(skip_baldr_basic_count))

# Bulk

with open(dirname+bulkname+"_annot.fa") as f:

    for name, seq in read_fasta(f):

        if "IGH" in name:

            cdr3_seq = name.split(' ')[9].split('=')[1]

            if cdr3_seq != "null":

                if name.split(' ')[3] == '*':
                    v_seq = "*"
                else:

                    v_first = name.split(' ')[3].split(',')[0]

                    v_length_ref = int(v_first.split('(')[1].split(')')[0])

                    v_end = int(v_first.split('):(')[1].split('-')[1])
                    cdr3_start = int(name.split(' ')[9].split('(')[1].split('-')[0])
                    v_ref_end = int(v_first.split('-')[-1].split(')')[0])

                    if v_end >= cdr3_start:
                        v_length = v_length_ref - (v_end - cdr3_start + 1) 
                    else:
                        v_length = v_length_ref
                    v_length = v_length - (v_length_ref - v_ref_end - 1)

                    left_side_trust4 = seq.split(cdr3_seq)[0]
                    v_seq = rev(rev(left_side_trust4)[0:v_length])

                if name.split(' ')[5] == '*':
                    j_seq = "*"
                else:

                    j_first = name.split(' ')[5].split(',')[0]

                    j_length_ref = int(j_first.split('(')[1].split(')')[0])

                    j_start = int(j_first.split('):(')[1].split('-')[0])
                    cdr3_end = int(name.split(' ')[9].split('-')[1].split(')')[0])
                    j_ref_start = int(j_first.split('-')[-2].split('(')[1])

                    if j_start <= cdr3_end:
                        j_length = j_length_ref - (cdr3_end - j_start + 1)
                    else:
                        j_length = j_length_ref
                    j_length = j_length - j_ref_start

                    right_side_trust4 = seq.split(cdr3_seq)[1]
                    j_seq = right_side_trust4[0:j_length]


                bulk_IGH.add((v_seq, cdr3_seq, j_seq))


cdr3_count = 0
v_count = 0
j_count = 0

file_count = 0

# Evaluation

for sep_item in sep_IGH:
    cdr3_status = False
    v_status = False
    j_status = False
    for bulk_item in bulk_IGH:       
        if sep_item[1] == bulk_item[1]:
            cdr3_status = True
            bulk_v_status = False
            if bulk_item[0].endswith(sep_item[0]):
                v_status = True
                j_status = False
                bulk_v_status = True
            if bulk_item[2].startswith(sep_item[2]):
                j_status = True
                if bulk_v_status == False and v_status == True:
                    j_status = False
        if cdr3_status == True and v_status == True and j_status == True:
            break               


    if v_status == True:
        v_count += 1
        out_file.write('V\t')
    else:
        out_file.write('.\t')        
    if cdr3_status == True:
        cdr3_count += 1
        out_file.write('CDR3\t')
    else:
        out_file.write('.\t')  
    if j_status == True:
        j_count += 1  
        out_file.write('J\t')
    else:
        out_file.write('.\t')  
    
    out_file.write("\t" + " " + file_list[file_count] + "\n")
    file_count += 1


print ('\n-TRUST4 Sequence recovery-\n')             

print("V Seq recovered:")
print(v_count)
print("CDR3 Seq recovered:")
print(cdr3_count)
print("J Seq recovered")
print(j_count)               


print("\nSingle Cell Length:")
print(len(sep_IGH))
print("Bulk Length")
print(len(bulk_IGH))

