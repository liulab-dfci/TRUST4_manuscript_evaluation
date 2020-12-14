# Evaluations of TRUST4 
This repository contains the scripts used in the evaluation of TRUST4 manuscript.

## Commands of running immune repertoire reconstruction methods
Mixcr v3.0.12:

	$mixcr align -t 8 -p rna-seq -s hsa -OallowPartialAlignments=true $R1 $R2 alignments.vdjca
	$mixcr assemblePartial alignments.vdjca alignments_rescued_1.vdjca
	$mixcr assemblePartial alignments_rescued_1.vdjca alignments_rescued_2.vdjca
	$mixcr extend alignments_rescued_2.vdjca alignments_rescued_2_extended.vdjca
	$mixcr assemble --write-alignments alignments_rescued_2_extended.vdjca clones.clna
	$mixcr assembleContigs --report report.txt clones.clna full_clones.clns
	$mixcr exportClones full_clones.clns full_clones.txt

CATT:
	
	rm ~/CATT/config.jl
	~/CATT/catt -t 8 --f1 $R1 --f2 $R2 -o $output --bowt 8 --chain {IGH, TRB, TRA}

TRUST3:
	trust -f $dataPath/rna_seq_bam_star_${l}bp/FZ-${s}.Aligned.sortedByCoord.out.bam -n 8 -g hg38 {-B} {-B -L}	

TRUST4:

For BAM input:
	
	./run_trust4 -t 8 -b $BAM -f bcrtcr.fa --ref IMGT+C.fa

For FASTQ input:

	./run_trust4 -t 8 -1 $R1 -2 $R2 -f bcrtcr.fa --ref IMGT+C.fa

For 10X Genomics data input:
	
	./run_trust4 -t 8 -b $BAM -f bcrtcr.fa --ref IMGT+C.fa --barcode CB

For SMART-seq data (each input file is for a cell)
	
	./run_trust4 -t 8 -1 $R1 -2 $R2 -f bcrtcr.fa --ref IMGT+C.fa --skipMateExtension

## Scripts 
