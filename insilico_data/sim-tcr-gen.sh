#!/bin/bash

# Read length for RNA-Seq simulation
export READ_LENGTHS="50 75 100 150"

# Number of threads for GNU parallel that may be used when running STAR or 
# other memory consuming job (e.g each STAR thread occupies ~50 Gb of RAM)
export NTHREADS_MEM_RESTICTED=1
# Command for running STAR
export STAR_COMMAND="STAR"

# Output paths for STAR, MiXCR, TRUST, IMREP and VDJER
export STAR_OUTPUT="star"
export MIXCR_OUTPUT="mixcr"
export MIXCR_DEFAULT_OUTPUT="mixcr-default"
export TRUST_OUTPUT="trust"

# path to binaries
export BINARIES="binaries"
# where to build star referenve
export STAR_REFERENCE="${BINARIES}/star_reference"

export REPSEQIO_EXECTUABLE="repseqio"
export MITOOLS_EXECTUABLE="mitools"
export MIXCR_EXECUTABLE="mixcr"
export TRUST_BINARIES="run-trust4"

# tmp directory
export TMP_DIR="tmp"

# mkdirs
for out in $BINARIES $STAR_REFERENCE $STAR_OUTPUT $MIXCR_OUTPUT $TRUST_OUTPUT;
do
	mkdir $out 2>/dev/null || true
done

# Log function for pretty printing
function log(){
	if [ -t 1 ] ; then
		# if tty output -> use colors
		printf "\e[1m\e[31m[%s]" $(date +%T)
		printf "\e[32m %s \e[0m" "$*"
		printf "\n"
	else 
		echo "$@"
	fi
}

# Default seed is fixed, though one can specify seed with --seed
export RANDOM_SEED=3141592
# Illumina profile for art_illumina
export ILLUMINA_PROFILE="HS25"

while [[ $# > 0 ]]
do
	key="$1"
	shift
	case $key in
		#
		--art-seq-sys)
			#
			export ILLUMINA_PROFILE=$1
			shift
			;;
		# random seed for mixcr-test and art_illumina
		--seed)
			#
			export RANDOM_SEED=$1
			shift
			;;
		# number of memory tests
		--m-threads)
			#
			export NTHREADS_MEM_RESTICTED=$1
			shift
			;;
		# don't re-generate synthetic data
		--no-generate)
			#
			NO_GENERATE_SYNTH=true
			;;
		# don't re-run STAR
		--no-star)
			#
			NO_STAR=true
			;;
		# don't run mixct
		--no-mixcr)
			#
			NO_MIXCR=true
			;;
        # don't run trust
        --no-trust)
			#
			NO_TRUST=true
			;;
        # test script for correct work
		--test)
			#
			# Just test the script for correctness, without generating a huge sample of data
			export TEST_RUN=true
			;;
		# clear all generated data and software
		--clear-all)
			#
			ls -1 | grep -v ".sh$" | parallel rm -r -f {}
			exit 1
			;;
		# clear all except STAR reference
		--clear)
			#
			ls -1 | grep -v ".sh$" | grep -v "${BINARIES}" | parallel rm -r -f {}
			exit 1
			;;
		# clear all except STAR reference and BAM files
		--clear-results)
			#
			rm -r -f ${MIXCR_OUTPUT}
			rm -r -f ${TRUST_OUTPUT}
			exit 1
			;;
		#
		*)
			#
			;;
#
esac
#
done


##############################################################################################################################
################################## Downloading reference transcriptome sequences (non V-D-J) #################################
##############################################################################################################################

export HUMAN_TRANSCRIPTS=gencode.v19.pc_transcripts.fa

if [ ! -f ${HUMAN_TRANSCRIPTS} ]; then
	log "Donwloading human $HUMAN_TRANSCRIPTS.gz ..."
	wget --quiet ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/$HUMAN_TRANSCRIPTS.gz

	log "Unzipping $HUMAN_TRANSCRIPTS.gz ..."
	gzip -d ${HUMAN_TRANSCRIPTS}.gz
fi


if [[ ${HUMAN_TRANSCRIPTS} == *"v25"* ]]; then 
	export HUMAN_TRANSCRIPTS_FILTERED=$(echo $HUMAN_TRANSCRIPTS | sed 's/.fa/.filtered.fa/' )
	log "Removing true V-D-J records from $HUMAN_TRANSCRIPTS; filtered file: $HUMAN_TRANSCRIPTS_FILTERED "
	cat ${HUMAN_TRANSCRIPTS} | grep -n -E -A1 'TRAC-204|TRDC-201|TRBC2-20' | sed -n 's/^\([0-9]\{1,\}\).*/\1d/p' | sed -f - $HUMAN_TRANSCRIPTS > $HUMAN_TRANSCRIPTS_FILTERED
else 
	export HUMAN_TRANSCRIPTS_FILTERED=$HUMAN_TRANSCRIPTS
fi 


##############################################################################################################################
###################################### Downloading reference and building index for STAR #####################################
##############################################################################################################################

STAR_HG38_FASTA="GRCh38.p10.genome.fa"
STAR_HG38_GTF="gencode.v26.chr_patch_hapl_scaff.annotation.gtf"



if [ ! -f ${STAR_REFERENCE}/${STAR_HG38_FASTA} ]; then
	log "Downloading hg38 reference for STAR"
	wget --quiet -O ${STAR_REFERENCE}/${STAR_HG38_FASTA}.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/${STAR_HG38_FASTA}.gz
	wget --quiet -O ${STAR_REFERENCE}/${STAR_HG38_GTF}.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/${STAR_HG38_GTF}.gz
	gzip -d ${STAR_REFERENCE}/${STAR_HG38_FASTA}.gz
	gzip -d ${STAR_REFERENCE}/${STAR_HG38_GTF}.gz
fi

STAR_HG38_PATH="${STAR_REFERENCE}/hg38"


declare -A STAR_HG_NAMES
STAR_HG_NAMES[${STAR_HG38_PATH}]="hg38"

mkdir $STAR_HG38_PATH 2>/dev/null || true

log "Generating STAR index for hg38"

${STAR_COMMAND} --runMode genomeGenerate --runThreadN 32 --genomeDir ${STAR_HG38_PATH} --genomeFastaFiles ${STAR_REFERENCE}/${STAR_HG38_FASTA} --sjdbGTFfile ${STAR_REFERENCE}/${STAR_HG38_GTF} --outTmpDir ${TMP_DIR}_1

rm -r -f ${TMP_DIR}*

##############################################################################################################################
################################################### Generating in-silico data ################################################
##############################################################################################################################


export RNASEQ_PREFIX="in_silico_RNA_Seq"
export IN_SILICO_TRB_DATA="in_silico_TRB.fasta"

if [ -z "${NO_GENERATE_SYNTH}" ]; then
	if [ -z "${RANDOM_SEED}" ]; then
		export RANDOM_SEED=$RANDOM
	fi
	log "Random seed: ${RANDOM_SEED}"

	log "Generating in-silico V-D-J records ${IN_SILICO_TRB_DATA}"
	${REPSEQIO_EXECTUABLE} generateClones --seed ${RANDOM_SEED} -a -b -c 1000 murugan | ${REPSEQIO_EXECTUABLE} normalizeClones | ${REPSEQIO_EXECTUABLE} exportCloneSequence -q 10000 -g 'VDJTranscript+CExon1' -d NFeature[CDR3] -d AAFeature[CDR3] > ${IN_SILICO_TRB_DATA}

	export TEST_FCONV=0.1
	# Runs art_illumina with the default options
	function run_art_illumina(){
		input=$1
		output=$2
		options=${@:3}
		if [ ! -z ${TEST_RUN} ]; then
			options="$options --fcov ${TEST_FCONV}"
		fi	

		art_illumina --rndSeed ${RANDOM_SEED} --seqSys ${ILLUMINA_PROFILE} --noALN --paired --mflen 200 --sdev 30  $options --in $input --out $output
	}
	export -f run_art_illumina


	# generate synthetic rna-seq data
	log "Generating RNA-Seq simulated (without any V-D-J records) reads with art_illumina"
	parallel 'run_art_illumina ${HUMAN_TRANSCRIPTS_FILTERED} ${RNASEQ_PREFIX}_{}bp_ --len {} --fcov $((35*{}/50))' ::: ${READ_LENGTHS}
	ls -1 *.fq | parallel mv {} '{= s:_([0-9]).fq$:_R\1.fastq:; =}'


	export TEST_FCONV=1
	log "Generating RNA-Seq simulated reads from in-silico V-D-J records"
	parallel run_art_illumina ${IN_SILICO_TRB_DATA} ${RNASEQ_PREFIX}_TRB_{}bp_ --len {} --fcov 8 ::: ${READ_LENGTHS}
	# randomizing reads order in paired-end fastq files
	parallel -j${NTHREADS_MEM_RESTICTED} ${MITOOLS_EXECTUABLE} -Xmx80g -Xms40g randomize --tmp-dir ${TMP_DIR} --seed ${RANDOM_SEED} ${RNASEQ_PREFIX}_TRB_{}bp_1.fq ${RNASEQ_PREFIX}_TRB_{}bp_2.fq ${RNASEQ_PREFIX}_TRB_{}bp_R1.fastq ${RNASEQ_PREFIX}_TRB_{}bp_R2.fastq ::: ${READ_LENGTHS}
	# remove art-generated data 
	ls -1 *.fq | parallel rm {}

	parallel "cp ${RNASEQ_PREFIX}_{1}bp_R{2}.fastq ${RNASEQ_PREFIX}_no_VDJ_{1}bp_R{2}.fastq" ::: ${READ_LENGTHS} ::: 1 2
	parallel -j1 "cat ${RNASEQ_PREFIX}_TRB_{1}bp_R{2}.fastq | head -n 4000 >> ${RNASEQ_PREFIX}_{1}bp_R{2}.fastq" ::: ${READ_LENGTHS} ::: 1 2

	log "Randomizing reads"
	parallel -j${NTHREADS_MEM_RESTICTED} ${MITOOLS_EXECTUABLE} -Xmx80g -Xms40g randomize --tmp-dir ${TMP_DIR} --seed ${RANDOM_SEED} ${RNASEQ_PREFIX}{2}_{1}bp_R1.fastq ${RNASEQ_PREFIX}{2}_{1}bp_R2.fastq ${RNASEQ_PREFIX}{2}_{1}bp_R1.fastq ${RNASEQ_PREFIX}{2}_{1}bp_R2.fastq ::: ${READ_LENGTHS} ::: "" "_no_VDJ"

    for READ_LENGTH in 50 75 100 150;
        do
            mv ${RNASEQ_PREFIX}_${READ_LENGTH}bp_R1.fastq ${RNASEQ_PREFIX}_${READ_LENGTH}bp_gclone_R1.fastq
            mv ${RNASEQ_PREFIX}_${READ_LENGTH}bp_R2.fastq ${RNASEQ_PREFIX}_${READ_LENGTH}bp_gclone_R2.fastq
            python rename-gclone.py ${RNASEQ_PREFIX}_${READ_LENGTH}bp_gclone_R1.fastq ${RNASEQ_PREFIX}_${READ_LENGTH}bp_R1.fastq
            python rename-gclone.py ${RNASEQ_PREFIX}_${READ_LENGTH}bp_gclone_R2.fastq ${RNASEQ_PREFIX}_${READ_LENGTH}bp_R2.fastq
        done

fi


##############################################################################################################################
############################################### Running STAR to produce BAM files ############################################
##############################################################################################################################

if [ -z "${NO_STAR}" ]; then
	for STAR_HG_PATH in ${STAR_HG38_PATH};
	do
		export STAR_HG_PATH=$STAR_HG_PATH
		export HG_PREFIX=${STAR_HG_NAMES[$STAR_HG_PATH]}

		function run_star(){
			input=$1
			output=$2
			thread=$3
			${STAR_COMMAND} --alignEndsType EndToEnd --outTmpDir ${TMP_DIR}${thread} --runThreadN 16 --genomeDir ${STAR_HG_PATH}/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outStd BAM_SortedByCoordinate --readFilesIn ${input} --outFileNamePrefix ${output} > ${STAR_OUTPUT}/${output}
			rm ${output}Log* 2>/dev/null || true
			rm ${output}SJ* 2>/dev/null || true
			rm -r -f ${TMP_DIR}${thread} 2>/dev/null || true
		}
		export -f run_star

		log "Running STAR with ${HG_PREFIX} for paired-end analysis"
			#parallel --line-buffer -j${NTHREADS_MEM_RESTICTED} 'run_star "${RNASEQ_PREFIX}{2}_{1}bp_R1.fastq ${RNASEQ_PREFIX}{2}_{1}bp_R2.fastq" ${RNASEQ_PREFIX}{2}_{1}bp.${HG_PREFIX}.paired.sorted.bam _{1}_{2}_' ::: ${READ_LENGTHS} ::: "" "_no_VDJ" 
        	parallel --line-buffer -j${NTHREADS_MEM_RESTICTED} 'run_star "${RNASEQ_PREFIX}_no_VDJ_{1}bp_R1.fastq ${RNASEQ_PREFIX}_no_VDJ_{1}bp_R2.fastq" ${RNASEQ_PREFIX}_no_VDJ_{1}bp.${HG_PREFIX}.paired.sorted.bam _{1}__no_VDJ_' ::: ${READ_LENGTHS}

        	parallel --line-buffer -j${NTHREADS_MEM_RESTICTED} 'run_star "${RNASEQ_PREFIX}_{1}bp_R1.fastq ${RNASEQ_PREFIX}_{1}bp_R2.fastq" ${RNASEQ_PREFIX}_{1}bp.${HG_PREFIX}.paired.sorted.bam _{1}__' ::: ${READ_LENGTHS}

		log "Running STAR with ${HG_PREFIX} for single-end analysis"
		parallel --line-buffer -j${NTHREADS_MEM_RESTICTED} 'run_star ${RNASEQ_PREFIX}{2}_{1}bp_R1.fastq ${RNASEQ_PREFIX}{2}_{1}bp.${HG_PREFIX}.single.sorted.bam _{1}_{2}_' ::: ${READ_LENGTHS} ::: "" "_no_VDJ"
	done

	#${STAR_COMMAND} --genomeLoad Remove --genomeDir ${STAR_HG_PATH}/

	log "Indexing BAM files with samtools"
	ls -1 ${STAR_OUTPUT}/* | parallel "samtools index {}"

		log "Running STAR with ${HG_PREFIX} for pair-end analysis"

			parallel --line-buffer -j${NTHREADS_MEM_RESTICTED} 'run_star "${RNASEQ_PREFIX}_{1}bp_R1.fastq ${RNASEQ_PREFIX}_{1}bp_R2.fastq" ${RNASEQ_PREFIX}_{1}bp.${HG_PREFIX}.paired.sorted.bam _{1}__' ::: ${READ_LENGTHS}

		log "Running STAR with ${HG_PREFIX} for single-end analysis"
		parallel --line-buffer -j${NTHREADS_MEM_RESTICTED} 'run_star ${RNASEQ_PREFIX}_{1}bp_R1.fastq ${RNASEQ_PREFIX}_{1}bp.${HG_PREFIX}.single.sorted.bam _{1}__' ::: ${READ_LENGTHS}
	done

	log "Indexing BAM files with samtools"
	ls -1 ${STAR_OUTPUT}/in_silico_RNA_Seq_[0-9]*bp.hg38.*.sorted.bam | parallel "samtools index {}"

fi

##############################################################################################################################
######################################################### Running MiXCR ######################################################
##############################################################################################################################

if [ -z "${NO_MIXCR}" ]; then
	log "Running MiXCR"
	parallel --line-buffer -j1 "${MIXCR_EXECUTABLE} align -OsaveOriginalReads=true -s hsa -f -p rna-seq -OallowPartialAlignments=true --report ${MIXCR_OUTPUT}/${RNASEQ_PREFIX}_{1}bp.paired.alignReport ${RNASEQ_PREFIX}_{1}bp_R1.fastq ${RNASEQ_PREFIX}_{1}bp_R2.fastq ${MIXCR_OUTPUT}/${RNASEQ_PREFIX}_{1}bp.paired.vdjca" ::: ${READ_LENGTHS} 
	parallel --line-buffer -j1 "${MIXCR_EXECUTABLE} align -OsaveOriginalReads=true -s hsa -f -p rna-seq -OallowPartialAlignments=true --report ${MIXCR_OUTPUT}/${RNASEQ_PREFIX}_{1}bp.single.alignReport ${RNASEQ_PREFIX}_{1}bp_R1.fastq ${MIXCR_OUTPUT}/${RNASEQ_PREFIX}_{1}bp.single.vdjca" ::: ${READ_LENGTHS} 
	ls -1 $MIXCR_OUTPUT | grep "vdjca$" | grep -v "rescued" | parallel --line-buffer "${MIXCR_EXECUTABLE} assemblePartial -f --report ${MIXCR_OUTPUT}/{.}.assemblePartialReport ${MIXCR_OUTPUT}/{} ${MIXCR_OUTPUT}/{.}_rescued1.vdjca"
	ls -1 $MIXCR_OUTPUT | grep "vdjca$" | grep -v "rescued" | parallel --line-buffer "${MIXCR_EXECUTABLE} assemblePartial -f --report ${MIXCR_OUTPUT}/{.}.assemblePartialReport ${MIXCR_OUTPUT}/{.}_rescued1.vdjca ${MIXCR_OUTPUT}/{.}_rescued2.vdjca"
	ls -1 $MIXCR_OUTPUT | grep "vdjca$" | grep -v "rescued" | parallel --line-buffer "${MIXCR_EXECUTABLE} extend -f --report ${MIXCR_OUTPUT}/{.}.extendReport ${MIXCR_OUTPUT}/{.}_rescued2.vdjca ${MIXCR_OUTPUT}/{.}_rescued2_extended.vdjca"
	ls -1 $MIXCR_OUTPUT | grep "vdjca$" | grep -v "rescued" | parallel --line-buffer "${MIXCR_EXECUTABLE} assemble -f -ObadQualityThreshold=0 --report ${MIXCR_OUTPUT}/{.}.assembleReport ${MIXCR_OUTPUT}/{.}_rescued2_extended.vdjca ${MIXCR_OUTPUT}/{.}.clns"
	ls -1 $MIXCR_OUTPUT | grep "clns$" | grep -v "rescued" | parallel --line-buffer "${MIXCR_EXECUTABLE} exportClones ${MIXCR_OUTPUT}/{.}.clns ${MIXCR_OUTPUT}/{.}.txt"
	
	ls -1 $MIXCR_OUTPUT | grep "vdjca$" | grep -v "rescued" | parallel --line-buffer "${MIXCR_EXECUTABLE} assemble -f --report ${MIXCR_DEFAULT_OUTPUT}/{.}.assembleReport ${MIXCR_OUTPUT}/{.}_rescued2_extended.vdjca ${MIXCR_DEFAULT_OUTPUT}/{.}.clns"
	ls -1 $MIXCR_DEFAULT_OUTPUT | grep "clns$" | grep -v "rescued" | parallel --line-buffer "${MIXCR_EXECUTABLE} exportClones ${MIXCR_DEFAULT_OUTPUT}/{.}.clns ${MIXCR_DEFAULT_OUTPUT}/{.}.txt"
fi



##############################################################################################################################
######################################################### Running TRUST ######################################################
##############################################################################################################################

if [ -z "${NO_TRUST}" ]; then
	log "Running TRUST"
	ls -1 ${STAR_OUTPUT}/in_silico_RNA_Seq_[0-9]*bp.hg38.*.sorted.bam | xargs -n 1 basename | parallel "perl ${TRUST_BINARIES} -b ${STAR_OUTPUT}/{} -f bcrtcr.fa --ref IMGT+C.fa -o ${TRUST_OUTPUT}/{}"
fi

# Linux readlink -f alternative for Mac OS X
function readlinkUniversal() {
    targetFile=$1

    cd `dirname $targetFile`
    targetFile=`basename $targetFile`

    # iterate down a (possible) chain of symlinks
    while [ -L "$targetFile" ]
    do
        targetFile=`readlink $targetFile`
        cd `dirname $targetFile`
        targetFile=`basename $targetFile`
    done

    # compute the canonicalized name by finding the physical path 
    # for the directory we're in and appending the target file.
    phys_dir=`pwd -P`
    result=$phys_dir/$targetFile
    echo $result
}


os=`uname`
dir=""

case $os in
    Darwin)
        dir=$(dirname "$(readlinkUniversal "$0")")
    ;;
    Linux)
        dir="$(dirname "$(readlink -f "$0")")"
    ;;
    FreeBSD)
        dir=$(dirname "$(readlinkUniversal "$0")")    
    ;;
    *)
       echo "Unknown OS."
       exit 1
    ;;
esac
