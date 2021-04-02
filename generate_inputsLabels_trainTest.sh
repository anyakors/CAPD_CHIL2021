#!/bin/sh

# has to be downloaded
# https://amp.pharm.mssm.edu/archs4/download.html
# "Expression (transcript level) Human", human_transcript_v8.h5 file
export ARCHS4_DATA=data/human_transcript_v8.h5
export COUNTS_DIR=data/counts_by_tissue/
export COUNTS_NORM_DIR=data/counts_by_tissue/norm
export TRANSCRIPT_LEN_DICT=data/transcript_length_dict.jsonl
export GENCODE_LIST=data/GENCODE.v26.hg38.comprehensive
export RBP_LIST=data/RBPDB_transcripts.csv
export RMOD_LIST=data/RNA_modif_LABOME_transcripts.csv
export TRAIN_DATA_DIR=data/main_inputs_train/
export TEST_DATA_DIR=data/main_inputs_test/
export GENES_KEEP_LIST=data/genes_to_keep.csv
# has to be downloaded
# wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
# gunzip hg38.fa.gz
export HG38_DIR=data/hg38.fa
export OUT_AUX=data/aux_inputs/


SECONDS=0

# =============================================================================================
# GENERATE COUNT TABLES PER TISSUE TYPE
# =============================================================================================
# generate the count tables by tissue source; set min/max sample limit
python tissue_type_parse.py \
     --h5file=$ARCHS4_DATA \
     --savedir=$COUNTS_DIR \
     --min_samples 100 \
     --max_samples 250 \
     --shuffle

# =============================================================================================
# NORMALIZE COUNT TABLES IN R
# =============================================================================================
filext="*.csv"
for f in $COUNTS_DIR$filext
do
  echo "Normalizing $f file"
  # takes one count file per time and puts it in a new dir under the same name;
  # -2.0 is a threshold parameter for outliers; parameters closer to 0 will cut more
  Rscript normalize_counts.r $f $COUNTS_NORM_DIR -2.0
done
echo "Plots for count files are saved in $COUNTS_NORM_DIR/plots folder"

# replacing count tables with new normalized files
rm -f $COUNTS_DIR$filext
mv $COUNTS_NORM_DIR/$filext $COUNTS_DIR

for f in $COUNTS_DIR$filext
do
  sed -i '' 's/""/"samples"/' $f
done

conda deactivate

# =============================================================================================
# GENERATE MAIN INPUT FILES
# =============================================================================================
# PUT TRAIN/TEST CHR MANUALLY
# + THE LONGEST GENES (above length_limit) WILL BE PUT IN TEST_OUT DIR
# train_set=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 \
#             chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 )
# test_set=( chr21 chr22 chrX chrY )

# Spliceai Train/Test split
train_set=( chr2 chr4 chr6 chr8 chr10 chr11 chr12 chr13 chr14 chr15 \
            chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY )
test_set=( chr1 chr3 chr5 chr7 chr9 )

# TRAIN SET
for i in "${train_set[@]}"
do
   python main_input.py \
     --gencode_list=$GENCODE_LIST \
     --hg38=$HG38_DIR \
     --genome_region $i \
     --genes_to_keep=$GENES_KEEP_LIST \
     --input=$COUNTS_DIR \
     --out_train=$TRAIN_DATA_DIR \
     --out_test=$TEST_DATA_DIR \
     --length_limit 50000 \
     --context 1000 \
     --test 0
done

# TEST SET
for i in "${test_set[@]}"
do
   python main_input.py \
     --gencode_list=$GENCODE_LIST \
     --hg38=$HG38_DIR \
     --genome_region $i \
     --genes_to_keep=$GENES_KEEP_LIST \
     --input=$COUNTS_DIR \
     --out_train=$TRAIN_DATA_DIR \
     --out_test=$TEST_DATA_DIR \
     --length_limit 50000 \
     --context 1000 \
     --test 1
done

# =============================================================================================
# GENERATE AUXILIARY INPUT FILES
# =============================================================================================
python aux_input_opt.py \
     --gencode_list=$GENCODE_LIST \
     --rbp_list=$RBP_LIST \
     --rmod_list=$RMOD_LIST \
     --input=$COUNTS_DIR \
     --tr_dict=$TRANSCRIPT_LEN_DICT \
     --output=$OUT_AUX

echo "OVERALL PROCESSING TIME: $SECONDS sec"