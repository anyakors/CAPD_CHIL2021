#!/bin/sh

export ARCHS4_DATA=../human_transcript_v8.h5
export COUNTS_DIR=data/counts_by_tissue/
export COUNTS_NORM_DIR=data/counts_by_tissue/norm
export TRANSCRIPT_LEN_DICT=data/transcript_length_dict.jsonl
export GENCODE_LIST=data/GENCODE.v26.hg38.comprehensive
export RBP_LIST=data/RBPDB_transcripts.csv
export RMOD_LIST=data/RNA_modif_LABOME_transcripts.csv
export TRAIN_DATA_DIR=data/main_inputs_train/
export TEST_DATA_DIR=data/main_inputs_test/
export GENES_KEEP_LIST=data/genes_to_keep.csv
export HG38_DIR=data/hg38.fa
export OUT_AUX=data/aux_inputs/

# generate the count tables by tissue source; set min/max sample limit
python tissue_type_parse.py \
     --h5file=$ARCHS4_DATA \
     --savedir=$COUNTS_DIR \
     --min_samples 10 \
     --max_samples 30

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
mv $COUNTS_NORM_DIR/$filext > $COUNTS_DIR

# genome_region: either one chromosome: 'chr1' or 'all'; for all always Killed on my pc ;(
python main_input.py \
     --gencode_list=$GENCODE_LIST \
     --hg38=$HG38_DIR \
     --genome_region chr1 \
     --genes_to_keep=$GENES_KEEP_LIST \
     --input=$COUNTS_DIR \
     --out_train=$TRAIN_DATA_DIR \
     --out_test=$TEST_DATA_DIR \
     --test_ratio 0.2 \
     --context 1000

python aux_input.py \
     --gencode_list=$GENCODE_LIST \
     --rbp_list=$RBP_LIST \
     --rmod_list=$RMOD_LIST \
     --input=$COUNTS_DIR \
     --tr_dict=$TRANSCRIPT_LEN_DICT \
     --output=$OUT_AUX
