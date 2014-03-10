#!/bin/bash
USAGE="./find_patterns_lemma.sh <input_folder> <output_folder>"
if [ "$#" -ne 2 ]; then
    echo $USAGE
    exit 1
fi
. ./init_env.sh
#PATTERN_MINER=$CWD/component/patterns/patterns_no_pos.py
PATTERN_MINER=$CWD/component/patterns/patterns_global_dict.py
mkdir .intermediates_l
mkdir pos_out_remme #avoid file descriptor thrash
${INTERFACE_PATH}parse.sh $INPUT_DIR .intermediates_l
${INTERFACE_PATH}pos.sh .intermediates_l pos_out_remme
${INTERFACE_PATH}lemmatize.sh pos_out_remme
mv pos_out_remme/* .intermediates_l
python $PATTERN_MINER .intermediates_l $OUTPUT_DIR/eclat_l
#$ECLAT_PATH -n$MAX_PATTERN_LENGTH -m$MIN_PATTERN_LENGTH -s-$MIN_PATTERN_SUPPORT $OUTPUT_DIR/eclat $OUTPUT_DIR/eclat_patterns
python $NGRAM_PATH $OUTPUT_DIR/eclat_l $OUTPUT_DIR/ngram_l_n$N_SIZE $N_SIZE $MIN_SUPPORT
rm -rf pos_out_remme
rm -rf .intermediates_l
