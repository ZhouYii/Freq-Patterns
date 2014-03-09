#!/bin/bash
USAGE="./find_patterns_lemma.sh <input_folder> <output_folder>"
if [ "$#" -ne 2 ]; then
    echo $USAGE
    exit 1
fi
PATTERN_MINER=$CWD/component/patterns/patterns_no_pos.py
ECLAT_PATH=$CWD/component/eclat/src/eclat
mkdir .intermediates
mkdir pos_out_remme #avoid file descriptor thrash
./parse.sh $INPUT_DIR .intermediates
./pos.sh .intermediates pos_out_remme
./lemmatize.sh pos_out_remme
mv pos_out_remme/* .intermediates
python $PATTERN_MINER .intermediates $OUTPUT_DIR/eclat 
$ECLAT_PATH -n$MAX_PATTERN_LENGTH -m$MIN_PATTERN_LENGTH -s-$MIN_PATTERN_SUPPORT $OUTPUT_DIR/eclat $OUTPUT_DIR/eclat_patterns
rm -rf pos_out_remme
rm -rf .intermediates
