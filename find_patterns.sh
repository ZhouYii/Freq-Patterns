#!/bin/bash
##### ECLAT parameters ###
MIN_PATTERN_LENGTH=2
MIN_PATTERN_SUPPORT=7 #Absolute. Can be configured to be percentage
MAX_PATTERN_LENGTH=2
##########################
USAGE="./find_patterns.sh <input_folder> <output_folder>"
if [ "$#" -ne 2 ]; then
    echo $USAGE
    exit 1
fi
CWD=`pwd`/
PATTERN_MINER=$CWD/component/patterns/patterns_no_pos.py
ECLAT_PATH=$CWD/component/eclat/src/eclat
INPUT_DIR=$1/
OUTPUT_DIR=$2/
mkdir .intermediates
./parse.sh $INPUT_DIR .intermediates
python $PATTERN_MINER .intermediates $OUTPUT_DIR/eclat
$ECLAT_PATH -n$MAX_PATTERN_LENGTH -m$MIN_PATTERN_LENGTH -s-$MIN_PATTERN_SUPPORT $OUTPUT_DIR/eclat $OUTPUT_DIR/eclat_patterns
rm -rf .intermediates
