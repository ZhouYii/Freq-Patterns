#!/bin/bash
USAGE="./pos <input_folder> <output_folder>"
if [ "$#" -ne 2 ]; then
    echo $USAGE
    exit 1
fi
CWD=`pwd`/
POS_HOME=$CWD/component/pos_tagger
INPUT_DIR=$1/
OUTPUT_DIR=$2/
OUTPUT_SUFFIX=_OUT
for file in $INPUT_DIR/*_raw
do
    filename=${file:${#INPUT_DIR}+1}
    OUTPUT=$CWD$OUTPUT_DIR$filename$OUTPUT_SUFFIX
    cd $POS_HOME
    ./stanford-postagger.sh ./models/wsj-0-18-bidirectional-nodistsim.tagger $CWD$file > $OUTPUT
    cd $CWD
done
