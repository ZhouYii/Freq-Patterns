#!/bin/bash
USAGE="./parse.sh <input_folder> <output_folder>"
if [ "$#" -ne 2 ]; then
    echo $USAGE
    exit 1
fi
CWD=`pwd`/
PARSER=$CWD/component/raw_parse/parser.py
INPUT_DIR=$1/
OUTPUT_DIR=$2/
for file in $INPUT_DIR/*
do
    OUTPUT=$CWD$OUTPUT_DIR$filename
    python $PARSER $file $OUTPUT
done
