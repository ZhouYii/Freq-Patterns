#!/bin/bash
#Lemmatization. Dependency : output format from STANFORD POS tagger is tsv
USAGE="./lemmatize.sh <input_folder>"
if [ "$#" -ne 1 ]; then
    echo $USAGE
    exit 1
fi

CWD=`pwd`/
LEMMATIZATION=$CWD/component/lemmatization/nltk_wordnet_lemmatizer.py
INPUT_DIR=$1/
for file in $INPUT_DIR/*_raw
do
    echo $file
    python $LEMMATIZATION $file
done

