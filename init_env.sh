#!/bin/bash
##### ECLAT parameters ###
MIN_PATTERN_LENGTH=2
MIN_PATTERN_SUPPORT=7 #Absolute. Can be configured to be percentage
MAX_PATTERN_LENGTH=2
##########################

##### Freq_ngrams parameters ###
MIN_SUPPORT=20
N_SIZE=2
#############################

CWD=`pwd`/
INTERFACE_PATH=$CWD/interface_scripts/
INPUT_DIR=$1/
OUTPUT_DIR=$2/

NGRAM_PATH=$CWD/component/ngram/freq_ngram.py
ECLAT_PATH=$CWD/component/eclat/src/eclat
