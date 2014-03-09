import sys
from random import randrange
from os import listdir
from os.path import isfile, join
PER_DICT_PATH=
ORG_DICT_PATH=
LOC_DICT_PATH=

def get_filename(filepath) :
    parts = filepath.split('/')
    filename = parts.pop()
    if filename == "" :
        filename = parts.pop()
    return filename

def get_text(filepath) :
    raw_file = open(filepath, "r")
    return raw_file.read()

def restore_dict(dict_path) :
    dict_file = open(dict_path, "r")
    label_set = set()
    dictionary = dict()
    for line in dict_file :
        line = line.strip()
        pair = line.split(',')
        dictionary[pair[0]] = pair[1]
        label_set.add(pair[1])
    return label_set, dictionary

def generate_global_dict() :
    def join_dict(from_dict, into_dict) :
        for key in from_dict.keys() :
            if key in into_dict.keys() :
                continue
            else :
                into_dict[key] = from_dict[key]
        return into_dict

    per_labels, per_dict = restore_dict(PER_DICT_PATH)
    loc_labels, loc_dict = restore_dict(LOC_DICT_PATH)
    org_labels, org_dict = restore_dict(ORG_DICT_PATH)
    full_labels = per_labels.union(loc_lables).union(org_labels)
    full_dict = dict()
    full_dict = join_dict(per_dict, full_dict)
    full_dict = join_dict(loc_dict, full_dict)
    full_dict = join_dict(org_dict, full_dict)
    return full_labels, full_dict

def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

def mine_patterns(sentences, label_set, output_file, seen_phrases, id_counter) :
    def commit_pattern(pattern, output_file) :
        '''Input is a pattern (list). Writes it to file'''
        #first format for eclat
        pattern_str = ""
        for token in pattern :
            pattern_str += str(token)+'\t'
        pattern_str = pattern_str.strip()
        output_file.write(pattern_str+'\n')

    for sentence in sentences :
        accumulator = ""
        pattern = []
        tokens = [t for t in sentence.split(" ") if len(t) > 0]
        for token in tokens :
            if token in label_set :
                if len(accumulator) != 0 :
                    if accumulator in seen_phrases.keys() :
                        pattern.append(seen_phrases[accumulator])
                    else :
                        seen_phrases[accumulator]=id_counter
                        pattern.append(id_counter)
                        id_counter += 1
                    accumulator = ""
                pattern.append(token)
            else :
                accumulator += " "+token
        if len(pattern) > 0 :
            commit_pattern(pattern, output_file)
    return id_counter, seen_phrases


if len(sys.argv) != 3 :
    print "Parser arguments : input_folder, output_file"
    sys.exit(1)
#Get all files in input dir
files = [ f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1],f)) ]
#Filter input dir files by raw text files only
raw_texts = [sys.argv[1]+'/'+f for f in files if f[-4:]=="_raw"]
seen_phrases = dict()
output_file = open(sys.argv[2], "w")
new_id = 0
label_set, dictionary = generate_global_dict()
phrases = sorted(dictionary.keys(), key=lambda x:len(x), reverse=True)
for filepath in raw_texts :
    text = get_text(filepath)
    for phrase in phrases : 
        text = text.replace(phrase, dictionary[phrase])
    sentences = text.split(". ")
    new_id, seen_phrases = mine_patterns(sentences, label_set, output_file, seen_phrases, new_id)
