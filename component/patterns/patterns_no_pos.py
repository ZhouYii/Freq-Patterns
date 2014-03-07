import sys
from random import randrange
from os import listdir
from os.path import isfile, join
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

files = [ f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1],f)) ]
raw_texts = [sys.argv[1]+'/'+f for f in files if f[-4:]=="_raw"]
seen_phrases = dict()
output_file = open(sys.argv[2], "w")
new_id = 0
for filepath in raw_texts :
    filename = get_filename(rreplace(filepath, "_raw", '', 1))
    dict_path = rreplace(filepath, "_raw", "_dict", 1)

    #recover dictionary
    label_set, dictionary = restore_dict(dict_path)
    #recover raw_text
    text = get_text(filepath)

    #substitute labelled phrases into raw_text
    phrases = sorted(dictionary.keys(), key=lambda x: len(x), reverse=True)
    for phrase in phrases : 
        text = text.replace(phrase, dictionary[phrase])

    sentences = text.split(". ")
    new_id, seen_phrases = mine_patterns(sentences, label_set, output_file, seen_phrases, new_id)
