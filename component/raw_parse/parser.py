import sys
TEXT_BEGIN = "<TEXT>"
TEXT_END = "</TEXT>"
LABEL_BEGIN = "<ENAMEX"
TYPE_BEGIN = "TYPE="
LABEL_END = "</ENAMEX>"
def skip_preamble(f) :
    '''Place parser market at start of text body'''
    if f.closed :
        print "Parser : file is not open"
        sys.exit(1)
    text = f.read()
    if text.count(TEXT_BEGIN) < 1 or text.count(TEXT_END) < 1 :
        print "Parser : can't markers not compatible for this text "
        sys.exit(1)
    else :
        start = text.find(TEXT_BEGIN)
        end = text.find(TEXT_END)
        return text[start+len(TEXT_BEGIN):end]
    print "Parser : text body not found "

def parse(text) :
    def extract_type(string) :
        start_type_index = string.find(TYPE_BEGIN)
        end_type_index = string.find('">', start_type_index)
        type_name = string[start_type_index+len(TYPE_BEGIN)+1:end_type_index]
        name = string[end_type_index+2:].rstrip('\n')
        name = ' '.join(name.split('\n'))
        return name, type_name

    dictionary = dict()
    while text.count(LABEL_BEGIN) > 0 :
        label_start = text.find(LABEL_BEGIN)
        label_end = text.find(LABEL_END)
        substr = text[label_start:label_end]
        name, label = extract_type(substr)
        text = text.replace(substr, name)
        text = text.replace(LABEL_END, '', 1)
        dictionary[name] = label
    return text, dictionary

if len(sys.argv) != 3 :
    print "Parser arguments : file_path, output_directory"
    sys.exit(1)
filepath = sys.argv[1]
#filepath = "APW19980213.1310.ne.txt"
f = open(filepath, "r")
body_text = skip_preamble(f)
raw_text, dictionary = parse(body_text)
raw_text = ' '.join(raw_text.split('\n'))

#write raw_text to file
parts = filepath.split('/')
filename = parts.pop()
if filename == "" :
    filename = parts.pop()
raw_filepath = sys.argv[2]+"//"+filename+"_raw"
raw_file = open(raw_filepath, "w")
raw_file.write(raw_text)
raw_file.close()

#write dictionary
dict_filepath = sys.argv[2]+"//"+filename+"_dict"
dict_file = open(dict_filepath, "w")
for key in dictionary.keys() :
    if key == '' or dictionary[key] == '' :
        continue
    dict_file.write(key+","+dictionary[key]+"\n")
dict_file.close()

f.close()
