'''
what about hyphenated words
'''
import sys
from nltk.stem.wordnet import WordNetLemmatizer
lemma = WordNetLemmatizer()

def to_wordnet_pos(typestr) :
    if typestr.count("JJ") >= 1 :
        return 'a'
    if typestr.count("RB") >= 1 :
        return 'r'
    if typestr[0] == 'N' :
        return 'n'
    if typestr[0] == 'V' :
        return 'v'
    '''If the type is not supported by NLTK morphy, then use the default 'n' POS
    option'''
    return 'n'


if len(sys.argv) != 2 :
    print "arguments : path_to_text_file"
    sys.exit(1)

f = open(sys.argv[1], 'r')
word_tokens = []
for line in f :
    line = line.strip()
    if len(line) == 0 :
        continue
    line = line.rsplit('\t', 1)
    lemmatized = lemma.lemmatize(line[0].strip(), to_wordnet_pos(line[1].strip()))
    word_tokens.append(lemmatized)

lemmatized = " ".join(word_tokens)
f.close
f = open(sys.argv[1], 'w') 
f.write(lemmatized) 
f.close()

