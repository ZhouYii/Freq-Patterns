'''
Things with colons typically noise (Wiki category or user)
Comma most used to separate two locations (for ex: state,country or city,state)
Parenthetical things are noisy : often times may not be location ex : (steamship) 
Replace &amp; with ampersand
'''
import sys
filename = sys.argv[1]
f = open(filename, "r")
word_set = set()
for line in f :
    tokens = line.split(",")
    for token in tokens :
        token = token.strip()
        token = token.replace("&amp;", "&")
        if len(token) == 0 or token.count(":") > 0 :
            continue
        # parenthetical things noisy : sometimes not locations
        if token.count("(") > 0 or token.count(")") > 0 :
            continue
        word_set.add(token)
entry_list = list(word_set)
o = open(filename+"_OUT", "w")
for loc in entry_list :
    o.write(loc+'\n')
o.close()

