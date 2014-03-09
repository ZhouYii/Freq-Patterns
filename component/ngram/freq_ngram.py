#Interface : <input patterns (1 per line)> <n value> <how many top-occurring to return >
from collections import Counter
import sys
if len(sys.argv) != 5 :
    print "arguments : input_pattern_file, out_file,  window_size_n, min_support"
    sys.exit(1)

f = open(sys.argv[1], 'r')
f_out = open(sys.argv[2], 'w')
n_size = int(sys.argv[3])
min_sup = int(sys.argv[4])
global_counter = Counter()
for line in f :
    tok = [t for t in line.split('\t') if len(t) != 0]
    local_counter = Counter(" ".join(tok[i:i+n_size]) for i in range(len(tok)-n_size))
    global_counter += local_counter
all_grams = global_counter.most_common()
for entry in all_grams :
    if entry[1] >= min_sup :
        f_out.write(str(entry[0])+ " SUP:"+str(entry[1])+'\n')
f.close()
