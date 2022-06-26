import treeswift as ts
import sys
from subprocess import Popen, PIPE, call
import tempfile

from itertools import groupby

# $1 200K tree
# $2 10K tree
# $3 16S seqs
# $4 occupancy report


def tc_parser(tc_outfp):
    with open(tc_outfp, "r") as tc_output:
        tc_output.readline()
        lines = map(lambda x: x.strip().split('\t'), tc_output.readlines())
    lines_sorted = sorted(lines, key=lambda x: x[1])
    clusters = [(key, [i[0] for i in list(group)]) for key, group in groupby(lines_sorted, lambda x: x[1])]
    return clusters

#t200=ts.read_tree_newick(sys.argv[1])
#t200labels = set([i.label for i in t.traverse_postorder(internal=False)])
t10=ts.read_tree_newick(sys.argv[2])
t10labels = set([i.label for i in t10.traverse_postorder(internal=False)])
with open(sys.argv[3]) as f:
        amps  = set(list(map( lambda x: x.strip(), f.readlines())))
amps_or_t10 = amps.union(t10labels)
amps_and_t10 = amps.intersection(t10labels)

occups = dict()
with open(sys.argv[4]) as f:
    for line in f.readlines():
        oc, tax = line.strip().split(" ")
        occups[tax] = int(oc)
treecluster_out = tempfile.NamedTemporaryFile(delete=False, mode='w+t').name
nldef = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
for tcur in [0.3]:
    s = ["TreeCluster.py", "-i", sys.argv[1], "-m", "max", "-t", str(tcur), "-o", treecluster_out]
    call(s, stdout=nldef, stderr=nldef)
    clusters = tc_parser(treecluster_out)
    select = []
    for n, clus in clusters:
        if n == -1 or n == "-1":
            select += [tag for tag in clus if tag in amps or tag in t10labels]
        else:
            sclus = set(clus).intersection(amps_or_t10)
            if len(sclus) == 0:
                continue
            old = sclus.intersection(t10labels)
            select += list(old)
            old_with_amp = sclus.intersection(amps_and_t10)
            if len(old_with_amp) == 0:
                new = sclus.difference(t10labels)
                if len(new) > 0:
                    select += [sorted([(occups[i], i) for i in new])[-1][1]]

    for i in select:
        print(i)
#checkpoint()

