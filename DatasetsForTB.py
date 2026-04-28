import sys
import ete3
import optparse
from ete3 import Tree
import random
import re
from collections import defaultdict
import statistics

parser=optparse.OptionParser()
parser.add_option('-t', '--tree', help='', type='str')
parser.add_option('-i', '--infile', help='', type='str')
parser.add_option('-o', '--outsuffix', help='', type='str')

options, args=parser.parse_args()

group_hypers = defaultdict(list)
with open (options.infile) as table:  #hyper lists for each metal
    header = table.readline()
    for line in table:
        line = line.rstrip()
        group = line.split("\t")[0]
        sps = line.split("\t")[1]
        sps = sps.split(",")
        for sp in sps:
            group_hypers[group].append(sp)


leaves = list()
tree = Tree(options.tree,format=1)
for leaf in tree:
    name = leaf.name
    leaves.append(name)

for group in group_hypers:
    outf = open(group+options.outsuffix,"w")
    hypers = group_hypers[group]
    for leaf in leaves:
        flag = 0
        if leaf in hypers:
            flag = 1
        outf.write(leaf+"\t"+str(flag)+"\n")
    outf.close()

