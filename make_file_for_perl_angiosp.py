import sys
import ete3
import optparse
from ete3 import Tree
import random
import re
from collections import defaultdict
import statistics

tree = ete3.Tree("FinalTree.GENERA.newick", format=0)
leaves = list()
for leaf in tree:
    name = leaf.name
    leaves.append(name)

outf1 = open("ElementHyper.GENERA.ANGIOSPERMS.tsv","w")
with open ("HyperGeneraList.1.Seedplants.tsv") as file1:
    first_line = file1.readline()
    for line in file1:
        line = line.rstrip()
        group = line.split("\t",maxsplit=2)[0]
        sp_list = line.split("\t",maxsplit=2)[1]
        species = sp_list.split(",")
        if group != "any_metal":
            for sp in species:
                if sp in leaves:
                    outf1.write(group+"\t"+sp+"\n")
outf1.close()

        
