import sys
import ete3
#import optparse
from ete3 import Tree
import random
import re
from collections import defaultdict
import statistics

with open ("TreeOfFamilies_Angiosperms.newick") as infile:
    tree = infile.readline()
    pattern = r"\)[a-zA-Z0-9.]+\:"
    matches = re.findall(pattern,tree)
    print(matches)
    tree = re.sub(pattern,"):",tree)
    pattern = r"\)[a-zA-Z0-9.]+\.rn\.d8s\.tre\:"
    print(matches)
    tree = re.sub(pattern,"):",tree)
    tree = tree.replace("_sp.","")



outf = open("TreeOfFamilies_Angiosperms.cleaned.newick","w")
outf.write(tree)
outf.close()
