import sys
import ete3
import optparse
from ete3 import Tree
import random
import re
from collections import defaultdict
import statistics

parser=optparse.OptionParser()
parser.add_option('-l', '--line_number', help='', type='int') 
options, args=parser.parse_args()
line_number = options.line_number

tree = ete3.Tree("FinalTree.GENERA.newick", format=0)

k = 0
group = "NA"
sp_list = list()
with open ("HyperGeneraList.1.Seedplants.tsv") as file1:
    for line in file1:
        k += 1
        if k == line_number:
            line = line.rstrip()
            group = line.split("\t",maxsplit=2)[0]
            sp_list = line.split("\t",maxsplit=2)[1]
            


species = sp_list.split(",")
species1 = list()
for sp in species:
    if (len(tree.search_nodes(name=sp))==1):
        species1.append(sp)

outf = open("DistToClosestHyper."+group+".1.genera.tsv","w")
for sp in tree:
    sp_name = sp.name
    min_dist = 1000000000
    closest_hyper = "NA"
    for sp1 in species1:
        if sp1 != sp_name:
            nd = tree&sp1
            dist = sp.get_distance(nd)
            if dist < min_dist:
                min_dist = dist
                closest_hyper = sp1
    outf.write(group+"\t"+sp_name+"\t"+closest_hyper+"\t"+str(min_dist)+"\n")
    
outf.close()





        
