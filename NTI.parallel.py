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

group_nsp = len(species1)
group_MNTD = "NA"
chisl = 0
znam = 0
for sp in species1:
    nd = tree&sp
    min_dist = 1000000000
    for sp1 in species1:
        if sp1 != sp:
            dist = nd.get_distance(sp1)
            if dist < min_dist:
                min_dist = dist
    chisl = chisl + min_dist
    znam = znam + 1
group_MNTD = chisl/znam
    

    
leaves = list()
for nd in tree:
    sp = nd.name
    leaves.append(sp)

null_outf = open("NTI."+group+".1.genera.NullDistr.tsv","w")
outf = open("NTI."+group+".1.genera.tsv","w")
mntd = group_MNTD
nsp = group_nsp
null_distr = list()
for i in range(1000):
    leaves1 = leaves
    random.shuffle(leaves1)
    species11 = leaves1[0:nsp]
    chisl = 0
    znam = 0
    for sp in species11:
        nd = tree&sp
        min_dist = 1000000000
        for sp1 in species11:
            if sp1 != sp:
                dist = nd.get_distance(sp1)
                if dist < min_dist:
                    min_dist = dist
        chisl = chisl + min_dist
        znam = znam + 1
    mntd0 = chisl/znam
    null_outf.write(str(mntd0)+"\n")
    null_distr.append(mntd0)
    
mean_null = statistics.mean(null_distr)
sd_null = statistics.stdev(null_distr)
print (group+"\t"+str(nsp)+"\t"+str(len(null_distr))+"\t"+str(mean_null)+"\t"+str(sd_null)+"\n")
NTI = (mntd-mean_null)/sd_null

sorted_null = sorted(null_distr)
pv_left = 1
k = 0
flag = 0
while flag==0:
    null = sorted_null[k]
    if null <= mntd:
        k = k+1
        if k==len(sorted_null):
            flag = 1
    elif k > 0:
        flag = 1
        pv_left = k/len(sorted_null)
    else:
        flag = 1
        pv_left = "less_0.001"

sorted_null.reverse()
pv_right = 1
k = 0
flag = 0
while flag==0:
    null = sorted_null[k]
    if null >= mntd:
        k = k+1
        if k==len(sorted_null):
            flag = 1    
    elif k > 0:
        flag = 1
        pv_right = k/len(sorted_null)
    else:
        flag = 1
        pv_right = "less_0.001"
    

outf.write(group+"\t"+str(mntd)+"\t"+str(mean_null)+"\t"+str(nsp)+"\t"+str(NTI)+"\t"+str(pv_left)+"\t"+str(pv_right)+"\n")
outf.close()





        
