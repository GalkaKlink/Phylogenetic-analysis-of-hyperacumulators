###formatting of table
setwd("C:/WorkBackupFromFlash/VANYA/Hyperaccumulators project/WORK/TO_XGB")
data_hyp = read.table(file="HyperDB.txt",header=TRUE)
data_nonhyp = read.table(file="NonHyperDB.txt",header=TRUE)
data = rbind(data_hyp,data_nonhyp)
colnames(data) = c("PlantSpecies","Symbol","Status")

data$plant_element = paste(data$PlantSpecies,data$Symbol,sep="*")
data = subset(data,data$PlantSpecies != "Not_provided")
data$Status = as.numeric(data$Status)

data$Symbol = ifelse(data$Symbol=="Cr6","Cr",data$Symbol)
data$Symbol = ifelse(data$Symbol=="Sc"|data$Symbol=="Y"|data$Symbol=="Gd"|data$Symbol=="Nd"|data$Symbol=="La"|data$Symbol=="Ce","REEs",data$Symbol)

nonseed_vect = c("Chara","Selaginella","Microsorum","Equisetum","Dicranopteris","Salvinia","Azolla","Pityrogramma","Pteris","Adiantum","Pteridium","Thelypteris","Athyrium","Blechnum","Nephrolepis","Microsporum","Leptochilus")
species = unique(data$PlantSpecies)
genera = character()
for (plant in species) {
	dt = subset(data,data$PlantSpecies == plant)
	str = dt[1,]
	sp = str[1,1]
	split_sp <- unlist(strsplit(sp, split = "_"))	
	gen = split_sp[1]
	if (!gen%in%nonseed_vect & gen%in%plants_taxonomy_table1$genus) {
	genera = c(genera,gen)
	}		
}
plants = unique(genera)

library("tidyr")
data$ID = data$PlantSpecies
data = separate(data,col="ID",into=c("genus","species"),sep="_")

#####
#####HYPERS
#####

st = 1 # flag that corresponds to hyperaccumulators    

    hyper_list_all = character()
    hypers = 0
    plants_els = data.frame()

    for (plant in plants) {
	  dt = subset(data,data$genus == plant)
	  dt = dt[order(dt$Status),]
	  str = dt[1,]
	  species = str[1,5]
	  state = as.numeric(str[1,3])
	  if (state == st) {
		hyper_list_all = c(hyper_list_all,species)
		hypers = hypers+1
	  }	
    }

    hyper_list_all = paste(hyper_list_all,collapse=",")
    all_hypers_vector = c("any_metal",hyper_list_all)

#####
#####hyper lists for each element
#####

    element_list = data.frame()
    for (element in (unique(data$Symbol))) {
 	  sub = subset(data,data$Symbol==element & data$Status == st)
 	  gens = unique(sub$genus)
	  gens1 =  gens[gens%in%plants] 
	  hyper_list = paste(gens1,collapse=",")
	  hyper_vect = c(element,hyper_list)
	  element_list = rbind(element_list,hyper_vect)
    }
    colnames(element_list) = c("element","hypers")
    element_list = rbind(element_list,all_hypers_vector)
    write.table(element_list,file=paste("HyperGeneraList.",st,".Seedplants.tsv",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}


##############################Tree for genera
nonseed_vect = c("Chara","Selaginella","Microsorum","Equisetum","Dicranopteris","Salvinia","Azolla","Pityrogramma","Pteris","Adiantum","Pteridium","Thelypteris","Athyrium","Blechnum","Nephrolepis","Microsporum","Leptochilus")

library(readxl)
library("V.PhyloMaker2")
species = unique(data$PlantSpecies)
plants_taxonomy_table = data.frame()
for (plant in species) {
	dt = subset(data,data$PlantSpecies == plant)
	str = dt[1,]
	sp = str[1,1]
	split_sp <- strsplit(sp, split = "_")	
	gen = unlist(split_sp)[1]
	new_sp = paste(gen,"sp.",sep="_")
	vect = c(new_sp,gen)
	plants_taxonomy_table = rbind(plants_taxonomy_table,vect)
}
colnames(plants_taxonomy_table) = c("species","genus")
plants_taxonomy_table = unique(plants_taxonomy_table)
plants_taxonomy_table = subset(plants_taxonomy_table,!plants_taxonomy_table$genus%in%nonseed_vect)
nrow(plants_taxonomy_table)

all_plants = read_excel("../AllPlantsAPG4.xlsx", sheet = 1)
all_plants = as.data.frame(all_plants)
all_plants = subset(all_plants,all_plants$group=="seedplant")
all_plants$species = paste(all_plants$genus,"sp.",sep="_")
all_plants = all_plants[,c(2,3,4)]
all_plants = unique(all_plants)
all_plants = subset(all_plants,all_plants$species != "X_sp.")

extra_plants = read.csv(file="FamsForOrphanGenera.txt",header=TRUE,sep="\t") #This file was manually created for genera from our database that were absent in AllPlantsAPG4.xlsx

gen_fam = all_plants[,c(2,3)]
gen_fam = rbind(gen_fam,extra_plants)

plants_taxonomy_table1 = merge(x=plants_taxonomy_table,y=gen_fam,by="genus",all.x=FALSE,all.y=FALSE)
orphs = plants_taxonomy_table[!plants_taxonomy_table$genus%in%plants_taxonomy_table1$genus,]

extra_sp = plants_taxonomy_table1[!plants_taxonomy_table1$species%in%all_plants$species,]
nrow(extra_sp)

nrow(plants_taxonomy_table1)
all_plants_plus_hyper = rbind(all_plants,extra_sp)
nrow(all_plants_plus_hyper)


result = phylo.maker(all_plants_plus_hyper, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, output.sp.list = TRUE, output.tree = FALSE, scenarios = "S3", r = 1)
tree=result$scenario.3
write.tree(tree,file="GBOTB_with_PlantsFromOurDataBase.SeedPlants.GENERA.newick") 

sp_list = as.data.frame(result$species.list)
sum(sp_list$status=="prune")
sum(sp_list$status=="bind")




