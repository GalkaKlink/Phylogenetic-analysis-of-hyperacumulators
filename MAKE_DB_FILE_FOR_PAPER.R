###formatting of table
setwd("C:/WorkBackupFromFlash/VANYA/Hyperaccumulators project/WORK/TO_XGB")
data_hyp = read.table(file="HyperDB.txt",header=TRUE)
data_nonhyp = read.table(file="NonHyperDB.txt",header=TRUE)
data = rbind(data_hyp,data_nonhyp)
colnames(data) = c("PlantSpecies","Symbol","Status")

data$plant_element = paste(data$PlantSpecies,data$Symbol,sep="*")
data = subset(data,data$PlantSpecies != "Not_provided")
#tab = as.data.frame(table(data$plant_element))
data$Status = as.numeric(data$Status)


data$Symbol = ifelse(data$Symbol=="Cr6","Cr",data$Symbol)
data$Symbol = ifelse(data$Symbol=="Sc"|data$Symbol=="Y"|data$Symbol=="Gd"|data$Symbol=="Nd"|data$Symbol=="La"|data$Symbol=="Ce","REEs",data$Symbol)

nonseed_vect = c("Chara","Selaginella","Microsorum","Equisetum","Dicranopteris","Salvinia","Azolla","Pityrogramma","Pteris","Adiantum","Pteridium","Thelypteris","Athyrium","Blechnum","Nephrolepis","Microsporum","Leptochilus")
species = unique(data$PlantSpecies)
plants = species

library("tidyr")
data$ID = data$PlantSpecies
data = separate(data,col="ID",into=c("genus","species"),sep="_")

#############

element_list = data.frame()
cols = character()
cols = c(cols,"species","genus")
for (element in (unique(data$Symbol))) {
    cols = c(cols,paste(element,"hyper",sep="_"))
    cols = c(cols,paste(element,"nonhyper",sep="_"))
}

for (plant in plants) {
    dt = subset(data,data$PlantSpecies == plant)
    gen = dt[1,5]
    vect = c(plant,gen)
    for (element in (unique(data$Symbol))) {
           sub = subset(dt,dt$Symbol==element)
           hyper_papers = sum(sub$Status == 1)
           nonhyper_papers = sum(sub$Status != 1)
           vect = c(vect,hyper_papers,nonhyper_papers)
     }  
     element_list = rbind(element_list,vect)
}
colnames(element_list) = c(cols)


library(readxl)
all_plants = read_excel("../AllPlantsAPG4.xlsx", sheet = 1)
all_plants = as.data.frame(all_plants)
all_plants = subset(all_plants,all_plants$group=="seedplant")
all_plants = all_plants[,c(3,4)]
all_plants = unique(all_plants)

extra_plants = read.csv(file="FamsForOrphanGenera.txt",header=TRUE,sep="\t")
gen_fam = rbind(all_plants,extra_plants)

element_list = merge(x=element_list,y=gen_fam,by="genus",all.x=FALSE,all.y=FALSE)
element_list = element_list[,c(2,1,35,3:34)]
write.table(element_list,file="C:/WorkBackupFromFlash/VANYA/Hyperaccumulators project/WORK/PhylogeneticPaper/Species_Hypers_Nonhypers.DATABASE_SUMMARY.tsv",sep="\t",quote=FALSE,row.names=FALSE)
