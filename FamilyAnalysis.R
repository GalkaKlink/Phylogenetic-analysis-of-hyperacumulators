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
	if (!gen%in%nonseed_vect) {
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

st = 1
    

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

##############################Tree for genera
nonseed_vect = c("Chara","Selaginella","Microsorum","Equisetum","Dicranopteris","Salvinia","Azolla","Pityrogramma","Pteris","Adiantum","Pteridium","Thelypteris","Athyrium","Blechnum","Nephrolepis","Microsporum","Leptochilus")

library(readxl)
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

extra_plants = read.csv(file="FamsForOrphanGenera.txt",header=TRUE,sep="\t")

gen_fam = all_plants[,c(2,3)]
gen_fam = rbind(gen_fam,extra_plants)

plants_taxonomy_table1 = merge(x=plants_taxonomy_table,y=gen_fam,by="genus",all.x=FALSE,all.y=FALSE)
orphs = plants_taxonomy_table[!plants_taxonomy_table$genus%in%plants_taxonomy_table1$genus,]

extra_sp = plants_taxonomy_table1[!plants_taxonomy_table1$species%in%all_plants$species,]
nrow(extra_sp)

nrow(plants_taxonomy_table1)
all_plants_plus_hyper = rbind(all_plants,extra_sp)
nrow(all_plants_plus_hyper)

data = read.table(file = paste("DistToClosestHyper.any_metal.1.genera.tsv",sep="")) #this is an output file from one of phylogenetic analyses that is used here to filter non-Angiosperm families (because non-Angiosperms were excluded manually from the tree)
  
colnames(data) = c("metal","plant","closest_hyper","dist")

all_plants_plus_hyper = subset(all_plants_plus_hyper,all_plants_plus_hyper$genus%in%data$plant)



##############ANALYSIS OF FAMILIES
element_list = element_list[order(element_list$element),]
element_list <- rbind(element_list[-3, ], element_list[3, ])

sign_df = data.frame(family = sort(unique(all_plants_plus_hyper$family)))
big_df <- data.frame(family = sort(unique(all_plants_plus_hyper$family)))
results = data.frame()
hyperaccumulation_fams = c()
nonhyperaccumulation_fams = c()
for (i in c(1:nrow(element_list))) {
    metal = element_list[i,1]
    sps = element_list[i,2]
    sp_list = strsplit(sps, ",")[[1]]
    plants_table = all_plants_plus_hyper
    plants_table$is_hyper = ifelse(plants_table$genus%in%sp_list,"yes","no")
    hyper_fraction = sum(plants_table$is_hyper == "yes")/nrow(plants_table)

    fam_hyp_nonhyp = data.frame()
    for (fam in unique(plants_table$family)) {
        sub = subset(plants_table,plants_table$family==fam)
        hyps = sum(sub$is_hyper=="yes")
        nonhyps = sum(sub$is_hyper=="no")
        vect = c(fam,hyps,nonhyps)
        fam_hyp_nonhyp = rbind(fam_hyp_nonhyp,vect)
    }
    colnames(fam_hyp_nonhyp) = c("fam","hyp","nonhyp")
    fam_hyp_nonhyp$hyp = as.numeric(fam_hyp_nonhyp$hyp)
    fam_hyp_nonhyp$nonhyp = as.numeric(fam_hyp_nonhyp$nonhyp)
    rownames(fam_hyp_nonhyp) <- fam_hyp_nonhyp[, 1]
    fam_hyp_nonhyp = fam_hyp_nonhyp[, -1]
    
    #test_result <- chisq.test(fam_hyp_nonhyp)
    test_result <- chisq.test(fam_hyp_nonhyp, simulate.p.value = TRUE, B = 10000)
    pv = test_result$p.value
    tr = as.data.frame(test_result$stdres)

    p_values = 2*pnorm(abs(tr$hyp),lower.tail = FALSE)
    p_adj = p.adjust(p_values, method = "BH")
    tr$p_BH = p_adj          

    signif1 = subset(tr,tr$p_BH < 0.01 & tr$hyp > 0)
    signif2 = subset(tr,tr$BH < 0.01 & tr$hyp < 0)
    hyperaccumulation_fams = c(hyperaccumulation_fams,rownames(signif1))
    nonhyperaccumulation_fams = c(nonhyperaccumulation_fams,rownames(signif2))
    hyper_fams = paste(rownames(signif1),collapse=",")
    nonhyper_fams = paste(rownames(signif2),collapse=",")
    vect = c(metal,pv,hyper_fraction,nrow(signif1),nrow(signif2),hyper_fams,nonhyper_fams)
    results = rbind(results,vect)
   
    tr$family = rownames(tr)
    tr1 = tr[,c(4,1)]
    colnames(tr1) = c("family",metal)
    big_df = merge(big_df,tr1,by="family",all.x=TRUE,all.y=TRUE)
    tr2 = tr[,c(4,3)]
    colnames(tr2) = c("family",metal)
    sign_df = merge(sign_df,tr2,by="family",all.x=TRUE,all.y=TRUE)


}

colnames(results) = c("element","p_value","hyper_fraction","hyper_fams_number","nonhyper_fams_number","hyper_fams","nonhyper_fams")
results = results[order(results$element),]
results$bh_pv = p.adjust(results$p_value,method="BH")
results$p_value = round(as.numeric(results$p_value),digits=4)
results$bh_pv = round(as.numeric(results$bh_pv),digits=4)
results$hyper_fraction = round(as.numeric(results$hyper_fraction),digits=4)

res = results[,c(1,2,3,4,5,8)]
write.table(results,file="HypersFamilies_withFamilyLists.ANGIOSPERMS.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(res,file="HypersFamilies.ANGIOSPERMS.tsv",sep="\t",quote=FALSE,row.names=FALSE)



sign_fams = c(hyperaccumulation_fams ,nonhyperaccumulation_fams)
sign_fams = unique(sign_fams)
filtered_df = subset(big_df,big_df$family%in%sign_fams)

rownames(filtered_df) <- filtered_df[, 1]
filtered_df = filtered_df[, -1]
filtered_df = data.frame(filtered_df)
mean = rowSums(filtered_df[,c(1:16)])/16
filtered_df$mean = mean
filtered_df = filtered_df[order(-filtered_df$mean),]

my_colors <- colorRampPalette(c("blue", "white", "red"))(100)

library(pheatmap)
filtered_df = filtered_df[, -ncol(filtered_df)]
filtered_df = data.matrix(filtered_df)
quantile_breaks <- quantile(filtered_df, probs = c(0.05, 0.95)) #чтобы выбросы не делали карту бесцветной
limit <- round(max(abs(quantile_breaks)),digits=0)
my_breaks <- seq(-limit, limit, length.out = 101)
pheatmap(filtered_df, 
         main = "",
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         #show_rownames = FALSE,
         color = my_colors,
         breaks = my_breaks,
         legend_breaks = c(-limit, 0, limit),
         legend_labels = c(paste("<", -limit), "0", paste(">", limit)),
         fontsize = 16, 
         filename = "BHsignFamilies_heatmap.ANGIOSPERMS.png", width = 8, height = 10,res=300)


filtered_df = subset(sign_df,sign_df$family%in%sign_fams)
filtered_df = data.frame(filtered_df)
filtered_df$metals_over= rowSums(filtered_df[, -ncol(filtered_df)] < 0.01, na.rm = TRUE)
#last_cols <- (ncol(filtered_df) - 1):ncol(filtered_df)
#filtered_df$metals_under <- rowSums(filtered_df[, -last_cols] < -1.96, na.rm = TRUE)

library(ggplot2)
library(gridExtra)

my_theme <- theme_classic() + 
  theme(
    text = element_text(size = 20),      # Общий шрифт
    #plot.title = element_text(hjust = 0.5, face = "bold"), 
    #panel.grid.major = element_line(color = "grey90"),    
   # axis.line = element_line(colour = "black"),           
    #legend.position = "bottom"                              )

ggp <- ggplot(filtered_df, aes(x=metals_over)) + geom_histogram()+ xlab("number of metals with\noverrepresented hyperaccumulators")+ylab("")+my_theme

ggsave(
  filename = "histogram_families.SignBH.ANGIOSPERMS.png", 
  plot = ggp,         
  width = 15,                 
  height = 15,                
  units = "cm",               
  dpi = 300                   
)

filtered_df = filtered_df[order(-filtered_df$metals_over),]
write.table(filtered_df, file="SupplTableFAMS.ANGIOSPERMS.tsv",sep="\t",quote=FALSE,row.names=FALSE)


write.table(big_df, file="FamsResidualsForRegression.tsv",sep="\t",quote=FALSE,row.names=FALSE)

df = filtered_df[,c(1,19)]
df = df[order(-df$metals_over),]


##############
##############Tree for families
##############

sp_ge_fam = data.frame()
for (fam in unique(all_plants_plus_hyper$family)) {
 species = paste(fam,"sp.",sep="_")
 genus = paste(fam,"gen.",sep="_")
 vect = c(species,genus,fam)
 sp_ge_fam = rbind(sp_ge_fam,vect)
}
colnames(sp_ge_fam) = c("species","genus","family")

library("V.PhyloMaker2")
result = phylo.maker(sp_ge_fam, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, output.sp.list = TRUE, output.tree = FALSE, scenarios = "S3", r = 1)
tree=result$scenario.3
write.tree(tree,file="Tree_Of_Families.newick")

sp_list = as.data.frame(result$species.list)
sum(sp_list$status=="prune")
sum(sp_list$status=="bind")

sign_fams_list = paste(sign_fams,collapse=",")  #for NR(T)I calculation

###############
###############Reading results of MPD/MNTD-based analysis for families
###############
setwd("PhylogeneticPaper")
data=read.table(file="NTI.FAMS.tsv")
colnames(data) = c("metal","MPD","MPDnull","genera","NTI","p-left","p-right")
data$MPD = round(data$MPD,digits=1)
data$MPDnull = round(data$MPDnull,digits=1)
data$NTI = round(data$NTI,digits=1)
data = data[order(data$metal),]
data$MPD_p_left_BH = p.adjust(data$"p-left",method="BH")
data$MPD_p_right_BH = p.adjust(data$"p-right",method="BH")
tableNTI = data[,c(1,2,3,5,8,9)]

data=read.table(file="NRI.FAMS.tsv")
colnames(data) = c("metal","MNTD","MNTDnull","genera","NRI","p-left","p-right")
data$MNTD = round(data$MNTD,digits=1)
data$MNTDnull = round(data$MNTDnull,digits=1)
data$NRI = round(data$NRI,digits=1)
data = data[order(data$metal),]
data$MNTD_p_left_BH = p.adjust(data$"p-left",method="BH")
data$MNTD_p_right_BH = p.adjust(data$"p-right",method="BH")
tableNRI = data[,c(1,2,3,5,8,9)]

table_both = merge(x=tableNTI,y=tableNRI,by="metal",all.x=TRUE,all.y=TRUE)
table_both$MPD_p_left_BH = round(table_both$MPD_p_left_BH,digits=3)
table_both$MPD_p_right_BH = round(table_both$MPD_p_right_BH,digits=3)
table_both$MNTD_p_left_BH = round(table_both$MNTD_p_left_BH,digits=3)
table_both$MNTD_p_right_BH = round(table_both$MNTD_p_right_BH,digits=3)

write.table(table_both,file="NTRI.families.BH.tsv")




