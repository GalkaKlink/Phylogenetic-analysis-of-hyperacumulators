setwd("CorrelationOfMetals")
data = read.table(file="GenusMetalCorrelation_Hypers_100000H0.Angiosperms")
colnames(data) = c("el1","el2","hypers1","hypers2","common","p_left","p_right","min_common_H0","max_common_H0")
data=subset(data,data$el1 != "element" & data$el2 != "element")
data = subset(data,data$el1 != "any_metal" & data$el2 != "any_metal")
data$p_left_bonf = data$p_left*nrow(data)
data$p_left_bh = p.adjust(data$p_left,method="BH")
sign = subset(data,data$p_left_bh < 0.01)

data$p_right_bh = p.adjust(data$p_right,method="BH")
antisign = subset(data,data$p_right_bh < 0.01)


metals = c(data$el1,data$el2)
metals = unique(metals)


metals = c(sign$el1,sign$el2)
metals = unique(metals)

mt_table = data.frame()
for (mt in metals) {
	sub = subset(sign,sign$el1 == mt | sign$el2 == mt)
	numb = nrow(sub)
	mtls = c(sub$el1,sub$el2)
	mtls = unique(mtls)
	mtls = mtls[mtls != mt]
	mtls_list = paste(mtls,collapse=",")
	vect = c(mt,numb,mtls_list)
	mt_table = rbind(mt_table,vect)
}
colnames(mt_table) = c("metal","correlated_number","correlated_metals")
mt_table$correlated_number = as.numeric(mt_table$correlated_number)
mt_table = mt_table[order(-mt_table$correlated_number),]
print(mt_table,row.names=FALSE)
write.table(mt_table,"GeneraHypers.CorrelatedMetals.ANGIOSPERMS.tsv",quote=FALSE,row.names=FALSE)
