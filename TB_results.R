library('ape')
setwd("C:/WorkBackupFromFlash/VANYA/Hyperaccumulators project/WORK/PhylogeneticPaper/TreeBreaker")
metals = c("Ag","Al","As","Cd","Co","Cr","Cu","Hg","Mn","Mo","Ni","Pb","REEs","Se","Tl","Zn","any_metal")
#metals = c("Ag","Al","As","Cd","Co","Cr","Cu","Hg","Mn","Mo","Ni","Pb","Se","Tl","Zn","any_metal")
met_tb = data.frame()
for (mt in metals) {
	outfile = paste(mt,".TB_results.genera_1.tsv",sep="")
	tree=read.tree(outfile)
	ntips=length(tree$tip.label)
	l=strsplit(tree$tip.label,'[{}|=]',perl=T)
	edge_index=rep(0,nrow(tree$edge))
	tip_pheno=rep(0,length(l))
	edge_posterior=rep(0,nrow(tree$edge))
	for (i in 1:length(l)) {
 	 	tree$tip.label[i]=l[[i]][1]
  		edge_index[which(tree$edge[,2]==i)]=as.numeric(l[[i]][3])
  		tip_pheno[i]=as.numeric(l[[i]][5])
  		edge_posterior[which(tree$edge[,2]==i)]=as.numeric(l[[i]][7])
	}	
	l=strsplit(tree$node.label,'[{}|=]',perl=T)
	node_posterior=rep(0,length(l))
	for (i in 1:length(l)) {
  		edge_index[which(tree$edge[,2]==i+ntips)]=as.numeric(l[[i]][3])
 		edge_posterior[which(tree$edge[,2]==(i+ntips))]=as.numeric(l[[i]][5])
	}


	post1 = sum(edge_posterior > 0.1)
	post2 = sum(edge_posterior > 0.25)
	post3 = sum(edge_posterior > 0.5)
	post4 = sum(edge_posterior > 0.9)

	vect = c(mt,post1,post2,post3,post4)
	met_tb = rbind(met_tb,vect)
}
	
colnames(met_tb) = c("metal","post10","post25","post50","post90")
print(met_tb,row.names=FALSE)

write.table(met_tb,file="TB_TABLE.genera_1.tsv",sep="\t",quote=FALSE,row.names=FALSE)




#Read the tree and information attached to it
outfile='any_metal.TB_results.genera_1.tsv'
tree=read.tree(outfile)
ntips=length(tree$tip.label)
l=strsplit(tree$tip.label,'[{}|=]',perl=T)
edge_index=rep(0,nrow(tree$edge))
tip_pheno=rep(0,length(l))
edge_posterior=rep(0,nrow(tree$edge))
for (i in 1:length(l)) {
  tree$tip.label[i]=l[[i]][1]
  edge_index[which(tree$edge[,2]==i)]=as.numeric(l[[i]][3])
  tip_pheno[i]=as.numeric(l[[i]][5])
  edge_posterior[which(tree$edge[,2]==i)]=as.numeric(l[[i]][7])
}
l=strsplit(tree$node.label,'[{}|=]',perl=T)
node_posterior=rep(0,length(l))
for (i in 1:length(l)) {
  edge_index[which(tree$edge[,2]==i+ntips)]=as.numeric(l[[i]][3])
  edge_posterior[which(tree$edge[,2]==(i+ntips))]=as.numeric(l[[i]][5])
}

#hist(edge_posterior)
sum(edge_posterior > 0.1)
sum(edge_posterior > 0.25)
sum(edge_posterior > 0.5)

tab = cbind(edge_index,edge_posterior)
tab = as.data.frame(tab)
sub = subset(tab,tab$edge_posterior > 0.5)
print(sub,row.names=FALSE)



#Read rest of file
t=read.table(outfile,comment.char='(')
states =as.matrix(t[,2:(ncol(t)-1)])
lambdas=as.vector(t[,ncol(t)])
#Calculate Bayes Factor
bf=length(which(lambdas>0))/length(which(lambdas==0))
cat('Bayes factor of model with one or more change points to model with no change point is: ',bf)

lambdas1 = lambdas[lambdas > 0]
hist(lambdas1)
mean(lambdas1)




#Plot tree showing phenotype and posterior
par(mfrow=c(1,1))
ec=edge_posterior
w=which(tree$edge[,1]==(ntips+1));if (length(w)==2) ec[w]=max(ec[w])
plot.phylo(tree,show.tip.label = F,edge.color=rgb(ec,0,0),edge.width=1+ec*10)
ncols=length(unique(tip_pheno))
tiplabels(NULL,pch=16,col=rainbow(2*ncols)[ncols+tip_pheno])
stop("Stopping here")

#Plot some MCMC traces
#par(mfrow=c(1,2))
#plot(lambdas,type = 'l',xlab='Sampled iteration',ylab='lambda')
#plot(rowSums(states)-1,type='l',xlab='Sampled iteration',ylab='Number of changepoints')

#Plot ten MCMC states
par(mfrow=c(2,5))
for (s in seq(nrow(states)/10,nrow(states),nrow(states)/10)) {
  ec=rep(0,nrow(tree$edge))
  for (i in 1:nrow(tree$edge)) {
    ec[i]=states[s,edge_index[i]+1]
  }
  w=which(tree$edge[,1]==(ntips+1));if (length(w)==2) ec[w]=max(ec[w])
  plot.phylo(tree,show.tip.label = F,edge.color=rgb(ec,0,0),edge.width=1+ec*10)
  tiplabels(NULL,pch=16,col=rainbow(2*ncols)[ncols+tip_pheno])
}


#plot the correlation between the change points
par(mfrow=c(1,1))
image(t(states[,-ncol(states)]), axes=FALSE)
axis(1, at=seq(0,1,length.out=10), labels= floor(seq(0,ncol(states)-1,length.out=10) ))
axis(2, at=seq(0,1,length.out=10), labels= floor(rev( seq(0,nrow(states),length.out=10)) ))