require(ggplot2)
require(reshape2)
require(scales)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

setwd("~/Workspace/btol/")

acc_sca <- read.table("data/results_apples.csv", header=F)
colnames(acc_sca) <- c("error", "query")
q <- ggplot(aes(x=error),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
  theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
ggsave("./figures/cdf-nuc.pdf",width=3.5, height=3.5)

 ###########################

acc_sca <- read.table("data/allres.csv", header=F)
colnames(acc_sca) <- c("error", "query","trim")
head(acc_sca)
q <- ggplot(aes(x=error,color=trim),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
  theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
  scale_colour_manual(name = "", values=cbPalette[c(2,3,4,9)])
ggsave("./figures/cdf-nuc.pdf",width=4, height=4)

###########################



df <- read.table("data/results_rax_bme.csv", sep="\t", header=F)
for( i in 1:5){
  df[,i] <- df[,i] - df[,6]
}
colnames(df)=c("OLS", "BME", "BE","FM", "PPLACER", "EST")
head(df)
z=melt(df[,1:5])
z=z[z$variable %in% c("FM","BME","PPLACER"),]
ggplot(z, aes(x=value,color=variable)) + stat_ecdf(geom = "line", pad = FALSE, alpha=0.33) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)+
  xlab("algorithm")+ylab("delta RF distance") + theme_classic()
#ggplot(melt(df[,1:10]), aes(variable, value)) + geom_boxplot()+xlab("algorithm")+ylab("delta RF distance") + labs(title=args[2])
ggsave("figures/bme-in-ds1.pdf", width=4, height=4)

#################################################
acc_sca <- read.table("data/allres_num_gene.csv", header=F)
colnames(acc_sca) <- c("error", "query","numgenes")
head(acc_sca)
q <- ggplot(aes(x=error,color=factor(numgenes)),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
  theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
  scale_colour_manual(name = "", values=cbPalette[c(2,3,4,7,9)])
ggsave("./figures/cdf-nuc-genes.pdf",width=4, height=4)

###########################

acc_sca <- read.table("data/results_apples_extended.csv", header=F)
colnames(acc_sca) <- c("error", "query","obj","distal","pendant")
b <- ggplot(acc_sca, aes(x = error, y = obj))
b + geom_point()+ geom_smooth(method = "lm") + theme_minimal() + scale_y_continuous(trans='log10',limits = c(0.5,1200) )
ggsave("./figures/obj-err-scatter.pdf",width=5, height=5)


###########################

acc_sca <- read.table("data/results_apples_extended.csv", header=F)
colnames(acc_sca) <- c("error", "query", "noname", "branch","obj","distal","pendant", "strategy", "numgenes", "replicate")
head(acc_sca)
acc_sca$d0 <- acc_sca$distal<=0
acc_sca$p0 <- acc_sca$pendant<=0
acc_sca$q0 <- acc_sca$obj<=0

zbn <- read.table("data/dp.txt", header=F)
colnames(zbn) <- c("error", "query", "strategy", "numgenes", "replicate", "dp")
head(zbn)
acc_sca = merge(acc_sca,zbn) 
acc_sca$x <- paste(acc_sca$d0,acc_sca$p0,acc_sca$q0,acc_sca$dp)
acc_sca$x = as.factor(acc_sca$x)

levels(acc_sca$x) = c("1,1", "0,0", "0,1",  "1,0", "0,0", "-1")



length(acc_sca[acc_sca$q0==T,]$error)
#acc_sca$x <- paste(acc_sca$d0,acc_sca$p0)
q <- ggplot(aes(x=error,color=x),data=acc_sca[acc_sca$q0==F & (acc_sca$strategy == "random" | acc_sca$numgenes == 381) ,]) + 
  stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
  theme(legend.position="none", text = element_text(size=16)) + #scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
  scale_colour_manual(name = "DISTAL==0, PENDANT==0", values=cbPalette[c(2,1,3,9,4,7)]) +
  scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) +
  annotate(geom="text",label="d > 0 and p > 0",color=cbPalette[2],x=0.45,y=0.68)+
  annotate(geom="text",label="d = 0 and p > 0",color=cbPalette[9],x=0.45,y=0.52)+
  annotate(geom="text",label="d > 0 and p = 0",color=cbPalette[3],x=0.45,y=0.82)+
  annotate(geom="text",label="d = 0 and p = 0",color=cbPalette[1],x=0.45,y=0.12)+
  guides(color=guide_legend(nrow=3, byrow=T, title.position = "top")) 

ggsave("./figures/triple-combo-ecdf-updated.pdf",width=5, height=5)

wilcox.test(acc_sca[(acc_sca$strategy == "random" | acc_sca$numgenes == 381) & acc_sca$x == "0,0" ,]$error, acc_sca[ (acc_sca$strategy == "random" | acc_sca$numgenes == 381) & !(acc_sca$x == "0,0") ,]$error)
wilcox.test(acc_sca[acc_sca$q0==F & acc_sca$p0==F,]$error,acc_sca[acc_sca$q0==F & acc_sca$p0==T,]$error)
wilcox.test(acc_sca[acc_sca$q0==F & acc_sca$d0==T & acc_sca$p0==F,]$error,acc_sca[acc_sca$q0==F & acc_sca$d0==F & acc_sca$p0==F,]$error)

#length(acc_sca[acc_sca$q0==F & acc_sca$p0==F,]$error)
##################################################################

mp <- read.table("~/Workspace/btol/data/misplacement_pattern.csv", header=F)
colnames(mp) <- c("query", "replicate", "numgenes")
mp$included_mp <- TRUE

require(plyr)
joined_mp <- join_all(list(acc_sca,mp), by = c("query", "replicate", "numgenes"), type = 'full')
head(joined_mp)
joined_mp$included_mp <- ! is.na(joined_mp$included_mp)

joined_mp$x <- paste(joined_mp$d0,joined_mp$p0)

joined_mp$dg02 <- joined_mp$pendant > 0.30
joined_mp$x <- paste(joined_mp$d0,joined_mp$p0)
joined_mp$xalt <- paste(joined_mp$dg02,joined_mp$d0)


q <- ggplot(aes(x=error,linetype=included_mp, color=x),data=joined_mp[joined_mp$q0==F & joined_mp$numgenes == 50,]) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
  theme(legend.position="bottom") + #scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
  scale_colour_manual(name = "distal==0, pendant==0", values=cbPalette[c(2,3,9,1,4,7)]) +
  scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,40)) +
#  annotate(geom="text",label="d > 0 and p > 0",color=cbPalette[2],x=3,y=1)+
#  annotate(geom="text",label="d = 0 and p > 0",color=cbPalette[9],x=3,y=0.82)+
#  annotate(geom="text",label="d > 0 and p = 0",color=cbPalette[3],x=3,y=0.57)+
#  annotate(geom="text",label="d = 0 and p = 0",color=cbPalette[1],x=3,y=0.09)+
 guides(color=guide_legend(nrow=3, byrow=T, title.position = "top"))

z<- joined_mp[joined_mp$error >= 10 & joined_mp$numgenes == 50 & joined_mp$included_mp==F,]
nrow(z)
ggsave("./figures/misplacement_pattern.pdf",width=5, height=5.5)

nrow(joined_mp[joined_mp$p0==T & joined_mp$dg02==T & joined_mp$included_mp == 50,])

#########################################################
acc_sca <- read.table("data/results_apples_extended.csv", header=F)
colnames(acc_sca) <- c("error", "query", "noname", "branch","obj","distal","pendant", "strategy", "numgenes", "replicate")
head(acc_sca)
acc_sca$d0 <- acc_sca$distal<=0
acc_sca$p0 <- acc_sca$pendant<=0
acc_sca$q0 <- acc_sca$obj<=0
acc_sca$x <- paste(acc_sca$d0,acc_sca$p0,acc_sca$q0)


q <- ggplot(aes(x=error,color=qbin),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
  theme(legend.position="bottom") + #scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
  scale_colour_manual(name = "MLSE Error (Q)", values=cbPalette[c(7,2,3,9,1,4)]) +
  scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) +
  guides(color=guide_legend(nrow=3, byrow=T, title.position = "top"))

ggsave("./figures/all-Q-bin.pdf",width=5, height=5.5)



acc_sca$qbin <- cut(acc_sca$obj, c(-Inf, 0, 0.025, 0.05, 0.2, 0.8,  Inf))
summary(acc_sca$obj)

q <- ggplot(aes(x=error,color=qbin),data=acc_sca[ acc_sca$q0 == F,]) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
  theme(legend.position=c(0.75,0.15)) + #scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
  scale_colour_manual(name = "MLSE Error (Q)", values=cbPalette[c(7,2,3,9,1,4)]) +
  scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) +
  guides(color=guide_legend(nrow=3, byrow=T, title.position = "top"))

ggsave("./figures/all-Q-bin-numgenes-gt-40-new.pdf",width=5, height=5)

qt=as.vector(quantile(acc_sca[acc_sca$obj > 0,]$obj,(1:11)/11))
qt[1]=0
qt[11] = qt[11]+1
qt=round(qt,digits=2)

ggplot(acc_sca[acc_sca$obj > 0,], aes(x = cut(obj,breaks = qt), y = error, group = strategy, color=strategy)) + 
  stat_summary(position = position_dodge(width=0.3)) + theme_classic() +  #+ 
  stat_summary(alpha=0.3,  size=0.3, position = position_dodge(width=0.3),fun.min = function(x) quantile(x,0.1),fun.max = function(x) quantile(x,0.90),fun =  function(x) mean(x))+
  stat_summary(geom="line") + labs( x = "MLSE Error (Q)", y = "Error") +
  theme(legend.position=c(0.2,0.8),text = element_text(size=16),axis.text.x = element_text(angle=45, hjust=1)) 

ggsave("./figures/all-Q-bin-strategy-line.pdf",width=5, height=5)

##############################################################

  #stat_smooth() +
  #geom_vline(xintercept = 0.05)+
  #geom_vline(xintercept = 0.2)
  #scale_x_continuous(trans='log2',limits = c(0.0001,20))+
  #scale_y_continuous(trans='sqrt')



ggplot(acc_sca[acc_sca$error>-1,], aes(x = distal, y = pendant-distal, color=cut(error,c(-1,4,30)))) + 
  geom_jitter(height = 0, width = 0.001, alpha=0.5) + theme_classic() 
  
 #ggplot(acc_sca, aes(x = x, y = error)) +1 geom_violin(scale = "count", trim = TRUE) + theme_classic()
 ggplot(acc_sca, aes(x = x, y = error)) + geom_boxplot(outlier.alpha=0.33,alpha=0, varwidth = T ) + 
   #geom_point(position = position_jitterdodge(dodge.width=0.75,jitter.width=0.3),alpha=0.33 )+
   theme_classic()
 
 ggsave("./figures/distal-pendant-correlation-distal-vs-pendant.pdf",width=5, height=5)
 
 #b + geom_point()+ geom_smooth(method = "lm") + theme_minimal() + scale_y_continuous(trans='log10',limits = c(0.5,1200) )

 
 ggplot(acc_sca,aes(x=error, y=x)) +   geom_point(alpha=0.75, size=0.7, position = position_jitter(width = 0, height = 0.1)) +
    geom_boxplot(outlier.alpha=0.33,alpha=0,)
 
 ###########################
 
 acc_sca <- read.table("data/hybrid_vs_mlse.csv", header=F)
 colnames(acc_sca) <- c("error", "query", "criteria")
 q <- ggplot(aes(x=error,color=criteria),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
 q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
   theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
 ggsave("./figures/hybrid_vs_mlse.pdf",width=3.5, height=3.5)
 
 ###########################
 
 acc_sca <- read.table("data/third_vs_nothird.csv", header=F)
 colnames(acc_sca) <- c("error", "query", "thirdcodon")
 q <- ggplot(aes(x=error,color=thirdcodon),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
 q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
   theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
 
 ggsave("./figures/third_vs_nothird.pdf",width=3.5, height=3.5)

 ###########################
 
 acc_sca <- read.table("data/neg_vs_noneg.csv", header=F)
 colnames(acc_sca) <- c("error", "query", "negativity")
 q <- ggplot(aes(x=error,color=negativity),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
 q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
   theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
 
 ggsave("./figures/neg_vs_noneg.pdf",width=3.5, height=3.5)

 
 ###########################
 
 acc_sca <- read.table("data/03_vs_infty.csv", header=F)
 colnames(acc_sca) <- c("error", "query", "threshold")
 q <- ggplot(aes(x=error,color=threshold),data=acc_sca[acc_sca$threshold=="0.3-200",]) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
 q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
   theme(legend.position="bottom") + ylim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),limits = c(0.5,1))
 
 ggsave("./figures/03.pdf",width=3.5, height=3.5)     
 
 
 #######################
 
 df <- read.table("data/thr_exp_results.csv", header=F)
 x <- df[,c(1,3,4)]
 head(x)
 colnames(x) <- c("error", "filterThr", "baseThr")
 aggdataAcc <- aggregate(x, by=list(x$filterThr, x$baseThr), FUN= function(c)sum(c==0)/sum(c>=0)) 
 aggdataAvg <- aggregate(x, by=list(x$filterThr, x$baseThr), FUN= mean) 
 
 agg <- merge(aggdataAcc, aggdataAvg, by=c("Group.1", "Group.2"))
 
 agg <- agg[,c(1,2,3,6)]
 colnames(agg) <- c("filterThr", "baseThr", "inaccuracy", "meanerror")
 agg$inaccuracy <- 100* agg$inaccuracy
 head(agg)
 
 
 p <- ggplot(data=agg, aes(x=meanerror, y=inaccuracy, colour = as.factor(filterThr), shape=as.factor(baseThr), fill=as.factor(filterThr)))
 p + geom_point(size=2.0,alpha = 0.5)  + theme_classic() + theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1)) +
   scale_colour_manual(name = "", values=cbPalette[c(4,2,9,1,5,6,8)]) +
   scale_shape_manual(name = "", values = c(21, 24, 22, 4, 25, 23)) + 
   scale_fill_manual(na.value=NA, guide="none",values=cbPalette[c(4,2,9,1,5,6,8)]) +
   scale_y_continuous( breaks = c(40,50,60,70,80,90,100), limits = c(60,80))+
   xlab("Average delta error (edges)") + ylab("Placement Accuracy")

 ggsave("./figures/thr-exp.pdf",width=6.0, height=6.0)  
 
 ####################### 
 
 df <- read.table("data/aln_vs_est.csv", header=F)
 x <- df[,c(1,3,4,5)]
 head(x)
 colnames(x) <- c("error", "filterThr", "baseThr","method")
 
 aggdataAcc <- aggregate(x, by=list(x$filterThr, x$baseThr, x$method), FUN= function(c)sum(c==0)/sum(c>=0)) 
 aggdataAvg <- aggregate(x, by=list(x$filterThr, x$baseThr, x$method), FUN= mean)
 agg <- merge(aggdataAcc, aggdataAvg, by=c("Group.1", "Group.2", "Group.3"))
 agg <- agg[agg$Group.2 != 500,]

  agg$threshold <- paste(agg$Group.1, agg$Group.2)
  res <- agg[c("Group.3","error.x","error.y", "threshold")]  
  colnames(res) <- c("alignment", "accuracy", "meanerror", "threshold") 
  res$accuracy <- 100* res$accuracy
  head(res)
  p <- ggplot(data=res, aes(x=meanerror, y=accuracy, colour = as.factor(alignment), shape=as.factor(threshold), fill=as.factor(alignment)))
  p + geom_point(size=2.0,alpha = 0.5)  + theme_classic() + theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_colour_manual(name = "", values=cbPalette[c(7,4,2,9,1,5,6,8)]) +
    scale_shape_manual(name = "", values = c(21, 24, 22, 4, 25, 23)) + 
    scale_fill_manual(na.value=NA, guide="none",values=cbPalette[c(7,4,2,9,1,5,6,8)]) +
    scale_y_continuous( breaks = c(40,50,60,70,80,90,100), limits = c(40,100))+
    xlab("Average delta error (edges)") + ylab("Placement Accuracy")
  
  ggsave("./figures/est-vs-true-aln.pdf",width=6.0, height=6.0)
  
  ####################### 
  
  df <- read.table("data/prot_vs_nuc.csv", header=F)
  x <- df[,c(1,3)]
  colnames(x) <- c("error", "molecule")
  head(x)
  x$error <- as.numeric(x$error)
  aggdataAcc <- aggregate(x, by=list(x$molecule), FUN= function(c)sum(c==0)/sum(c>=0)) 
  aggdataAvg <- aggregate(x, by=list(x$molecule), FUN= mean) 
  
  agg_first <- merge(aggdataAcc, aggdataAvg, by=c("Group.1"))
head(agg_first)
agg_first <- agg_first[,c(1,2,4)]
  colnames(agg_first) <- c("molecule", "inaccuracy", "meanerror")
  agg_first$inaccuracy <- 100* agg_first$inaccuracy
  head(agg_first)
  
  df <- read.table("data/aln_vs_est.csv", header=F)
  x <- df[,c(1,3,4,5)]
  head(x)
  colnames(x) <- c("error", "filterThr", "baseThr","method")
  
  aggdataAcc <- aggregate(x, by=list(x$filterThr, x$baseThr, x$method), FUN= function(c)sum(c==0)/sum(c>=0)) 
  aggdataAvg <- aggregate(x, by=list(x$filterThr, x$baseThr, x$method), FUN= mean)
  agg <- merge(aggdataAcc, aggdataAvg, by=c("Group.1", "Group.2", "Group.3"))
  head(agg)
  agg <- agg[agg$Group.2 == 50 & agg$Group.1 == 0.2 & agg$Group.3 == "true",]
  
  agg <- agg[,c(3,4,8)]
  colnames(agg) <- c("molecule", "inaccuracy", "meanerror")
  agg$inaccuracy <- 100* agg$inaccuracy
  
  y=rbind(agg,agg_first)
  agg=y
  head(agg)
  
   
  p <- ggplot(data=agg, aes(x=meanerror, y=inaccuracy, colour = molecule, fill=molecule))
  p + geom_point(size=2.0,alpha = 0.5)  + theme_classic() + theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_colour_manual(name = "", values=cbPalette[c(4,2,9,1,5,6,8)]) +
    #scale_shape_manual(name = "", values = c(21, 24, 22, 4, 25, 23)) + 
    scale_fill_manual(na.value=NA, guide="none",values=cbPalette[c(4,2,9,1,5,6,8)]) +
    scale_y_continuous( breaks = c(40,50,60,70,80,90,100), limits = c(60,80))+
    xlab("Average delta error (edges)") + ylab("Placement Accuracy") + xlim(0.5,2)
  
  ggsave("./figures/prot-vs-nuc.pdf",width=6.0, height=6.0)
  
  
  ######################
  
  
  df <- read.table("data/tc_parameter_compare.csv", header=F)
  x <- df[,c(1,3)]
  colnames(x) <- c("error", "tcparam")
  head(x)
  x$error <- as.numeric(x$error)
  aggdataAcc <- aggregate(x, by=list(x$tcparam), FUN= function(c)sum(c==0)/sum(c>=0)) 
  aggdataAvg <- aggregate(x, by=list(x$tcparam), FUN= mean) 
  
  agg_first <- merge(aggdataAcc, aggdataAvg, by=c("Group.1"))
  head(agg_first)
  agg_first <- agg_first[,c(1,2,4)]
  colnames(agg_first) <- c("threshold", "accuracy", "meanerror")
 #agg_first$inaccuracy <- 100* agg_first$inaccuracy
  head(agg_first)
  
  p <- ggplot(data=agg_first, aes(x=meanerror, y=accuracy, colour = as.factor(threshold), shape=as.factor(threshold), fill=as.factor(threshold)))
  p + geom_point(size=4.0,alpha = 1)  + theme_classic() + 
    theme(legend.position = c(0.2,0.8) ,
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_colour_manual(name = "Multiplier", values=cbPalette[c(4,2,9,1,5,6,8)]) +
    scale_shape_manual(name = "Multiplier", values = c(21, 24, 22, 4, 25, 23)) + 
    scale_fill_manual(name="Multiplier",values=cbPalette[c(4,2,9,1,5,6,8)]) + coord_cartesian(xlim = c(0.47,0.65), ylim=c(0.6,1))+
    #scale_y_continuous( breaks = c(40,50,60,70,80,90,100))+
    xlab("Average delta error (edges)") + ylab("Placement Accuracy")
  
  ggsave("./figures/treecluster_compare.pdf",width=5.0, height=5.0)  
  
  ########
  
  acc_sca <- read.table("data/tc_on_vs_off.csv", header=F)
  colnames(acc_sca) <- c("error", "query","treecluster")
  head(acc_sca)
  q <- ggplot(aes(x=error,color=treecluster),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom") + xlim(0, 25) +  scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
    scale_colour_manual(name = "", values=cbPalette[c(2,3,4,9)])
  ggsave("./figures/tc_on_off.pdf",width=4, height=4)
  
  
  ########
  
  df <- read.table("data/all_jobsizes.txt", sep="\t", header=F)
  head(df)
  colnames(df) <- c("threshold", "parname", "gene", "size")
  

  # Histogram with density plot
  ggplot(df, aes(x=size)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins=50) +
    geom_density(alpha=.2, fill="#FF6666") + facet_grid(.~threshold) + theme_bw()
  ggsave("./figures/jobsizes.pdf",width=10, height=4)

  ggplot(df, aes(x=size)) + 
    geom_histogram(colour="black", fill="white", bins=50) +
    geom_density(alpha=.2, fill="#FF6666") + facet_grid(threshold~.) + theme_bw()
  ggsave("./figures/jobsizes_density.pdf",width=4, height=10)  
  
  
  ##########################
  
  acc_sca <- read.table("data/numgenes.csv", header=F)
  colnames(acc_sca) <- c("error", "query","numgenes")
  head(acc_sca)
  q <- ggplot(aes(x=error,color=factor(numgenes)),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+ 
    scale_colour_manual(name = "", values=cbPalette[c(2,3,4,7,9)])
  ggsave("./figures/cdf-nuc-genes-50-400.pdf",width=4, height=4)
  
  
  ########################
  setwd("~/Workspace/btol/")
  
  acc_sca <- read.table("data/randomtest.csv", header=F)
  head(acc_sca)
  colnames(acc_sca) <- c("error", "query", "replicate","numgenes", "selection","numcopy")
  #acc_sca$selection=as.factor(acc_sca$selection)
  #levels(acc_sca$selection) <- c("best","random")
  rnds = acc_sca[acc_sca$selection != "best",]
  q <- ggplot(aes(x=error, color=as.factor(numgenes)),data=rnds) + stat_ecdf(geom = "line", pad = FALSE) + 
    stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
  #geom_linerange(aes(ymin=0.3,ymax=1.0, x=mean(error), color=method),linetype=1,size=0.5)
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom") + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
    scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) +
    scale_colour_manual(name = "", values=cbPalette[c(2,3,4,9,5,6,1)]) + 
    scale_linetype_manual(name="",values=c(1,2,1,3)) + #guides(color=guide_legend(nrow=1, byrow=T)) +
    guides(color=guide_legend(nrow=2, byrow=T), linetype=guide_legend(nrow=2, byrow=T)) +
    scale_shape_manual(name="", values = c(19, 17, 15, 4)) #+ theme(legend.position="bottom", legend.box="vertical")
  
  ggsave("./figures/randomtest-randomonly-cdf.pdf",width=5, height=5.5)
  
  wilcox.test(acc_sca[acc_sca$numgenes==50 & acc_sca$selection == "random",]$error,acc_sca[acc_sca$numgenes==50 & acc_sca$selection == "best",]$error)
  
  #########################################################
  setwd("~/Workspace/btol/")
  
  acc_sca <- read.table("data/randomtest.csv", header=F)
  head(acc_sca)
  colnames(acc_sca) <- c("error", "query", "replicate","numgenes", "selection","numcopy")
  #acc_sca$selection=as.factor(acc_sca$selection)
  #levels(acc_sca$selection) <- c("best","random")
  cmp = acc_sca[acc_sca$numgenes %in% c(10,25,50,381),]
  q <- ggplot(aes(x=error, color=as.factor(numgenes),shape=selection,linetype=selection),data=cmp) + stat_ecdf(geom = "line", pad = FALSE) + 
    stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
  #geom_linerange(aes(ymin=0.3,ymax=1.0, x=mean(error), color=method),linetype=1,size=0.5)
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom") + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
    scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) +
    scale_colour_manual(name = "", values=cbPalette[c(2,3,4,9,5,6,1)]) + 
    scale_linetype_manual(name="",values=c(1,2,1,3)) + #guides(color=guide_legend(nrow=1, byrow=T)) +
    guides(color=guide_legend(nrow=2, byrow=T), linetype=guide_legend(nrow=2, byrow=T)) +
    scale_shape_manual(name="", values = c(19, 17, 15, 4)) #+ theme(legend.position="bottom", legend.box="vertical")
  
  ggsave("./figures/randomtest-random-vs-best.pdf",width=5, height=5.5)
  
  
  wilcox.test(cmp[cmp$selection=="random" & cmp$numgenes == 25,]$error,cmp[cmp$selection=="best" & cmp$numgenes == 10,]$error)
  mean(cmp[cmp$selection=="random" & cmp$numgenes == 25,]$error)
  mean(cmp[cmp$selection=="best" & cmp$numgenes == 50,]$error)-mean(cmp[cmp$selection=="random" & cmp$numgenes == 50,]$error)
  sum(cmp[cmp$selection=="best" & cmp$numgenes == 381,]$error<=2)
  ####################traj1##############################
  
  acc_sca <- read.table("data/randomtest.csv", header=F)
  head(acc_sca)
  colnames(acc_sca) <- c("error", "query","replicate", "totgenes" , "selection", "numgenes")
  #acc_sca$selection=as.factor(acc_sca$selection)
  #levels(acc_sca$selection) <- c("best","random")
  
  
  
  
  
  
  x = acc_sca[acc_sca$numgenes >= 50 & acc_sca$error >=12 ,]
  prob_queries = x$query
  print(prob_queries)
  df_prob = acc_sca[acc_sca$query %in% prob_queries,]
  
  p <- ggplot(data=df_prob[df_prob$numgenes > 20,], aes(x=numgenes, y=error, color=as.factor(replicate))) + facet_wrap(vars(query),nrow = 5)
  p + geom_point(size=2) + geom_line(size=1) + scale_x_continuous(breaks=c(0,10,25,50,75,100,150,200)) 
  
  ggsave("./figures/compare-replicate-facet-query.pdf",width=12, height=16)
  
  print(df_prob)
  
  ######################traj2###########################
  
  acc_sca <- read.table("data/randomtestnumcopy.csv", header=F)
  head(acc_sca)
  colnames(acc_sca) <- c("error", "query","replicate","totgenes", "numcopy")
  acc_sca$numcopy=as.factor(acc_sca$numcopy)
  #levels(acc_sca$selection) <- c("best","random")
  #mean(acc_sca[acc_sca$totgenes==50 ,]$error)
  
  x = acc_sca[acc_sca$numcopy >= 20 & acc_sca$error >=15 ,]
  prob_queries = x$query
  print(prob_queries)
  df_prob = acc_sca[acc_sca$query %in% prob_queries,]
  
  p <- ggplot(data=df_prob, aes(x=numcopy, y=error, color=as.factor(totgenes)))
  p + geom_point(size=2) + geom_smooth(method="lm")+ scale_x_continuous(breaks=c(0,10,25,50,75,100,150,200))  + 
    theme_classic()
  
  ggsave("./figures/badqueries-trend-pertotal.pdf",width=6, height=6)
  

  x = acc_sca[acc_sca$numcopy >= 20 & acc_sca$error >=15 ,]
  prob_queries = x$query
  print(prob_queries)
  df_prob = acc_sca[acc_sca$query %in% prob_queries,]
  p <- ggplot(data=df_prob, aes(x=numcopy, y=error, color=query))
  p + geom_point(size=2) + geom_line(size=1) + scale_x_continuous(breaks=c(0,10,25,50,75,100,150,200))  + 
    theme_classic()
  
  ggsave("./figures/badqueries-trend-pertotal.pdf",width=6, height=6)
  
  ######################################################
  acc_sca <- read.table("data/randomtestnumcopy.csv", header=F, sep = "\t")
  head(acc_sca)
  colnames(acc_sca) <- c("error", "query","replicate" ,"totgenes", "numcopy")
  
  acc_sca$bin <- cut(acc_sca$numcopy, c(0, 20, 40, 80, 120, 150))
  head(acc_sca)
  
  ggplot(acc_sca, aes(bin, error)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5, alpha=0.5)+
    scale_y_continuous(trans='sqrt', breaks = c(0, 4, 16, 64)) + theme_classic()
  ggsave("./figures/numcopy-jitter.pdf",width=6, height=6)
  
  ###################################################
  
  rmd <- acc_sca[acc_sca$bin != "(0,5]"  & acc_sca$bin != "(5,10]", ]
  
  q <- ggplot(aes(x=error, color=bin),data=rmd) + stat_ecdf(geom = "line", pad = FALSE) + 
    stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
  #geom_linerange(aes(ymin=0.3,ymax=1.0, x=mean(error), color=method),linetype=1,size=0.5)
  
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom")  + 
    scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) 

  ggsave("./figures/numcopy-rmd.pdf",width=6, height=6)
  
  q <- ggplot(aes(x=error, color=bin),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + 
    stat_ecdf(geom = "point", size=0.75, pad = FALSE) 
  #geom_linerange(aes(ymin=0.3,ymax=1.0, x=mean(error), color=method),linetype=1,size=0.5)
  
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom") +
    scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) +
    scale_colour_manual(name = "", values=cbPalette[c(7,4,2,9,1,5,6,8)]) +
  
  
  ggsave("./figures/numcopy-all.pdf",width=6, height=6)
  
  
  ###################################################
  
  
  acc_sca$bin2 <- cut(acc_sca$numcopy, c(0, 20, 40, 80, 150))
  print(acc_sca[acc_sca$totgenes==150 & acc_sca$bin2 == "(20,40]",])
  
  q <- ggplot(aes(x=error, color=as.factor(totgenes)),data=acc_sca[acc_sca$bin2 != "(0,20]",]) + stat_ecdf(geom = "line", pad = FALSE) + 
    stat_ecdf(geom = "point", size=0.75, pad = FALSE) + facet_grid(.~bin2)
  #geom_linerange(aes(ymin=0.3,ymax=1.0, x=mean(error), color=method),linetype=1,size=0.5)
  
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom") +
    scale_x_continuous(trans='sqrt', breaks = c(0,1,4,9,16,25,36,49)) + coord_cartesian(c(0,20)) 
  
  
  ggsave("./figures/numcopy-ratio.pdf",width=9, height=4)

  ############################################
  
  df <- read.table("~/Workspace/btol/data/domain_marker_counts.txt", header=F)
  colnames(df) <- c("gene", "bcount", "acount")
  df$gene <- factor(df$gene, levels = df$gene)
  
  b <- ggplot(df, aes(x = gene, y = acount))
  b + geom_point() + theme_classic() + 
    theme(text = element_text(size=6), axis.text.x = element_text(angle=90, hjust=1)) +
    geom_vline(xintercept = 50, colour="#22BB00", linetype="dashed" ) +
    geom_vline(xintercept = 70, colour="#22BBBB", linetype="dashed" ) +
    geom_vline(xintercept = 85, colour="#BB6666", linetype="dashed" ) +
    labs( x = "Genes", y = "Archaea count")

  ggsave("./figures/domain_counts.pdf",width=16, height=4)
  
  ############################################
  head(df)
  
  summ <- 0.0
  nroww <- 0.0
  for (row in 1:nrow(df)) {
    nroww = nroww + 1
     summ = summ + df[row, "acount"]
     df[row, "archacc"] <- summ/669
     df[row, "archavgacc"] <- summ/669/nroww
  }
  head(df)
  
  ggplot(df[c(1:100),]) + geom_path( aes(x = gene, y = archavgacc, group=2))  + theme_classic() + 
    geom_path(aes(x = gene, y = archacc/100, group=1), color=cbPalette[4])+
    geom_vline(xintercept = 68, colour="#BB6666", linetype="dashed" ) +    
    theme(text = element_text(size=8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y.right = element_text(color = cbPalette[4])) + 
    scale_y_continuous("Average Archaea Occupation", sec.axis = sec_axis(~.*100, name="Cumulative Archaea Occupation"))
    
  ggsave("./figures/marker_set_choice.pdf",width=7, height=5)
  
  
  #############################################
  
  
  acc_sca <- read.table("data/normal_vs_tapered.csv", header=F)
  head(acc_sca)
  colnames(acc_sca) <- c("x","error","y", "filtering")
  q <- ggplot(aes(x=error, color=filtering, linetype=filtering, shape=filtering),data=acc_sca) + stat_ecdf(geom = "line", pad = FALSE) + stat_ecdf(geom = "point", size=0.75, pad = FALSE)
  q + theme_classic() + labs( x = "Delta Error (edges)", y = "Empirical CDF") + 
    theme(legend.position="bottom") + xlim(0, 25) + scale_y_continuous(breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

  ggsave("./figures/normal-vs-tapered.pdf",width=4.5, height=4.5)
  
  
  #################
  df <- read.table("data/dists_5000.txt", header=F)
  df2 <- read.table("data/dists_rand.txt", header=F)
  df$method <- "treecluster"
  df2$method <- "random"
  dfx<-rbind(df,df2)
  ggplot() + geom_density(data=dfx,aes(V1, linetype=V2, color=method)) + #geom_density(data=df, color="blue",aes(V1,linetype=V2)) +
    theme_classic()  + coord_cartesian(ylim=c(0,1.5*10^-7))
  ggsave("./figures/taxon_subsampling.pdf",width=5.5, height=4.5)
  
  
  df3 <- read.table("data/together-query-len.txt", header=F)
head(df3)  
ggplot() + geom_density(data=df3,aes(2*V2, color=V4, linetype=as.factor(V3))) + #geom_density(data=df, color="blue",aes(V1,linetype=V2)) +
  theme_classic() +scale_x_log10() + geom_vline(xintercept = 1) #+ coord_cartesian(ylim=c(0,1.5*10^-7))

ggsave("./figures/taxon_subsampling_density.pdf",width=7, height=5)

ggplot() + stat_ecdf(data=df3,aes(2*V2, color=V4, linetype=as.factor(V3))) + #geom_density(data=df, color="blue",aes(V1,linetype=V2)) +
  theme_classic() +scale_x_log10() + geom_vline(xintercept = 1) #+ coord_cartesian(ylim=c(0,1.5*10^-7))

ggsave("./figures/taxon_subsampling_ecdf.pdf",width=7, height=5)

###############################

df <- read.table("data/report.txt")
head(df)
ggplot() + geom_bar(aes(x=df$V1)) + theme_classic() + xlab("Occupancy")
ggsave("./figures/occupancy_200k.pdf",width=5, height=5)

###########################################################

df <- read.table("data/backbone_discordance_06.csv")
head(df)
colnames(df) <- c( "uDance", "FastTree2", "true", "truefull", "rep")
dfm <- melt(data=df, id=c("rep"))
dfm$variable=as.factor(dfm$variable)
dfm$variable = factor(dfm$variable, levels=c("truefull", "true", "FastTree2",  "uDance"))
#levels(dfm$variable) = c("true", "uDance", "FastTree2")
ggplot(dfm, aes(x = variable, y = value, fill=variable)) + 
  geom_violin(scale = "count", trim = TRUE ) + 
  geom_boxplot(width = 0.2, outlier.alpha = 0)+theme_classic() +
  scale_fill_manual(values = cbPalette[c(2,3,4,6)])

ggsave("./figures/violin_backbone_udance_vs_fasttree.pdf",width=5, height=5)

#########################################################
df <-read.table("data/contigcount.csv",sep=",")
colnames(df) <- c("numcontigs", "backbone")
ggplot(aes(x=numcontigs,color=backbone),data=df) + stat_ecdf(geom = "line", pad = FALSE) + geom_hline(yintercept = 0.95) +
  theme_classic() + labs( x = "Number of contigs", y = "Empirical CDF") + scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 0.95, 1))+
  theme(legend.position="bottom") + coord_cartesian(c(0, 500))

ggsave("./figures/contigcount.pdf",width=5, height=5)


df <-read.table("data/guncscores.csv",sep=",")
colnames(df) <- c("guncscore", "backbone")
ggplot(aes(x=guncscore,color=backbone),data=df) + stat_ecdf(geom = "line", pad = FALSE) + geom_vline(xintercept = 0.45) + geom_vline(xintercept = 0.95, color="red") +
  theme_classic() + labs( x = "GUNC score", y = "Empirical CDF") + scale_x_continuous(breaks=c(0,0.25, 0.45, 0.5, 0.75, 0.95, 1))+
  theme(legend.position="bottom") + coord_cartesian(c(0, 1))

ggsave("./figures/guncscores.pdf",width=5, height=5)

df <-read.table("data/chim_scores.csv",sep="\t")
colnames(df) <- c( "species","guncscore", "backbone")
ggplot(aes(x=guncscore,color=backbone),data=df) + stat_ecdf(geom = "line", pad = FALSE) + geom_vline(xintercept = 0.45) + geom_vline(xintercept = 0.95, color="red") +
  theme_classic() + labs( x = "GUNC score", y = "Empirical CDF") + scale_x_continuous(breaks=c(0,0.25, 0.45, 0.5, 0.75, 0.95, 1))+
  theme(legend.position="bottom") + coord_cartesian(c(0, 1))

ggsave("./figures/guncscores_updated.pdf",width=5, height=5)


df <-read.table("data/sing_chim_scores.csv",sep="\t")
colnames(df) <- c( "species","guncscore", "ANI", "backbone")
ggplot(aes(x=guncscore,color=backbone),data=df) + stat_ecdf(geom = "line", pad = FALSE) + geom_vline(xintercept = 0.45) + geom_vline(xintercept = 0.95, color="red") +
  theme_classic() + labs( x = "GUNC score", y = "Empirical CDF") + scale_x_continuous(breaks=c(0,0.25, 0.45, 0.5, 0.75, 0.95, 1))+
  theme(legend.position="bottom") + coord_cartesian(c(0, 1))


ggsave("./figures/gunc_ccs_sing.pdf",width=5, height=5)

ggplot(aes(x=ANI,color=backbone),data=df) + stat_ecdf(geom = "line", pad = FALSE)  + geom_vline(xintercept = 0.95, color="red") +
  theme_classic() + labs( x = "ANI score", y = "Empirical CDF") + scale_x_continuous(breaks=c(0,0.25, 0.45, 0.5, 0.75, 0.95, 1))+
  theme(legend.position="bottom") + coord_cartesian(c(0, 1))

ggsave("./figures/gunc_rrs_sing.pdf",width=5, height=5)

##########################
df <- read.table("data/gcorr.csv")

ggplot(aes(x=V1, y=V2), data=df[df$V1 < 0.45,]) + geom_point() + theme_bw() + geom_smooth(method="lm") + labs( x = "CCS", y = "RRS")

ggsave("./figures/gunc_corr.png",width=5, height=5)

ggplot(aes(x=V1, y=V2), data=df[df$V1 < 0.45 & df$V1 > 0,]) + geom_point() + theme_bw() + geom_smooth(method="lm") + labs( x = "CCS", y = "RRS")

ggsave("./figures/gunc_corr_nonzero.png",width=5, height=5)

head(df)



#####################
require(ggrepel)
df <- read.table("data/qc_comp.tsv")
head(df)
colnames(df) <- c("tag", "gtdbfail", "in10K", "CSS", "contamination", "dbIdentity")

ggplot(aes(x=CSS, y=contamination, color=gtdbfail, shape=in10K), data=df) + geom_point() + theme_bw()+
  geom_vline(xintercept = 0.45, color="red")

ggplot(aes(x=dbIdentity, y=contamination, color=gtdbfail, shape=in10K), data=df[df$CSS >=.45, ]) + geom_point(size=2) + theme_bw() +
  geom_text_repel(
    data = subset(df, tag %in% c("G008349105","G001068645","G001059375","G001055105","G900143315","G002006985","G000780515","G011754595","G900034435","G000508745","G009805725","G001228205","G001394855","G001318425","G001397415","G001394035","G001397415","G003049385","G001076125")),
    aes(label = tag),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +   
  scale_colour_manual(values=cbPalette[c(2,3,6,4)]) #+ coord_cartesian(ylim = c(0,1.1)) 
  


ggsave("./figures/chim_cont_vs_nov.pdf",width=7, height=6)

ggplot(aes(x=CSS, y=contamination, color=gtdbfail, shape=in10K), data=df[df$CSS >=.45, ]) + geom_point(size=2) + theme_bw() +
  geom_text_repel(
    data = subset(df, tag %in% c("G008349105","G001068645","G001059375","G001055105","G900143315","G002006985","G000780515","G011754595","G900034435","G000508745","G009805725","G001228205","G001394855","G001318425","G001397415","G001394035","G001397415","G003049385","G001076125")),
    aes(label = tag),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +   
  scale_colour_manual(values=cbPalette[c(2,3,6,4)]) #+ coord_cartesian(xlim = c(0,1.1)) 

ggsave("./figures/chim_cont_vs_css.pdf",width=7, height=6)

##############################
df <- read.table("data/qc_comp_n50.tsv")
head(df)
colnames(df) <- c("tag", "gtdbfail", "in10K", "CSS", "contamination", "dbIdentity", "length", "numcontig", "n50", "ambiguous")

ggplot(aes(x=1/numcontig, y=n50, color=gtdbfail, shape=in10K), data=df) + geom_point(size=2) + theme_bw() +
  geom_text_repel(
    data = subset(df, tag %in% c("G008349105","G001068645","G001059375","G001055105","G900143315","G000780515","G900034435","G000508745","G009805725","G001228205","G001394855","G001318425","G001397415","G001394035","G001397415","G003049385","G001076125")),
    aes(label = tag),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +   
  scale_colour_manual(values=cbPalette[c(2,3,6,4)]) + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10")+
  geom_vline(xintercept = 0.001, colour="#BB6666", linetype="dashed")+
  geom_hline(yintercept = 5000, colour="#BB6666", linetype="dashed")
#+ coord_cartesian(ylim = c(0,1.1)) 

ggsave("./figures/n50_numcontig_16k.pdf",width=8, height=7)

ggplot(aes(x=1/numcontig, y=n50, color=gtdbfail, shape=in10K), data=df[!(df$CSS >=.5 & df$contamination >= .25),]) + geom_point(size=2) + theme_bw() +
  #geom_text_repel(
  #  data = subset(df, tag %in% c("G008349105","G001068645","G001059375","G001055105","G900143315","G000780515","G900034435","G000508745","G009805725","G001228205","G001394855","G001318425","G001397415","G001394035","G001397415","G003049385","G001076125")),
  #  aes(label = tag),
  #  size = 2,
  #  box.padding = unit(0.35, "lines"),
  #  point.padding = unit(0.3, "lines")
  #) +   
  scale_colour_manual(values=cbPalette[c(2,3,6,4)]) + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
                                              #geom_vline(xintercept = 0.002, colour="#BB6666", linetype="dashed")+
  geom_hline(yintercept = 3000, colour="#BB6666", linetype="dashed")
#+ coord_cartesian(ylim = c(0,1.1)) 

ggsave("./figures/n50_numcontig_16k_nochim.pdf",width=8, height=7)

ggplot(aes(x=length, y=n50), data=df[df$numcontig>100,]) + 
  geom_point(aes(color=gtdbfail, shape=in10K), size=2) + geom_smooth(method = "lm") + theme_bw() + scale_colour_manual(values=cbPalette[c(2,3,6,4)]) #+
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10")

ggsave("./figures/n50_length_corr.pdf",width=8, height=7)


  ggplot(aes(x=1/numcontig, y=n50/length, color=gtdbfail, shape=in10K), data=df[!(df$CSS >=.5 & df$contamination >= .25),]) + geom_point(size=2) + theme_bw() +
    scale_colour_manual(values=cbPalette[c(2,3,6,4)]) + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
    #geom_vline(xintercept = 0.002, colour="#BB6666", linetype="dashed")+
    geom_hline(yintercept = 1/1500, colour="#BB6666", linetype="dashed")
  
ggsave("./figures/n50_numcontig_16k_nofrag.pdf",width=8, height=7)

##############################

df <- read.table("data/support_comp.csv")
colnames(df) <- c("support", "tree")
ggplot(aes(x=support,color=tree),data=df) + stat_ecdf(geom = "line", pad = FALSE) +
  theme_classic() + labs( x = "Support", y = "Empirical CDF") + scale_x_continuous(breaks=c(0,0.25, 0.45, 0.5, 0.75, 0.95, 1))+
  theme(legend.position="bottom") + coord_cartesian(c(0, 1))

ggsave("./figures/support_16K_comp.pdf",width=5, height=4)

ggplot(aes(x=support,color=tree,),data=df[df$support <.95,]) + stat_ecdf(geom = "line", pad = FALSE) +
  theme_classic() + labs( x = "Support", y = "Empirical CDF") + scale_x_continuous(breaks=c(0,0.25, 0.5, 0.75, 0.95, 1))+
  theme(legend.position="bottom") + coord_cartesian(c(0, 1))

ggsave("./figures/support_16K_comp_lt_0.95.pdf",width=5, height=4)


#####################

library(dplyr)
df <-read.table("./data/snakemake_timeline_100_10.csv", header=F)
colnames(df) <- c("type", "rule", "time", "jobid", "cpus", "act_cpus")
dfun <- df %>% group_by(time) %>% top_n(1, act_cpus)
#dfun <- distinct(df,time, .keep_all= TRUE)
#dfun <- df
ggplot(dfun, aes(x = time/60, y = act_cpus+1, group=1) ) + theme_classic()  + geom_path(color=cbPalette[4]) +
  scale_y_continuous(trans='log2', breaks=c(1,4, 16,64, 256, 1024, 2048)) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1)) +
  labs( x = "Time (minutes)", y = "CPU usage")


ggsave("./figures/udance_timeline_500.pdf",width=3, height=3)

ggplot(dfun, aes(x = time/60, y = act_cpus+1, group=1) ) + theme_classic()  + geom_path(color=cbPalette[4]) +
  #scale_y_continuous(trans='log2', breaks=c(1,4, 16,64, 256, 1024, 2048)) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1)) +
  labs( x = "Time (minutes)", y = "CPU usage")

ggsave("./figures/udance_timeline_500_lin.pdf",width=3, height=3)

#####################

df <- read.table("./data/score_vs_acc.csv", header=F, fill= TRUE)
colnames(df) <- c("qsUpdates", "qsIncremental", "rfUpdates", "rfIncremental")
head(df)
ggplot(aes(y=rfUpdates-rfIncremental, x=(qsUpdates-qsIncremental)/qsIncremental), data=df) + theme_classic()  + 
  geom_point() + geom_smooth(method = "lm") + geom_hline(yintercept = 0.0, colour="#BB6666", linetype="dashed")+
  geom_vline(xintercept = 0.0, colour="#BB6666", linetype="dashed")

ggsave("./figures/hgt_100_lin_corr.pdf",width=4, height=4)

#################################
df <- read.table("./data/gtee_4modelcond_stgt.csv", header=F)
colnames(df) <- c("rf", "hghs", "alnlenpar", "rep", "gene")
df$mc = paste(df$hghs, df$alnlenpar)
head(df)
ggplot(df, aes(x = mc, y = rf, fill=mc)) + 
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + 
  stat_summary()+theme_classic() +
  scale_fill_manual(values = cbPalette[c(2,3,4,7)]) + xlab("Gene Tree Discordance")+ylab("nRF")


ggsave("./figures/four_model_condition_hgt.pdf",width=5, height=5)

df <- read.table("./data/gtee_4modelcond_gtee_extended.csv", header=F)
colnames(df) <- c("rf", "hghs", "alnlenpar", "rep", "gene")
df$mc = paste(df$hghs, df$alnlenpar)
head(df)
ggplot(df, aes(x = mc, y = rf, fill=mc)) + 
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + 
  stat_summary()+theme_classic() +
  scale_fill_manual(values = cbPalette[c(2,3,4,7)]) + xlab("Gene Tree Estimation Error")+ylab("nRF")

ggsave("./figures/four_model_condition_hgt_gtee_ext.pdf",width=5, height=5)

#################################

df <- read.table("./data/astral_u_support_corr.txt")
colnames(df) <- c("supp", "bl")
#ggplot(aes(y=supp,x=bl), data=df) + geom_point()
qnt= quantile(df$bl,c(0:20)/20)
qnt = as.vector(qnt)
qnt[1] = 0
qnt_fixed = qnt 
ggplot(aes(x=cut(bl,qnt_fixed), y = supp),data=df) + 
  theme_classic() + labs( y = "Support", x = "Branch Length") + 
  #geom_vline(aes( xintercept=., color=method,linetype=datatype) , alpha=0.5,
  #           data = dcast(method+datatype+size~.,data=acc_sca[acc_sca$numcopy >= 20,],value.var = "error",fun.aggregate = mean))+
  #stat_summary(aes(y=error),orientation = "x")+
  #geom_vline(aes(xintercept=mean(error),color=method,linetype=datatype), alpha=0.3 , size=0.5) +
  #geom_boxplot(aes(group=numcopy.x), outlier.alpha=0.5, outlier.size = 0.75) + 
  stat_smooth(se=F, color="red") + #facet_wrap(.~datatype, scales = "free") +
  stat_summary(alpha=0.8,  size=0.3, fun.min = function(x) quantile(x,0.05),fun.max = function(x) quantile(x,0.95),fun =  function(x) mean(x))+
  #stat_summary(size=0.3) +
  theme(legend.position="bottom",text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1)) 

ggsave("./figures/astral_bl_u_support_corr.pdf",width=6, height=5)

################################
df <- read.table("./data/consistency_results.tsv")
colnames(df) <- c("tree", "rank", "count", "consistency")

ranks = df[df$tree=="16k.uDance" & df$count >= 15,]$rank
df$cat <- startsWith(df$rank, "d")

ggplot(aes(x=tree,y=rank, fill=consistency), data=df[df$rank %in% ranks,]) +  
  facet_grid(cat ~ ., scales = "free", space = "free_y") + geom_tile() + theme_classic()+
  theme(legend.position="top", 
        text = element_text(size=12), 
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  xlab("")+ylab("") +

ggsave("./figures/tree_comparison_consistency.pdf",width=8, height=7)


dfper <- merge(df[df$rank %in% ranks & df$tree %in% c("16k.uDance"),],
               df[df$rank %in% ranks & df$tree %in% c("10k.astral"),],
               by="rank")
dfper$diffperr = (dfper$count.x - dfper$count.y)/dfper$count.y
head(dfper)

ggplot(aes(y=reorder(rank, -diffperr), x=diffperr), data=dfper) + geom_bar(stat="identity") +
  theme_classic()+
  facet_grid(cat.y ~ ., scales = "free", space = "free_y") + xlab("Taxon sampling increase") + ylab("") +
  theme(strip.background = element_blank(),
                strip.text.y = element_blank() ) +
  scale_x_continuous(labels = scales::percent)

ggsave("./figures/taxon_sampling_increase_10_to_16.pdf",width=7, height=5)


df <- read.table("./data/consistency_results_fulltaxdump.tsv")
colnames(df) <- c("tree", "rank", "count", "consistency")

ranks = df[df$tree=="16k.uDance" & df$count >= 20,]$rank
df$cat <- startsWith(df$rank, "d")

df$tree <- as.factor(df$tree)
levels(df$tree) <- c("10k.concat", "10k.astral", "16k.uDance", "200k.uDance", "gtdb")
#acc_sca$selection=as.factor(acc_sca$selection)
#levels(acc_sca$selection) <- c("best","random")

ggplot(aes(x=tree,y=reorder(rank, count), fill=consistency), data=df[df$rank %in% ranks,]) +  
  facet_grid(reorder(cat,-count) ~ ., scales = "free", space = "free_y") + geom_tile() + theme_classic()+
  scale_fill_gradient2(low="#000000", mid="#FF0000",midpoint = .5 , high = "yellow")+
  theme(legend.position="top", 
        text = element_text(size=12), 
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  xlab("")+ylab("") 
  
ggsave("./figures/tree_comparison_consistency_fulltaxdump2.pdf",width=8, height=7)
  
############################

for( ctypeastr in c("rand", "cons", "astral")){
  df <-read.table(sprintf("./data/treedistcomp_sampled_%s.csv",ctypeastr), header=F)
  dim(df)
  df$type <- paste(df$V1, df$V2, sep = "-")
  ggplot(aes(color=type, x=V3, y=V4), data=df) + geom_point(alpha=0.1) + theme_classic() + 
    geom_abline(slope=1, intercept = 0, color="blue", linetype="dashed") +
    stat_smooth(se=F)+
    xlab("16k.uDance") + ylab(sprintf("10k.astral.%s", ctypeastr)) +
    scale_colour_manual(name = "", values=cbPalette[c(1,3,5,2,4,6)]) +
    theme(legend.position=c(0.2,0.8), 
          text = element_text(size=12)) + coord_cartesian(xlim = c(0,6.6), ylim=c(0,6.6))
    
  ggsave(sprintf("./figures/pairwise_distance_dotplot_withcpr_%s.png",ctypeastr),width=5, height=5)
  
  df2 <- read.table(sprintf("./data/treedistcomp_sampled_%s.csv",ctypeastr), header=F)
  df2[df2$V1 != "archaea",]$V1 = "bacteria"
  df2[df2$V2 != "archaea",]$V2 = "bacteria"
  df2$type <- paste(df2$V1, df2$V2, sep = "-")
  ggplot(aes(color=type, x=V3, y=V4), data=df2) + geom_point(alpha=0.1) + theme_classic() + 
    geom_abline(slope=1, intercept = 0, color="blue", linetype="dashed") +
    stat_smooth(se=F)+
    xlab("16k.uDance") + ylab(sprintf("10k.astral.%s", ctypeastr)) +
    scale_colour_manual(name = "", values=cbPalette[c(1,2,6)]) +
    theme(legend.position=c(0.2,0.8), 
          text = element_text(size=12)) + coord_cartesian(xlim = c(0,6.6), ylim=c(0,6.6))
  
  ggsave(sprintf("./figures/pairwise_distance_dotplot_%s.png", ctypeastr),width=5, height=5)
}


#######################################
df <-read.table("./data/ag-disag-16k-10krand.txt", header=F)
head(df)
colnames(df) <- c("tag", "rootd", "edgelen", "V4", "tree", "agreement")
df$rootdnorm = 0
for(i in c("tree1", "tree2")){
  m = max(df[df$tree==i,]$rootd)
  df[df$tree==i,]$rootdnorm <-  df[df$tree==i,]$rootd / m
}

ggplot(aes(color=tree, x=rootd, linetype=agreement),data=df) + geom_density() + theme_classic() +
  xlab("Depth") + ylab("Density") + 
  scale_colour_manual(name="",labels = c("16k.uDance", "10k.astral.rand"), values=cbPalette[c(1,2)]) +
  theme(legend.position=c(0.9,0.65)) 

ggsave("./figures/root2tip_16k_10k.pdf",width=9, height=3)

ggplot(aes(color=tree, x=rootdnorm, linetype=agreement),data=df) + geom_density() + theme_classic() +
  xlab("Depth") + ylab("Density") + 
  scale_colour_manual(name="",labels = c("16k.uDance", "10k.astral.rand"), values=cbPalette[c(1,2)]) +
  theme(legend.position=c(0.9,0.65)) 

ggsave("./figures/root2tip_16k_10k_norm.pdf",width=9, height=3)

ggplot(aes(color=tree, x=rootdnorm, linetype=agreement),data=df) + geom_density() + theme_classic() +
  xlab("Depth") + ylab("Density") + facet_wrap(tree~., ncol=1) + 
  scale_colour_manual(name="",labels = c("16k.uDance", "10k.astral.rand"), values=cbPalette[c(1,2)]) +
  theme(legend.position=c(0.9,0.8), 
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        strip.text.x = element_blank()) 

ggsave("./figures/root2tip_16k_10k_norm_facet.pdf",width=9, height=6)

ggplot(aes(x=tree, y=rootd, linetype=agreement),data=df) + geom_boxplot() + theme_classic() + 
  xlab("") + ylab("Depth") + 
  scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral.rand")) +
  theme(legend.position=c(0.8,0.8)) 

ggsave("./figures/root2tip_16k_10k_boxplot.pdf",width=4, height=4)

ggplot(aes(x=tree, y=rootdnorm, linetype=agreement),data=df) + geom_boxplot() + theme_classic() + 
  xlab("") + ylab("Depth") + 
  scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral.rand")) +
  theme(legend.position="right") 

ggsave("./figures/root2tip_16k_10k_boxplot_norm.pdf",width=5, height=4)


ggplot(aes(x=tree, y=edgelen, linetype=agreement, color=tree),data=df) + 
  geom_boxplot(outlier.alpha = 0.5) + theme_classic() + 
  xlab("") + ylab("Branch Length") + 
  scale_colour_manual(name="",labels = c("16k.uDance", "10k.astral.rand"), values=cbPalette[c(1,2)]) +
  scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral.rand")) +
  theme(legend.position="right") + coord_cartesian(ylim = c(0,0.9)) +
  scale_y_continuous(trans = 'sqrt', breaks = seq(0,0.9,0.1))  + guides(color="none") +
  theme(legend.position=c(0.85,0.85)) 

ggsave("./figures/bl_16k_10k_boxplot.pdf",width=5, height=5)

##########################

# df <- read.table("./data/aggonly-16k-10k.csv", header=T)
# head(df)
# df$dif = df$Len1 - df$Len2
# df$domain=""
# df[(1:551),]$domain="archaea"
# df[(552:dim(df)[1]),]$domain="bacteria"
# tail(df)
# ggplot(aes(x=reorder(Boot1,dif), y=dif),data = df) + geom_point(size=0.5) + 
#   facet_grid(domain ~ ., scales = "free", space = "free_x") +
#   geom_hline(yintercept = 0) + theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) 


##################
library(stringr)
df <-read.table("./data/gene_tree_discordance_16k.csv", header=F)
head(df)
colnames(df) <- c("quartet", "rf", "partition", "gene")
df$gene = as.factor(df$gene)
aggdf = aggregate(df, by=list(df$gene), FUN=mean)
x = aggdf[order(aggdf$quartet),]
#x=df[df$gene=="p0309",]

#df$gene <- factor(df$gene, levels=x$Group.1)

ggplot(aes(x=reorder(gene, quartet) , y=partition, fill=quartet), data=df) + geom_tile() +
  theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 

g2f <-read.table("data/gene_to_function.csv", sep=",")
head(g2f)
colnames(g2f) <- c("gene", "go")
mf <-read.delim("data/molecular_function_go.txt", sep="\t")
dim(mf)
colnames(mf) <- c("go", "goname")
mgd <- merge(g2f,mf, by="go", all.x = T)
dim(mgd)
aggdf$gene = aggdf$Group.1
mgddf <- merge(mgd,aggdf,by="gene", all.x=T)

mgddf$goname[ mgddf$goname %in% names(which(table(mgddf$goname) < 20)) ] = "other"

mgddf = mgddf[mgddf$go!="GO:0003674",]
ggplot(aes(y=1-quartet, x=reorder(str_wrap(goname, width=10),quartet,na.rm = TRUE) , color=goname), data=mgddf[mgddf$goname != "other",]) + 
  geom_boxplot() + theme_classic() + geom_jitter(width = 0.1)  +
  #scale_x_discrete(labels = wrap_format(10)) +
  ylab("Quartet distance") + xlab("Molecular function") +
  scale_color_brewer(palette = "Paired") +
  guides(colour = FALSE)

ggsave("./figures/gene_tree_discordance_16k.pdf",width=8, height=5)



###############################


dfall <- read.delim("data/stee_partition_all.csv",fill = T,na.strings = "", header=F)
head(dfall)

for (sz in c("500", "100")){
  for(i in c("QD", "nRF")){
    if (i == "QD") {
    colnames(dfall) <- c("rep", "cluster", "size", "uDance", "uDancer", 
                         "FastTree2", "FastTree2r", "concat", "concatr","mc")
    dfrf <- dfall[,c(1,2,3,4,6,8,10)]
    } 
    else {
      colnames(dfall) <- c("rep", "cluster", "size", "uDanceq", "uDance", 
                           "FastTreeq", "FastTree2", "concatq", "concat","mc")
      dfrf <- dfall[,c(1,2,3,5,7,9,10)]
    }
    
    dfm <- melt(data=dfrf, id=c("mc","rep","cluster", "size"))
    dfm$variable=as.factor(dfm$variable)
    levels(dfm$variable) = c( "uDance", "FT2-Astral", "concat")
    dfm$variable = factor(dfm$variable, levels=c( "uDance", "FT2-Astral", "concat"))
    #levels(dfm$variable) = c( "uDance", "FastTree2", "concat")
    dfm$mc = as.factor(dfm$mc)
    levels(dfm$mc) = c("mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500", "mc1", "mc2")
    dfm$mc = factor(dfm$mc,  levels=c("mc1", "mc2", "mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500"))
    
    
    x=dfm[dfm$size == sz & !is.na(dfm$value),]
    aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
    y=aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
    colnames(y) <-c("replicate", "mc", "method", "nvm")
    z=aggregate(y$replicate, by=list(y$mc, y$method), FUN=length)
    z$con = z$x
    colnames(z) <- c("mc", "variable", "x", "con")
    
    x=dfm[dfm$size == sz & is.na(dfm$value),]
    aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
    y=aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
    colnames(y) <-c("replicate", "mc", "method", "nvm")
    z2=aggregate(y$replicate, by=list(y$mc, y$method), FUN=length)
    z2$con = 10-z2$x
    colnames(z2) <- c("mc", "variable", "x", "con")
    z=rbind(z,z2)
    
    if(i == "QD"){
    ggplot(dfm[dfm$size == sz,], aes(x = variable, y = value, fill=variable)) + 
      facet_grid(.~mc) +
      geom_text(data=z,aes(y=0.0000001, x=variable, label=con))+
      geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + #facet_wrap(.~size)+
      stat_summary()+theme_classic() +
      scale_y_continuous(trans="log", breaks = c(0.0000001, 0.0000001, 0.000001,0.00001, 0.0001, 0.001, 0.01, 0.1, 1 ))+
      theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
      scale_fill_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
      #coord_cartesian(ylim=c(0.000001, 1)) + 
      guides(fill=F)
    }
    else{
      ggplot(dfm[dfm$size == sz,], aes(x = variable, y = value, fill=variable)) + 
        facet_grid(.~mc) +
        geom_text(data=z,aes(y=0.01, x=variable, label=con))+
        geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + #facet_wrap(.~size)+
        stat_summary()+theme_classic() +
        #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
        theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
        scale_fill_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
        #coord_cartesian(ylim=c(0.01, 0.8)) + 
        guides(fill=F)
    }
    
    ggsave(sprintf("./figures/violin_stee_partition_%s_%s.pdf", i,sz),width=8, height=4)
  }
}
i="nRF"
sz="100"

ggplot(dfm[dfm$size == sz & dfm$mc != "mc3-500",], aes(x = variable, y = value, fill=variable)) + 
  facet_grid(.~mc) +
  geom_text(data=z[z$mc!="mc3-500",],aes(y=0.01, x=variable, label=con))+
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + #facet_wrap(.~size)+
  stat_summary()+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
  #coord_cartesian(ylim=c(0.01, 0.8)) + 
  guides(fill=F)

ggsave(sprintf("./figures/violin_stee_partition_%s_%s_no500.pdf", i,sz),width=8, height=4)

###########################

dfall <- read.delim("data/stee_full_all.csv",fill = T,na.strings = "", header=F)
head(dfall)

for (sz in c("500", "100")){
  for(i in c("QD", "nRF")){
    if (i == "QD") {
      colnames(dfall) <- c("rep", "size", "uDance", "uDancer", 
                           "FastTree2", "FastTree2r", "concat", "concatr","mc")
      dfrf <- dfall[,c(1,2,3,5,7,9)]
    } 
    else {
      colnames(dfall) <- c("rep",  "size", "uDanceq", "uDance", 
                           "FastTreeq", "FastTree2", "concatq", "concat","mc")
      dfrf <- dfall[,c(1,2,4,6,8,9)]
    }

# colnames(dfall) <- c("rep", "size", "uDance", "uDancer", 
#                      "FastTree2", "FastTree2r", "concat", "concatr","mc")
# dfrf <- dfall[,c(1,2,3,5,7,9)]
  dfm <- melt(data=dfrf, id=c("mc","rep", "size"))
  dfm$variable=as.factor(dfm$variable)
  levels(dfm$variable) = c( "uDance", "FT2-Astral", "concat")
  dfm$variable = factor(dfm$variable, levels=c( "uDance", "FT2-Astral", "concat"))
  #levels(dfm$variable) = c( "uDance", "FastTree2", "concat")
  dfm$mc = as.factor(dfm$mc)
  levels(dfm$mc) = c("mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500", "mc1", "mc2")
  dfm$mc = factor(dfm$mc,  levels=c("mc1", "mc2", "mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500"))
  
  x=dfm[dfm$size == sz & !is.na(dfm$value),]
  aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
  y=aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
  colnames(y) <-c("replicate", "mc", "method", "nvm")
  z=aggregate(y$replicate, by=list(y$mc, y$method), FUN=length)
  z$con = z$x
  colnames(z) <- c("mc", "variable", "x", "con")
  
  x=dfm[dfm$size == sz & is.na(dfm$value),]
  aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
  y=aggregate(x$variable, by=list(x$rep,x$mc,x$variable), FUN=length)
  colnames(y) <-c("replicate", "mc", "method", "nvm")
  z2=aggregate(y$replicate, by=list(y$mc, y$method), FUN=length)
  z2$con = 10-z2$x
  colnames(z2) <- c("mc", "variable", "x", "con")
  z=rbind(z,z2)
  head(z)
  
  if(i == "QD"){
    ggplot(dfm[dfm$size == sz,], aes(x = variable, y = value)) + 
      facet_grid(.~mc) +
      geom_text(data=z,aes(y=0.000001, x=variable, label=con))+
      geom_point(aes(color=variable)) + #facet_wrap(.~size)+
      #geom_errorbar(aes(xmin = value - std.err, xmax = estimate + std.error), width = 0.3)+
      stat_summary( geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75)+
      theme_classic() +
      scale_y_continuous(trans="log", breaks = c(0.000001, 0.0000001, 0.000001,0.00001, 0.0001, 0.001, 0.01, 0.1, 1 ))+
      theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
      scale_color_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
      #coord_cartesian(ylim=c(0.000001, 1)) + 
      guides(color=F)
  }
  else {
    ggplot(dfm[dfm$size == sz,], aes(x = variable, y = value)) + 
      facet_grid(.~mc) +
      geom_text(data=z,aes(y=0.01, x=variable, label=con))+
      geom_point(aes(color=variable)) + #facet_wrap(.~size)+
      stat_summary( geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75)+ 
      theme_classic() +
      #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
      theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
      scale_color_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
      #coord_cartesian(ylim=c(0.01, 0.8)) + 
      guides(color=F)
  }
  
  
    ggsave(sprintf("./figures/violin_stee_full_%s_%s.pdf", i,sz),width=8, height=4)
  }
}

sz=100
i=nRF

ggplot(dfm[dfm$size == sz & dfm$mc != "mc3-500",], aes(x = variable, y = value)) + 
  facet_grid(.~mc) +
  geom_text(data=z[z$mc != "mc3-500",],aes(y=0.01, x=variable, label=con))+
  geom_point(aes(color=variable)) + #facet_wrap(.~size)+
  stat_summary( geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75)+ 
  theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
  #coord_cartesian(ylim=c(0.01, 0.8)) + 
  guides(color=F)

ggsave(sprintf("./figures/violin_stee_full_%s_%s_no500.pdf", i,sz),width=8, height=4)


#################################


dfall <- read.table("data/gene_discordance_all.csv")
head(dfall)
colnames(dfall) <- c("rep", "cluster", "gene", "size", "uDance", "FastTree2", "true", "uDance-GTEE", "FastTree2-GTEE","mc")
#dfall <- dfall[dfall$size == 100,]
df <- dfall[,c("mc", "rep", "cluster", "gene", "size", "uDance", "FastTree2", "true")]
dfm <- melt(data=df, id=c("mc", "rep","cluster", "gene", "size"))
dfm$variable=as.factor(dfm$variable)
dfm$variable = factor(dfm$variable, levels=c("true", "uDance", "FastTree2"))
levels(dfm$variable) = c("true", "uDance", "FastTree2")
dfm$mc = as.factor(dfm$mc)
levels(dfm$mc) = c("mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500", "mc1", "mc2")
dfm$mc = factor(dfm$mc,  levels=c("mc1", "mc2", "mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500"))


sz=500
i=nRF
ggplot(dfm[dfm$size == sz,], aes(x = variable, y = value, fill=variable)) + 
  facet_grid(.~mc) +
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + #facet_wrap(.~size)+
  stat_summary()+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = cbPalette[c(1,2,3)]) + xlab("")+ylab(sprintf("Gene Tree Discordance (%s)", i)) +
  #coord_cartesian(ylim=c(0.01, 0.8)) + 
  guides(fill=F)


ggsave("./figures/violin_gene_discordance_nRF_100.pdf",width=8, height=4)

dfd <- dfall[,c("mc", "rep", "cluster", "gene", "size", "uDance-GTEE", "FastTree2-GTEE")]
colnames(dfd) <- c("mc", "rep", "cluster", "gene", "size", "uDance", "FastTree2")


dfd <- melt(data=dfd, id=c("mc", "rep","cluster", "gene", "size"))

dfd$mc = as.factor(dfd$mc)
levels(dfd$mc) = c("mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500", "mc1", "mc2")
dfd$mc = factor(dfd$mc,  levels=c("mc1", "mc2", "mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500"))



ggplot(dfd[dfd$size == sz, ], aes(x = variable, y = value, fill=variable)) + 
  facet_grid(.~mc) +
  #geom_hline(yintercept = 0, color="red", linetype = 'dotted')+
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + #facet_wrap(.~size)+
  stat_summary()+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Gene Tree Estimation Error (%s)", i)) +
  #coord_cartesian(ylim=c(0.01, 0.8)) + 
  guides(fill=F)

# ggplot(dfd, aes(x = variable, y = value, fill=variable)) + 
#   geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) +  
#   geom_hline(yintercept = 0, color="red", linetype = 'dotted')+
#   facet_grid(.~mc)+
#   stat_summary()+theme_classic() +
#   scale_fill_manual(values = cbPalette[c(3,4)]) + xlab("Gene Tree Estimation Error")+ylab("nRF")

ggsave("./figures/violin_gene_delta_nRF_100.pdf",width=8, height=4)

######################################

df <-read.table("data/runtime_100_all.csv", sep=',')

colnames(df) <-c("Task", "cores", "time", "tottime", "size", "rep", "method")
goodreps <-df[df$method == "FT2+ASTRAL" & !is.na(df$time) & df$Task == "ASTRAL",]$rep

dfs <- df[df$rep %in% goodreps & df$size == 100,]

dfs$Task = as.factor(dfs$Task)
levels(dfs$Task) = c("APPLES-2", "ASTRAL", "ASTRAL", "ASTRAL", "FastTree-2", "IQTree", "IQTree", "RAxML-ng", "RAxML-ng" )
#fd$mc = factor(dfd$mc,  levels=c("mc1", "mc2", "mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500"))

dfs$method = as.factor(dfs$method)
dfs$method = factor(dfs$method, levels=c("uDance", "FT2+ASTRAL", "concat"))


ggplot(aes(x=reorder(as.factor(rep),tottime), y=tottime, fill=Task), data=dfs) + geom_bar(position="stack", stat="identity") +
  facet_wrap(.~method) + theme_classic() + scale_fill_brewer(palette="Set2") +
  ylab("Time (CPU seconds)") + xlab("Replicate ID") + #theme(axis.text.x = element_blank()) +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000))
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), 
                    breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000))

ggsave("./figures/hgt_runtime_100.pdf",width=5, height=4)

  
