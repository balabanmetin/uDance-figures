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

ggplot() + stat_ecdf(data=df3[df3$V4 != "max_clade",],aes(2*V2, color=V4, linetype=as.factor(V3))) + #geom_density(data=df, color="blue",aes(V1,linetype=V2)) +
  theme_classic() +scale_x_continuous(trans="log10", name="Novelty") + geom_vline(xintercept = 1) +#+ coord_cartesian(ylim=c(0,1.5*10^-7))
scale_linetype_manual(name="Backbone Size", values = c(1,3)) + scale_color_discrete(name="Selection Strategy") + 
  ylab("ECDF")

ggsave("./figures/taxon_subsampling_ecdf.pdf",width=6, height=5)

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
library(ggforce)
df <-read.table("./data/snakemake_timeline_500_mc-2-07.csv", header=F)
colnames(df) <- c("type", "rule", "time", "jobid", "cpus", "act_cpus", "method")
dfun <- df %>% group_by(time) %>% top_n(1, act_cpus)
#dfun <- distinct(df,time, .keep_all= TRUE)
#dfun <- df

ggplot(dfun, aes(x = time/60, y = act_cpus, color=method) )   + geom_path() + theme_classic() +
  scale_y_continuous(trans='log2', breaks=c(1,4,16, 64, 256, 1024), limits = c(1,1600)) +
  #facet_zoom(xlim = c(0, 850),zoom.size = 3)+
  annotate("rect", xmin = 0, xmax = 344, ymin = 0, ymax = 800, alpha = .1,fill = "blue")+
  annotate("text", x=172, y=1300, label="Backbone tree\ninference", size=2.5)+
  annotate("rect", xmin = 344, xmax = 880, ymin = 0, ymax = 800, alpha = .1,fill = "yellow")+
  annotate("text", x=622, y=1300, label="Backbone tree\nupdate", size=2.5)+  
  theme(legend.position = c(0.8,0.8),
        text = element_text(size=11),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.0))) + 
  labs( x = "Wall clock time (minutes)", y = "CPU usage") +
  scale_colour_manual(name = "", values=cbPalette[c(4,2,4,6)])

ggsave("./figures/udance_timeline_500.pdf",width=7.5, height=2.5)


ggplot(dfun, aes(x = time/60, y = act_cpus+1, group=1) ) + theme_classic()  + geom_path(color=cbPalette[4]) +
  #scale_y_continuous(trans='log2', breaks=c(1,4, 16,64, 256, 1024, 2048)) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1)) +
  labs( x = "Time (minutes)", y = "CPU usage")

ggsave("./figures/udance_timeline_500_lin.pdf",width=3, height=3)

dfc <-read.table("data/concat_timeline_data.csv", sep=",")
colnames(dfc) <- c("type", "rule", "time", "jobid", "cpus", "act_cpus", "method")
dfun <- rbind(dfun,dfc)

ggplot(dfun, aes(x = time/60, y = act_cpus, color=method) ) + theme_classic()  + geom_path() +
  scale_y_continuous(trans='log2', breaks=c(1,4, 16,64, 256, 1024, 2048)) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1)) +
  labs( x = "Time (minutes)", y = "CPU usage") 

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
  xlab("")+ylab("") 

ggsave("./figures/tree_comparison_consistency.pdf",width=8, height=7)



df <- read.table("./data/consistency_results_fulltaxdump_clean.tsv", sep = "\t")
colnames(df) <- c("tree", "rn", "rank", "count", "consistency")

ranks = df[df$tree=="16k.uDance" & df$count >= 20,]$rank
df$cat <- df$rn == "d"

df$tree <- as.factor(df$tree)
levels(df$tree) <- c("10k.concat", "10k.astral", "16k.uDance", "200k.uDance", "GTDB")
#acc_sca$selection=as.factor(acc_sca$selection)
#levels(acc_sca$selection) <- c("best","random")


dfper <- merge(df[df$rank %in% ranks & df$tree %in% c("16k.uDance"),],
               df[df$rank %in% ranks & df$tree %in% c("10k.astral"),],
               by="rank")
dfper$diffperr = (dfper$count.x - dfper$count.y)/dfper$count.y
dfperacc <- dfper
dfper <- merge(df[df$rank %in% ranks & df$tree %in% c("200k.uDance"),],
               df[df$rank %in% ranks & df$tree %in% c("10k.astral"),],
               by="rank")
dfper$diffperr = (dfper$count.x - dfper$count.y)/dfper$count.y
head(dfperacc)
dfperacc <- rbind(dfperacc,dfper)


ggplot(aes(fill=reorder(tree.x, diffperr), y=factor(rank, levels=levels(reorder(dfperacc$rank, dfperacc$diffperr))), x=diffperr),data=dfperacc) + 
  geom_bar(data=dfperacc[dfperacc$tree.x=="200k.uDance",], stat="identity", position=position_dodge(width = 0), color="black",size=0.2) +
  geom_bar(data=dfperacc[dfperacc$tree.x=="16k.uDance",], stat="identity", position=position_dodge(width = 0), color="black",size=0.2) +
  theme_classic()+
  facet_grid(reorder(cat.y, -count.x) ~ ., scales = "free", space = "free_y") + xlab("Taxon sampling increase") + ylab("") +
  theme(legend.position = c(0.65,0.2), strip.background = element_blank(),
        strip.text.y = element_blank() ) +
  scale_x_continuous() +
  scale_fill_grey(name="", labels=c("10k -> 16k", "10k -> 200k")) 

ggsave("./figures/taxon_sampling_increase_10_to_16and2000.pdf",width=150/40, height=200/40)

df$tree <- as.factor(df$tree)
levels(df$tree) <- c("10k.concat", "10k.astral", "16k.uDance", "200k.uDance", "GTDB")

ggplot(aes(x=tree,y=factor(rank, levels=levels(reorder(dfperacc$rank, dfperacc$diffperr))), fill=consistency), data=df[df$rank %in% ranks,]) +  
  facet_grid(reorder(cat,-count) ~ ., scales = "free", space = "free_y") + geom_tile(color="#BBBBBB") + theme_classic()+
  scale_fill_gradient2(name="",low="#000000", mid="#777777",midpoint = .5 , high = "#FFFFFF", breaks=c(0.25,0.5,0.75,1), labels=scales::percent)+
  theme(legend.position="top", 
        text = element_text(size=12), 
        axis.text.x = element_text(angle=25, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  xlab("")+ylab("") 

ggsave("./figures/tree_comparison_consistency_fulltaxdump2.pdf",width=110/25, height=225/25)
  
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


df<-read.table("data/ncbi_top30_rfqd_rename.csv", header=F, sep="\t")
colnames(df) <- c("qd", "rf", "method1", "method2", "phylum")


order_ref <-dfperacc[dfperacc$tree.x == "200k.uDance" & dfperacc$rank %in% df$phylum,]
ggplot(aes(fill=qd, x=method1, y=factor(phylum, levels=levels(reorder(levels(as.factor(df$phylum)), order_ref$diffperr)))), 
           data=df[df$method1 != "200k" & df$method2 == "200k",]) + geom_tile() + theme_classic() + #facet_grid(~method2) #+ 
  facet_grid(. ~ method2, scales = "free", space = "free_x") + geom_tile() + theme_classic()+
  theme(legend.position="top", 
        text = element_text(size=12), 
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  xlab("")+ylab("") +
scale_fill_gradient2(name= "QD   ", 
                     low = "#075AFF",
                     mid = "#FFFFCC",
                     high = "#FF0000") 


ggsave("./figures/qd_tree_compare_heat.pdf",width=4, height=6)

ggplot(aes(fill=qd, x=method1, y=factor(phylum, levels=levels(reorder(levels(as.factor(df$phylum)), order_ref$diffperr)))), 
       data=df) + geom_tile() + theme_classic() + #facet_grid(~method2) #+ 
  facet_grid(. ~ method2, scales = "free", space = "free_x") + geom_tile() + theme_classic()+
  theme(legend.position="top", 
        text = element_text(size=12), 
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  xlab("")+ylab("") +
  scale_fill_gradient2(name= "QD   ", 
                       low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 

ggsave("./figures/qd_tree_compare_heat_all.pdf",width=10, height=6)
#######################################
ctypeastr="astral"
df <-read.table(sprintf("./data/treedistcomp_sampled_%s.csv",ctypeastr), header=F)
dim(df)
df$type <- paste(df$V1, df$V2, sep = "-")
head(df)
colnames(df) <- c("o1", "o2", "16k.uDance", "10k.astral", "s1", "s2", "type")
dfm <- melt(data=df, measure=c("16k.uDance", "10k.astral"))
ggplot(aes(x=value, linetype=variable), data=dfm) + geom_density() + theme_classic() +
  facet_grid(type~.) + xlab("Root-to-tip distance") + ylab("Density") + theme(legend.title = element_blank())
ggsave(sprintf("./figures/pairwise_distance_density_%s.pdf", ctypeastr),width=8, height=7)


df2 <- read.table(sprintf("./data/treedistcomp_sampled_%s.csv",ctypeastr), header=F)
df2[df2$V1 != "archaea",]$V1 = "bacteria"
df2[df2$V2 != "archaea",]$V2 = "bacteria"
df2$type <- paste(df2$V1, df2$V2, sep = "-")
colnames(df2) <- c("o1", "o2", "16k.uDance", "10k.astral", "s1", "s2", "type")
dfm <- melt(data=df2, measure=c("16k.uDance", "10k.astral"))

ggplot(aes(x=value, linetype=variable), data=dfm) + geom_density() + theme_classic() +
  facet_grid(type~.) + xlab("Root-to-tip distance") + ylab("Density") + theme(legend.title = element_blank())

ggsave(sprintf("./figures/pairwise_distance_density_AB_%s.pdf", ctypeastr),width=8, height=4)


#######################################
#df <-read.table("./data/ag-disag-16k-10krand.txt", header=F)
df <-read.table("./data/ag-disag-16k-10kastralbl.txt", header=F)

head(df)
colnames(df) <- c("tag", "rootd", "edgelen", "V4", "tree", "agreement")
df$rootdnorm = 0
for(i in c("tree1", "tree2")){
  m = max(df[df$tree==i,]$rootd)
  df[df$tree==i,]$rootdnorm <-  df[df$tree==i,]$rootd / m
}

df2 <-read.table("./data/ag-disag-200k-16k.txt", header=F)

head(df2)
colnames(df2) <- c("tag", "rootd", "edgelen", "V4", "tree", "agreement")
df2$rootdnorm = 0
for(i in c("tree1", "tree2")){
  m = max(df2[df2$tree==i,]$rootd)
  df2[df2$tree==i,]$rootdnorm <-  df2[df2$tree==i,]$rootd / m
}

df$treebig="16k.uDance"
df2$treebig="200k.uDance"
df <- rbind(df2,df)

ggplot(aes(color=tree, x=rootd, linetype=agreement),data=df) + geom_density() + theme_classic() +
  xlab("Depth") + ylab("Density") + facet_wrap(.~treebig)+
  scale_colour_manual(name="",labels = c("Extended", "Backbone"), values=cbPalette[c(1,2)]) +
  theme(legend.position=c(0.9,0.65)) 

ggplot(aes(x=rootd, color=agreement),data=df[df$tree=="tree1",]) + stat_ecdf() + theme_classic() +
  xlab("Depth") + ylab("ECDF") + facet_wrap(.~treebig)+
  scale_y_continuous(trans="sqrt", breaks=c(0,0.01,0.04,0.09,0.16,0.25,0.36,0.49, 0.64,0.81, 1), 
                     label=percent) +
  scale_x_continuous(breaks = c(0,0.35,1,2,3,4))+
  #geom_vline(xintercept = 0.25, color="red", linetype="dashed ")+
  scale_colour_manual(name="", values=cbPalette[c(1,2)]) + coord_cartesian(xlim=c(0,3.6))+
  theme(legend.position=c(0.9,0.25),text = element_text(size=12)) 

ggsave("./figures/root2tip_200k_16k_10k.pdf",width=8, height=4)

ggplot(aes(x=rootd, color=agreement),data=df[df$tree=="tree1" & df$treebig == "16k.uDance",]) + stat_ecdf() + theme_classic() +
  xlab("Depth") + ylab("ECDF") + #facet_wrap(.~treebig)+
  scale_y_continuous(trans="sqrt", breaks=c(0,0.01,0.04,0.09,0.16,0.25,0.36,0.49, 0.64,0.81, 1), 
                     label=percent) +
  scale_x_continuous(breaks = c(0,0.35,1,2,3,4))+
  #geom_vline(xintercept = 0.25, color="red", linetype="dashed ")+
  scale_colour_manual(name="", values=cbPalette[c(1,2)]) + coord_cartesian(xlim=c(0,3.6))+
  theme(legend.position=c(0.8,0.25),text = element_text(size=10)) 

ggsave("./figures/root2tip_16k_10k_n.pdf",width=3.5, height=3.5)

ggplot(aes(x=rootdnorm, linetype=agreement),data=df[df$tree=="tree1",]) + geom_density() + theme_classic() +
  xlab("Depth") + ylab("Density") + facet_wrap(.~treebig)+
  scale_colour_manual(name="",labels = c("16k.uDance", "10k.astral"), values=cbPalette[c(1,2)]) +
  theme(legend.position=c(0.9,0.65)) 

ggsave("./figures/root2tip_200k_16k_10k_norm_density.pdf",width=9, height=3)

ggplot(aes(x=tree, y=rootd, linetype=agreement),data=df) + geom_boxplot() + theme_classic() + 
  xlab("") + ylab("Depth") + 
  scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral")) +
  theme(legend.position=c(0.8,0.8)) 

ggsave("./figures/root2tip_16k_10k_boxplot.pdf",width=4, height=4)

ggplot(aes(x=treebig, y=rootdnorm, linetype=agreement),data=df[df$tree=="tree1",]) + geom_boxplot() + theme_classic() + 
  xlab("") + ylab("Depth") + 
  #scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral")) +
  theme(legend.position="right") 

ggsave("./figures/root2tip_200k_16k_boxplot_norm.pdf",width=5, height=4)


ggplot(aes(x=edgelen, color=agreement),data=df[df$tree=="tree1",]) + 
  stat_ecdf() + theme_classic() + facet_wrap(.~treebig)+
  xlab("Branch Length") + ylab("ECDF") + 
  scale_colour_manual(name="", values=cbPalette[c(1,2)]) +
  #scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral")) +
  theme(legend.position="right") + coord_cartesian(xlim = c(0+0.000001,1)) +
  scale_x_continuous(trans = 'log10', breaks = c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1))  + #guides(color="none") +
  theme(legend.position=c(0.10,0.85)) 

ggsave("./figures/bl_200_16k_log_ecdf.pdf",width=8, height=4)

ggplot(aes(x=edgelen, color=agreement),data=df[df$tree=="tree1",]) + 
  stat_ecdf() + theme_classic() + facet_wrap(.~treebig)+
  xlab("Branch Length") + ylab("ECDF") + 
  scale_colour_manual(name="", values=cbPalette[c(1,2)]) +
  #scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral")) +
  theme(legend.position="right") + coord_cartesian(xlim = c(0+0.000001,1)) +
  scale_x_continuous(trans = 'sqrt', breaks=c(0,0.01,0.04,0.09,0.16,0.25,0.36,0.49, 0.64,0.81, 1))  + #guides(color="none") +
  theme(legend.position=c(0.9,0.2)) 

ggsave("./figures/bl_200_16k_sqrt_ecdf.pdf",width=8, height=4)

ggplot(aes(x=edgelen, color=agreement),data=df[df$tree=="tree1" & df$treebig == "16k.uDance",]) + 
  stat_ecdf() + theme_classic() + #facet_wrap(.~treebig)+
  xlab("Branch Length") + ylab("ECDF") + 
  scale_colour_manual(name="", values=cbPalette[c(1,2)]) +
  #scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral")) +
  theme(legend.position="right") + coord_cartesian(xlim = c(0+0.000001,1)) +
  scale_x_continuous(trans = 'sqrt', breaks=c(0,0.01,0.04,0.09,0.16,0.25,0.36,0.49, 0.64,0.81, 1))  + #guides(color="none") +
  theme(legend.position=c(0.8,0.2), text = element_text(size=10)) 

ggsave("./figures/bl_16k_sqrt_ecdf_n.pdf",width=3.5, height=3.5)


ggplot(aes(x=edgelen, linetype=tree, color=agreement),data=df) + 
  geom_density() + theme_classic() + 
  xlab("") + ylab("Branch Length") + 
  scale_colour_manual(name="",labels = c("16k.uDance", "10k.astral"), values=cbPalette[c(1,2)]) +
  #scale_x_discrete(name="",labels = c("16k.uDance", "10k.astral")) +
  theme(legend.position="right") + #coord_cartesian(ylim = c(0,0.9)) +
  scale_x_continuous(trans = 'sqrt', breaks = seq(0,0.9,0.1))  + #guides(color="none") +
  theme(legend.position=c(0.65,0.65)) 

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

######################
df <-read.table("./data/gene_tree_discordance_16k.csv", header=F)
head(df)
colnames(df) <- c("quartet", "rf", "partition", "gene")

ggplot(aes(x=reorder(gene, quartet) , y=as.factor(partition), fill=quartet), data=df[df$partition != 12, ]) + geom_tile() +
  theme_classic() + 
  theme(legend.position="top",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_gradient2(name = 'QD    ',low="#FFFFCC", mid="#FF0000",midpoint = .5 , high = "#000000")+
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5 )) +ylab("")

ggsave("./figures/qdheatmap_16k.pdf",width=12, height=4)
  

ggplot(aes(x=reorder(gene, rf) , y=reorder(partition,rf), fill=rf), data=df) + geom_tile() +
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_gradient2(low="yellow", mid="#FF0000",midpoint = .5 , high = "#000000")

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
mgddf2 <- mgddf

mgddf$goname[ mgddf$goname %in% names(which(table(mgddf$goname) < 20)) ] = "other"
mgddf2$go="ALL"
mgddf2$goname="All"
unique(mgddf2)
mgddf3=rbind(unique(mgddf2),mgddf)

mgddf3 = mgddf3[mgddf3$go!="GO:0003674",]
require(RColorBrewer)
xpal=brewer.pal(12, "Set3")
xpal[2]=xpal[12] # replace light yellow with normal yellow

ggplot(aes(y=1-quartet, x=reorder(str_wrap(goname, width=10),quartet,na.rm = TRUE) , color=goname), data=mgddf3[mgddf3$goname != "other",]) + 
  geom_boxplot() + theme_classic() + geom_jitter(width = 0.1)  +
  #scale_x_discrete(labels = wrap_format(10)) +
  ylab("Quartet similarity") + xlab("Molecular function") +
  scale_color_manual(values = xpal)+
  guides(colour = FALSE)

ggsave("./figures/gene_tree_discordance_16k_all.pdf",width=8, height=5)



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

colnames(dfall) <- c("rep", "size", "uDance-QD", "uDance-nRF", 
                     "FastTree2-QD", "FastTree2-nRF",
                     "ASTRID-QD", "ASTRID-nRF", "concat-QD", "concat-nRF","mc")
dfall$mc=paste(dfall$mc, dfall$size,sep="-")

dfd <- melt(data=dfall, id=c("mc", "rep", "size"))
head(dfd)

library(dplyr)
library(tidyr)

dfd = dfd %>% separate(variable, c("method", "measure"), "-") 

head(dfd)
colnames(dfd) <- c("mc", "rep", "size", "method","measure", "error")


dfd$mc = as.factor(dfd$mc)
levels(dfd$mc) = c("HD-P1", "HD-P2", "HD-P3", "HD-P4", "HD-P5", "LD-100", "MD-100", "MD-500", "HD-100", "HD-500")
dfd$mc = factor(dfd$mc, levels= c("LD-100", "MD-100", "MD-500", "HD-100", "HD-500","HD-P1", "HD-P2", "HD-P3", "HD-P4", "HD-P5"))

dfd$method = as.factor(dfd$method)
levels(dfd$method) = c("D", "C", "A","u" )
dfd$method = factor(dfd$method, levels=c("u","A","D","C"))

dfd = dfd[!(dfd$mc == "LD-100" & dfd$method=="A"),]

x=dfd
x$na=!is.na(x$error)
z = aggregate(x$na, by=list(x$mc,x$method, x$measure), FUN=sum)
colnames(z) <- c("mc", "method", "measure", "con")
z$con = 10-z$con
z[z$mc == "LD-100",]$con = z[z$mc == "LD-100",]$con-2



ggplot(dfd[dfd$measure == "nRF",], aes(x = method, y = error)) + 
  facet_grid(.~mc,scales = "free_x") +
  geom_text(data=z[z$measure == "nRF",],aes(y=0.000001, x=method, label=con, group=method),position = position_dodge(width = 0.7),size=3)+
  geom_point(alpha=0.6,aes(color=method), position = position_dodge(width = 0.7), size=0.75) + #facet_wrap(.~size)+
  #geom_errorbar(aes(xmin = value - std.err, xmax = estimate + std.error), width = 0.3)+
  stat_summary( aes(group=method), position = position_dodge(width = 0.7), geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .5)+
  stat_summary( aes(group=method), position = position_dodge(width = 0.7), geom = "errorbar", width = .35, size=0.5)+
  theme_classic() +
  coord_cartesian(ylim = c(0,0.25))+
  #scale_y_continuous(trans="log", breaks = c(0.000001, 0.0000001, 0.000001,0.00001, 0.0001, 0.001, 0.01, 0.1, 1 ))+
  #theme(legend.position = c(0.85,0.3), text = element_text(size=11), )+#axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(labels = c("uDance (u)", "FT2+ASTRAL (A)", "FT2+ASTRID (D)", "Concat (C)"), name="", values = cbPalette[c(2,3,7,4)]) + xlab("")+ylab("Species Tree Estimation Error (nRF)") +
  #coord_cartesian(ylim=c(0.000001, 1)) #+ 
  guides(color=guide_legend(nrow = 2),) +
  theme(legend.position = c(0.2,0.85),legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        text = element_text(size=11),
        panel.border = element_rect(fill="NA"),
        panel.spacing = unit(1.5,"pt"))

ggsave("./figures/violin_stee_full_all_nRF.pdf",width=9, height=3.5)

ggplot(dfd[dfd$measure == "QD",], aes(x = method, y = error)) + 
  facet_grid(.~mc,scales = "free_x") +
  geom_text(data=z[z$measure == "QD",],aes(y=0.000001, x=method, label=con, group=method),position = position_dodge(width = 0.7),size=3)+
  geom_point(alpha=0.6,aes(color=method), position = position_dodge(width = 0.7), size=0.75) + #facet_wrap(.~size)+
  #geom_errorbar(aes(xmin = value - std.err, xmax = estimate + std.error), width = 0.3)+
  stat_summary( aes(group=method), position = position_dodge(width = 0.7), geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .5)+
  stat_summary( aes(group=method), position = position_dodge(width = 0.7), geom = "errorbar", width = .35, size=0.5)+
  theme_classic() +
  #coord_cartesian(ylim = c(0,0.25))+
  scale_y_continuous(trans="log", breaks = c(0.000001, 0.0000001, 0.000001,0.00001, 0.0001, 0.001, 0.01, 0.1, 1 ))+
  #theme(legend.position = c(0.85,0.5)  )+#axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(labels = c("uDance (u)", "FT2+ASTRAL (A)", "FT2+ASTRID (D)", "Concat (C)"), name="", values = cbPalette[c(2,3,7,4)]) + xlab("")+ylab("Species Tree Estimation Error (QD)") +
  #coord_cartesian(ylim=c(0.000001, 1)) #+ 
  #guides(color=guide_legend(nrow = 2),) +
  guides(color=F,) +
  theme(legend.position = c(0.82,0.2),legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        text = element_text(size=11),
        panel.border = element_rect(fill="NA"),
        panel.spacing = unit(1.5,"pt"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  scale_x_discrete(position = "top") 


ggsave("./figures/violin_stee_full_all_QD.pdf",width=9, height=3)



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


sz=100
i="nRF"
ggplot(dfm[dfm$size == sz,], aes(x = mc, y = value, fill=variable)) + 
  #facet_grid(.~mc) +
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) , position = position_dodge(width = 0.75)) + #facet_wrap(.~size)+
  stat_summary(position = position_dodge(width = 0.75))+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(legend.position = "top", text = element_text(size=10))+#, axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(name="", values = cbPalette[c(1,2,3)]) + xlab("")+ylab(sprintf("Gene Tree Discordance (%s)", i)) #+
  #coord_cartesian(ylim=c(0.01, 0.8)) + 
  #guides(fill=F)


ggsave(sprintf("./figures/violin_gene_discordance_%s_%s.pdf", i,sz),width=8, height=4)


sz=500
i="nRF"
ggplot(dfm[dfm$size == sz,], aes(x = mc, y = value, fill=variable)) + 
  #facet_grid(.~mc) +
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), position = position_dodge(width = 0.8) ) + #facet_wrap(.~size)+
  stat_summary(position = position_dodge(width = 0.8))+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(legend.position = "top", text = element_text(size=10))+#, axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(name="", values = cbPalette[c(1,2,3)]) + xlab("")+ylab(sprintf("Gene Tree Discordance (%s)", i)) #+
#coord_cartesian(ylim=c(0.01, 0.8)) + 
#guides(fill=F)


ggsave(sprintf("./figures/violin_gene_discordance_%s_%s.pdf", i,sz),width=3, height=4)

##################################################
df <- read.table("./data/gtee_all_hgt.csv", header=F, sep="\t")
head(df)
colnames(df)=c("QD", "nRF", "size", "rep","cluster","gene","method","mc")
df$mc=paste(df$mc, df$size,sep="-")
#head(df[df$mc=="mc-5-100-500",])
     
dfd <- melt(data=df, id=c("mc", "rep","cluster", "gene", "size","method"))
head(dfd)
colnames(dfd) <- c("mc", "rep", "cluster", "gene", "size", "method","measure", "error")


dfd$mc = as.factor(dfd$mc)
levels(dfd$mc) = c("LD-100", "MD-100", "MD-500", "HD-100", "HD-500", "HD-P1", "HD-P2", "HD-P3", "HD-P4", "HD-P5")
dfd$method = as.factor(dfd$method)
levels(dfd$method) = c("A","u" )
dfd$method = factor(dfd$method, levels=c("u","A"))

dfd[dfd$measure=="QD" & dfd$error >= 0.67,]$error = 0.67
ggplot(dfd[dfd$measure == "nRF",], aes(x = method, y = error, fill=method)) + 
  facet_grid(.~mc) +
  #geom_hline(yintercept = 0, color="red", linetype = 'dotted')+
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_violin( trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), ) + #facet_wrap(.~size)+
  stat_summary(size=0.4, position = position_dodge(width = 0.9))+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(legend.position = c(0.05,0.85), 
        text = element_text(size=10),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=6))+#, axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(name="",labels=c("RAxML-NG", "FT2"), values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Gene Tree Estimation Error (%s)", "nRF"))  +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))#+


ggsave("./figures/violin_gene_delta_nRF.pdf",width=8, height=3)

ggplot(dfd[dfd$measure == "QD",], aes(x = method, y = error, fill=method)) + 
  facet_grid(.~mc) +
  #geom_hline(yintercept = 0, color="red", linetype = 'dotted')+
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_violin( trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), ) + #facet_wrap(.~size)+
  stat_summary(size=0.4, position = position_dodge(width = 0.9))+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(legend.position = c(0.05,0.85), 
        text = element_text(size=10),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=6))+#, axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(name="",labels=c("RAxML-NG", "FT2"), values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Gene Tree Estimation Error (%s)", "QD"))  +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))#+

ggsave("./figures/violin_gene_delta_QD.pdf",width=8, height=3)


################################################
df <-read.table("data/runtime_100_500_all.csv", sep=',')

colnames(df) <-c("Task", "cores", "time", "tottime", "size", "rep", "method")
goodreps <-df[df$method == "FT2+ASTRAL" & !is.na(df$time) & df$Task == "ASTRAL" & df$size == 100,]$rep
goodreps2 <-df[df$method == "concat" & !is.na(df$time) & df$Task == "FastTree-2" & df$size == 500,]$rep
goodreps <- goodreps[goodreps %in% goodreps2]

dfs <- df[df$rep %in% goodreps ,]

dfs$Task = as.factor(dfs$Task)
levels(dfs$Task) = c("APPLES-2", "ASTRAL", "ASTRAL", "ASTRAL", "ASTRID", "FastTree-2", "IQTree", "IQTree", "RAxML-ng", "RAxML-ng" )
#fd$mc = factor(dfd$mc,  levels=c("mc1", "mc2", "mc3-100", "mc3-200", "mc3-300", "mc3-400", "mc3-500"))

dfs$method = as.factor(dfs$method)
levels(dfs$method) = c("C", "A", "D", "u")
dfs$method = factor(dfs$method, levels=c("u", "A", "D", "C"))


dfm <- read.table("./data/hgt_memory.csv", header=F)
colnames(dfm) <- c("method", "rep", "size", "memory")

dfm$method = as.factor(dfm$method)
levels(dfm$method) = c("C", "A", "D", "u")
dfs$method = factor(dfs$method, levels=c("u", "A", "D", "C"))

ggplot(aes(x=method, y=tottime/3), data=dfs) + geom_bar(aes(fill=Task), position="stack", stat="identity") +
   theme_classic() + scale_fill_brewer(palette="Set2") + 
  geom_point(aes(y=memory/20000),data=dfm, alpha=0.5)+ 
  stat_summary(aes(y=memory/20000), data = dfm, geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75)+
  #stat_summary(aes(y=memory*5*10^4), data = dfm, geom = "errorbar", width = .25, size=0.3)+
  facet_wrap(.~size)+
  geom_text(data=dfs[dfs$method == "FT2+ASTRAL" & dfs$size == 500 ,],aes(y=7*10^5, x=method, label="Fail"), color="orange")+
  ylab("Cumulative time (CPU seconds)") + xlab("") + theme(legend.position = "left") +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000))
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), 
                     breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000),
                     sec.axis = sec_axis(~.*20000, name="Peak memory", labels = scales::label_bytes()))


ggsave("./figures/hgt_runtime_100_norep.pdf",width=4.5, height=3)


df<-read.table("/Users/metin/Workspace/btol/wol2_cpureport.txt", sep=",")

colnames(df) <-c("tree", "Task", "cores", "time", "tottime")
head(df)
dfs <- df
dfs$Task = as.factor(dfs$Task)
levels(dfs$Task) = c("APPLES-2", "ASTRAL", "IQTree", "RAxML-ng")

my_colors <- RColorBrewer::brewer.pal(5, "Set2")[c(1,2,4,5)]


ggplot(aes(x=tree, y=tottime/60/60), data=dfs) + geom_bar(aes(fill=Task), position="stack", stat="identity") +
  theme_classic() + scale_fill_manual(values = my_colors) + 
  #stat_summary(aes(y=memory*5*10^4), data = dfm, geom = "errorbar", width = .25, size=0.3)+
  ylab("Cumulative time (CPU hours)") + xlab("") + theme(legend.position = "left", axis.text.x = element_text(angle=25, hjust=1)) +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000))
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))#, 

ggsave("./figures/wol2_runtime.pdf",width=3.5, height=5)


############################################

df <- read.table("quartet_scores_root2tip_16k.csv", header=F)
colnames(df) <- c("qs", "r2t", "qs2", "bl"  )

#ggplot(data=df, aes(x=r2t, y=qs)) + geom_point(alpha=0.2) + geom_smooth(method = "loess", color="red") + theme_classic()

qnt= quantile(df$r2t,c(0:30)/30)
qnt = as.vector(qnt)
qnt[0] = 0
qnt_fixed = qnt 
ggplot(aes(x=cut(r2t,qnt_fixed, include.lowest = T), y = qs),data=df) + 
  theme_classic() + labs( y = "Quartet Score", x = "Depth") + 
  geom_hline(yintercept = 100/3, linetype="dotted")+
  #geom_vline(aes( xintercept=., color=method,linetype=datatype) , alpha=0.5,
  #           data = dcast(method+datatype+size~.,data=acc_sca[acc_sca$numcopy >= 20,],value.var = "error",fun.aggregate = mean))+
  #stat_summary(aes(y=qs),orientation = "x")+
  #geom_vline(aes(xintercept=mean(error),color=method,linetype=datatype), alpha=0.3 , size=0.5) +
  #geom_boxplot(aes(group=numcopy.x), outlier.alpha=0.5, outlier.size = 0.75) + 
  stat_smooth(se=F, color="red") + #facet_wrap(.~datatype, scales = "free") +
  stat_summary(alpha=0.8,  size=0.3, fun.min = function(x) quantile(x,0.05),fun.max = function(x) quantile(x,0.95),fun =  function(x) mean(x))+
theme(legend.position="bottom",text = element_text(size=12),axis.text.x = element_text(angle=90, hjust=1)) 


ggsave("./figures/depth_vs_qs_16k.pdf",width=5, height=5.5)


qnt= quantile(df$bl,c(0:30)/30)
qnt = as.vector(qnt)
qnt[1] = 0
qnt_fixed = qnt 
ggplot(aes(x=cut(bl,qnt_fixed, include.lowest = T), y = qs),data=df) + 
  theme_classic() + labs( y = "Quartet Score", x = "Branch length") + 
  geom_hline(yintercept = 100/3, linetype="dotted")+
  #geom_vline(aes( xintercept=., color=method,linetype=datatype) , alpha=0.5,
  #           data = dcast(method+datatype+size~.,data=acc_sca[acc_sca$numcopy >= 20,],value.var = "error",fun.aggregate = mean))+
  #stat_summary(aes(y=qs),orientation = "x")+
  #geom_vline(aes(xintercept=mean(error),color=method,linetype=datatype), alpha=0.3 , size=0.5) +
  #geom_boxplot(aes(group=numcopy.x), outlier.alpha=0.5, outlier.size = 0.75) + 
  stat_smooth(se=F, color="red") + #facet_wrap(.~datatype, scales = "free") +
  stat_summary(alpha=0.8,  size=0.3, fun.min = function(x) quantile(x,0.05),fun.max = function(x) quantile(x,0.95),fun =  function(x) mean(x))+
  theme(legend.position="bottom",text = element_text(size=12),axis.text.x = element_text(angle=90, hjust=1)) 

ggsave("./figures/bl_vs_qs_16k.pdf",width=5, height=5.5)

###################
df <-read.table("./data/vsize_run_mem.csv", header=F)
head(df)
colnames(df) <- c("CPUsec", "peakmem", "method", "size")

ggplot(data=df, aes(x=size, y=CPUsec, color=method)) + geom_point(aes(shape=method)) + 
  stat_smooth(aes(linetype="Runtime"),data=df, method="lm",se=F) + 
  scale_x_continuous(trans="log2",labels=comma) +  
  scale_y_continuous(trans="log2",sec.axis = sec_axis(~./32,name="Memory"))+
  #annotate(geom="text",label=format(lm(log(CPUsec)~log(size),df[ df$method == "concat",])[[1]][[2]],digits=2),color=cbPalette[2],x=35000,y=300000)+
  #annotate(geom="text",label=format(lm(log(CPUsec)~log(size),df[ df$method == "FT2+ASTRAL",])[[1]][[2]],digits=2),color=cbPalette[3],x=2200,y=1200000)+
  #annotate(geom="text",label=format(lm(log(CPUsec)~log(size),df[ df$method == "uDance",])[[1]][[2]],digits=2),color=cbPalette[4],x=85000,y=15000000)+
  xlab("Output tree size") + ylab("CPU time (seconds)") + theme_classic() +
  scale_colour_manual(name = "", values = cbPalette[c(2,3,4)])+
  scale_shape_manual(name="",values = c(19, 18, 15, 4)) + 
  stat_smooth( aes(size, peakmem*32, colour = method, shape=method, linetype="Memory"),
               data=df,method="lm",se=F) + 
  geom_point(aes(size, peakmem*32, color = method, shape=method),
             data=df,stat="summary",shape=5,fill="white",size=1.6) + 
  annotation_logticks(sides="lr") + 
  guides(color=guide_legend(nrow=1,byrow=T))+
  theme(legend.position = "top")+
  scale_linetype_manual(name="", values = c(1,3,1))
  
  
################################

df <-read.table("data/lgrates_all.csv", header=F)
colnames(df) <- c("catid", "rate", "partition", "gene", "tree")


df2 = df[df$tree == "16k" & df$catid == 1 & df$partition != 12,]
df2$partition = as.factor(df2$partition)
 
ggplot(aes(x=reorder(gene,-rate),y=partition, fill=rate-min(df2$rate)), data=df2) +  
   geom_tile(color="#BBBBBB") + theme_classic()+
  scale_fill_gradient2(name= "Rate   ", 
                      low = "#FFFFCC",
                      mid = "#FF0000",
                      midpoint=0.4,
                      high = "#000000",
                      breaks=c(0,0.2,0.4,0.6,0.8)-min(df2$rate),
                      labels=function(x) x+min(df2$rate))+   
  #scale_fill_gradient2(name="",low="#FFFFFF", 
  #                     mid="#777777",midpoint = .5 , 
  #                     high = "#000000")+
  theme(legend.position="top", 
        text = element_text(size=12), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5 )) +
  xlab("")+ylab("") 



ggsave("./figures/lgrates1_new.pdf",width=12, height=4)


df2 = df[df$tree == "16k" & df$catid == 4 & df$partition != 12,]
df2$partition = as.factor(df2$partition)

 
#ggplot(aes(x=reorder(gene,-rate),y=partition, fill=rate), data=df2) +  
#  geom_tile(color="#BBBBBB") + theme_classic()+
#  scale_fill_gradient(name= "Rate   ", 
#                       low = "#FFFFCC",
#                       high = "#FF0000") +
#  #scale_fill_gradient2(name="",low="red", 
#  #                     mid="yellow",midpoint = .5 , 
#  #                     high = "#007f00")+
#  theme(legend.position="top", 
#        text = element_text(size=12), 
#        axis.text.x = element_blank(),
#        strip.background = element_blank(),
#        strip.text.y = element_blank()) + 
#  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
#  xlab("")+ylab("") 


ggplot(aes(x=reorder(gene,rate),y=partition, fill=rate), data=df2) +  
  geom_tile(color="#BBBBBB") + theme_classic()+
  scale_fill_gradient2(name= "Rate   ", 
                       high = "#FFFFCC",
                       mid = "#FF0000",
                       midpoint=1.7,
                       low = "#000000",
                       breaks=c(1.2,1.6,2.0,2.4,2.8, 3.2),
                       labels=function(x) x)+   
  #scale_fill_gradient2(name="",low="#FFFFFF", 
  #                     mid="#777777",midpoint = .5 , 
  #                     high = "#000000")+
  theme(legend.position="top", 
        text = element_text(size=12), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5 )) +
  xlab("")+ylab("") 


ggsave("./figures/lgrates4_new.pdf",width=12, height=4)


##############33 

df <- read.table("data/lpps_noblrees.csv", header=F)
colnames(df) <- c("lpp", "Tree")

ggplot(data=df, aes(color=Tree,x=lpp)) + stat_ecdf() + theme_classic() +
  xlab("Local posterior probability") + ylab("ECDF")

df <-read.table('data/lpp_200k.csv', header=F)
colnames(df) <- c("lpp", "Tree", "cluster")

ggplot(data=df, aes(color=as.factor(cluster%%10),x=lpp)) + stat_ecdf() + theme_classic() +
  facet_wrap(.~cluster%/%10)+scale_color_brewer(palette="Spectral")+
  xlab("Local posterior probability") + ylab("ECDF")

summr <- dcast(data=df, cluster~.)

ggsave("./figures/lpp_undance_compare.pdf",width=5, height=4.5)


df2 <-read.table("data/200k_diversity.csv", header=F)
head(df2)
colnames(df2) <- c("diversity", "cluster")
mgd <- merge(merge(df,df2,by="cluster"),dcast(data=df, cluster~'count'),by="cluster")
head(mgd)

ggplot(data=mgd, aes(linetype=cut(diversity,5),color=as.factor(cluster%/%10),x=lpp)) + stat_ecdf() + theme_classic() +
  facet_wrap(.~cluster%%10)+#scale_color_brewer(palette="Set2")+
  xlab("Local posterior probability") + ylab("ECDF")


ggsave("./figures/lpp_each_cluster_ecdf.pdf",width=10, height=6)

ggplot(data=mgd, aes(x=diversity,y=lpp,color=cut(count,breaks=c(1200,1800,2500,6000)))) + 
  stat_summary(alpha=0.7, fun.data = function(x) data.frame(ymin=quantile(x,0.25), ymax=quantile(x,0.75), y=median(x))) + 
  theme_classic() + 
  scale_x_continuous(trans = "log10")+
  stat_smooth(se=F, color="red") +
  scale_color_brewer(name="Size", palette="Dark2")+
  xlab("Average branch length") + ylab("Local posterior probability") #+ geom_text(position=position_jitter(width = 0),aes(label=cluster,x=diversity,y=0.4), color="black",data=df2)

ggplot(data=mgd, aes(x=lpp,group = cluster, color=as.factor(cluster))) + 
  stat_ecdf()+
  scale_color_discrete(name="Partition")+
  theme_classic() +
  ylab("ECDF") + xlab("Local posterior probability")


ggsave("./figures/lpp_16k_ecdf.pdf",width=6, height=5.5)

ggplot(data=mgd, aes(x=as.factor(cluster),y=lpp, color=count)) + 
  geom_boxplot()+
  scale_color_binned(name="Partition")+
  theme_classic() 

ggsave("./figures/lpp_average_bl.pdf",width=8, height=4)


#############TEMP-ACCURACY###################
df <-read.table("./data/vsize_notime.csv", header=F)
colnames(df) <- c("qd","rf","method","size")
head(df)

ggplot(aes(color=method, x=size/1000, y=rf), data=df) + geom_path() + theme_classic() +
  xlab("Size (x 1000)") + ylab("nRF") + scale_x_continuous(trans="log2", breaks = c(0.5,1,2,4,8,16,32,64))+
  scale_color_manual(name="", values = cbPalette[c(4,3,2)])+ coord_cartesian(ylim = c(0,0.22))+
  theme(legend.position = c(0.8,0.8))

ggsave("./figures/vsize_nRF.pdf",width=3.5, height=3.5)

ggplot(aes(color=method, x=size/1000, y=qd), data=df) + geom_path() + theme_classic() +
  xlab("Size (x 1000)") + ylab("QD")+ scale_x_continuous(trans="log2", breaks = c(0.5,1,2,4,8,16,32,64))+
  scale_color_manual(name="", values = cbPalette[c(4,3,2)])+ coord_cartesian(ylim = c(0,0.025))+
  theme(legend.position = "None")

ggsave("./figures/vsize_qd.pdf",width=3.5, height=3.5)

##############RUNTIME######################

df <-read.table("./data/vsize_incremental_combined.csv", header=F)
colnames(df) <- c("qd","rf","method","size","time","memory")
head(df)

ggplot(aes(color=method, x=size/1000, y=time/60/60), data=df) + geom_path() + theme_classic() +
  xlab("Size (x 1000)") + ylab("Time (CPU hours)")+ scale_x_continuous(trans="log2", breaks = c(0.5,1,2,4,8,16,32,64))+
  theme(legend.position = "None") + scale_y_continuous(trans="log2") +
  scale_color_manual(name="", values = cbPalette[c(4,3,2)])+
  annotate(geom="text",label=format(lm(log(time/60/60)~log(size),df[df$method == "concat",])[[1]][[2]],digits=3),color=cbPalette[4],x=32000/1000,y=30) +
  annotate(geom="text",label=format(lm(log(time/60/60)~log(size),df[df$method == "FT2+ASTRAL",])[[1]][[2]],digits=3),color=cbPalette[3],x=2000/1000,y=320)+
  annotate(geom="text",label=format(lm(log(time/60/60)~log(size),df[df$method == "uDance",])[[1]][[2]],digits=3),color=cbPalette[2],x=64000/1000,y=2600)

  

ggsave("./figures/vsize_runtime.pdf",width=3.5, height=3.5)



ggplot(aes(color=method, x=size/1000, y=memory/1000), data=df) + geom_path() + theme_classic() +
  xlab("Size (x 1000)") + ylab("Memory (GB)")+ scale_x_continuous(trans="log2", breaks = c(0.5,1,2,4,8,16,32,64))+
  theme(legend.position = "None") + scale_y_continuous() +
  scale_color_manual(name="", values = cbPalette[c(4,3,2)])

ggsave("./figures/vsize_memory.pdf",width=3.5, height=3.5)

dfc <- read.table("./data/clocktimes.csv", header=F)
colnames(dfc) <- c("wtime","size","method")

ggplot(aes(color=method, x=size/1000, y=wtime/60/60), data=dfc) + geom_path() + theme_classic() +
  xlab("Size (x 1000)") + ylab("Time (Wall-clock hours)")+ scale_x_continuous(trans="log2", breaks = c(0.5,1,2,4,8,16,32,64))+
  theme(legend.position = "None") + scale_y_continuous(trans="log2") +
  scale_color_manual(name="", values = cbPalette[c(4,3,2)]) +
  annotate(geom="text",label=format(lm(log(wtime/60/60)~log(size),dfc[dfc$method == "concat",])[[1]][[2]],digits=3),color=cbPalette[4],x=32000/1000,y=26) +
  annotate(geom="text",label=format(lm(log(wtime/60/60)~log(size),dfc[dfc$method == "FT2+ASTRAL",])[[1]][[2]],digits=3),color=cbPalette[3],x=2000/1000,y=20)+
  annotate(geom="text",label=format(lm(log(wtime/60/60)~log(size),dfc[dfc$method == "uDance",])[[1]][[2]],digits=3),color=cbPalette[2],x=64000/1000,y=19)

 
ggsave("./figures/vsize_wallclock.pdf",width=3.5, height=3.5) 


ggplot(data=df, aes(x=size/1000, y=time/60/60, color=method)) + geom_point(aes(shape=method)) + 
  stat_smooth(aes(linetype="CPUtime"),data=df, method="lm",se=F) + 
  scale_x_continuous(trans="log2", breaks = c(0.5,1,2,4,8,16,32,64)) +  
  annotate(geom="text",label=format(lm(log(wtime/60/60)~log(size),dfc[dfc$method == "concat",])[[1]][[2]],digits=3),color=cbPalette[4],x=64000* 1.3/1000,y=50) +
  annotate(geom="text",label=format(lm(log(wtime/60/60)~log(size),dfc[dfc$method == "FT2+ASTRAL",])[[1]][[2]],digits=3),color=cbPalette[3],x=2000*1.4/1000,y=20)+
  annotate(geom="text",label=format(lm(log(wtime/60/60)~log(size),dfc[dfc$method == "uDance",])[[1]][[2]],digits=2),color=cbPalette[2],x=64000*1.3/1000,y=19)+ 
  annotate(geom="text",label=format(lm(log(time/60/60)~log(size),df[df$method == "concat",])[[1]][[2]],digits=3),color=cbPalette[4],x=64000*1.3/1000,y=160) +
  annotate(geom="text",label=format(lm(log(time/60/60)~log(size),df[df$method == "FT2+ASTRAL",])[[1]][[2]],digits=3),color=cbPalette[3],x=2000*1.4/1000,y=320)+
  annotate(geom="text",label=format(lm(log(time/60/60)~log(size),df[df$method == "uDance",])[[1]][[2]],digits=3),color=cbPalette[2],x=64000*1.3/1000,y=4096)+
  xlab("Size (x 1000)") + ylab("Time (CPU hours)") + theme_classic() +
  scale_colour_manual(name = "", values = cbPalette[c(4,3,2)])+
  scale_y_continuous(trans="log2",sec.axis = sec_axis(~./1,name="Time (Wall-clock hours)")) +
  stat_smooth( aes(size/1000, 1*wtime/60/60, colour = method, linetype="Clocktime"),
               data=dfc,method="lm",se=F) + 
  geom_point(aes(size/1000, 1*wtime/60/60, color = method, shape=method),
             data=dfc,stat="summary",fill="white",size=1.6) + 
  scale_shape_manual(name="",values = c(19, 18, 15, 4)) + 
  #annotation_logticks(sides="lr") + 
  guides(color=guide_legend(nrow=1,byrow=T))+
  theme(legend.position = "top")+
  scale_linetype_manual(name="", values = c(1,3,1)) 

ggsave("./figures/vsize_cpu_wallclock.pdf",width=5, height=4) 


df <-read.table("./data/vsize_incremental_rep2.csv", header=F)
colnames(df) <- c("qd","rf","size","method")
head(df)

ggplot(aes(color=method, x=size, y=rf), data=df) + geom_path() + theme_classic() +
  xlab("Size") + ylab("nRF") + scale_x_continuous(trans="log2", breaks = c(500,1000,2000,4000,8000))+
  scale_color_manual(name="", values = cbPalette[c(2,4)])+
  theme(legend.position = c(0.8,0.8))

ggsave("./figures/vsize_nRF_rep2.pdf",width=5, height=5)


ggplot(aes(color=method, x=size, y=qd), data=df) + geom_path() + theme_classic() +
  xlab("Size") + ylab("QD") + scale_x_continuous(trans="log2", breaks = c(500,1000,2000,4000,8000,16000,32000))+
  scale_color_manual(name="", values = cbPalette[c(2,4)])+
  theme(legend.position = c(0.8,0.2))


ggsave("./figures/vsize_QD_rep2.pdf",width=5, height=5)



df <-read.table("data/vsize_10rep.csv", header=F)
head(df)
colnames(df) <- c("qd","rf","size","method","rep")
ggplot(aes(color=method, x=as.factor(size), y=qd), data=df) + geom_boxplot() + theme_classic() +
  xlab("Size") + ylab("QD")+ #+ scale_x_continuous(trans="log2", breaks = c(500,1000,2000,4000,8000,16000,32000))+
  scale_color_manual(name="", values = cbPalette[c(4,2)])+
  theme(legend.position = "right")

ggsave("./figures/vsize_QD_10rep.pdf",width=6, height=5)

ggplot(aes(color=method, x=as.factor(size), y=rf), data=df) + geom_boxplot() + theme_classic() +
  xlab("Size") + ylab("nRF")+ #+ scale_x_continuous(trans="log2", breaks = c(500,1000,2000,4000,8000,16000,32000))+
  scale_color_manual(name="", values = cbPalette[c(4,2)])+
  theme(legend.position = "right")



ggsave("./figures/vsize_nRF_10rep.pdf",width=6, height=5)
###############################

df <- read.table("data/ppscomp.csv", header=F)


colnames(df) <- c("pp","method")


ggplot(aes(color=method, x=pp), data=df) + stat_ecdf() + theme_classic()


###########################


df <-read.table("data/smallinsert_runtime.csv", header=F)
head(df) 

colnames(df) = c("time", "memory", "method", "qsize", "mode")

df2 <-read.table("data/smallinsert_acc.csv", header=F)
head(df2)

colnames(df2) = c("delta", "error", "bberror", "dk","mode","qsize", "method")
df2$method = as.factor(df2$method)
levels(df2$method) = c("APPLES2","udance.incremental","udance.maxqs")

df3 <- merge(df, df2, by=c("mode", "qsize", "method"))

head(df3)

df3$mode = as.factor(df3$mode)
levels(df3$mode) = c("auto", "fast")
df3$fmethod = paste(df3$method, df3$mode,sep="-")
head(df3)
df3 = df3[df3$fmethod != "APPLES2-fast",]

ggplot(data=df3, aes(y=delta,x=time/60/60, color=fmethod))+ 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y="mean", geom="line", aes(group=fmethod))+
  scale_x_continuous(trans="log", breaks=c(10, 40, 160, 640,2560))+
  xlab("CPU hours")+
  #facet_wrap(nrow=1, .~qsize) + 
  theme_classic()

ggplot(data=df3, aes(x=as.factor(qsize),y=delta, color=fmethod))+ 
  #geom_boxplot()+
  stat_summary(aes(size=time/60/60),fun.y = mean, alpha=0.5, position = position_dodge(width = 0.7)) +
  #stat_summary(fun.y="mean", geom="line", aes(group=fmethod))+
  scale_size_area(max_size = 3)+
  #scale_x_continuous(trans="log", breaks=c(10, 40, 160, 640,2560))+
  xlab("CPU hours")+
  #facet_wrap(nrow=1, .~qsize) + 
  theme_classic()


ggplot(data=merge(
  dcast(df3,qsize+fmethod~"delta",value.var = "delta",fun.aggregate = mean),
  dcast(df3,qsize+fmethod~"time",value.var = "time",fun.aggregate = mean)), 
  aes(x=as.factor(qsize),y=time/60/60, shape=fmethod,
      fill=cut(delta,breaks=c(4,2,1,0,0.5,-1,-4,-8,-16,-32,-64)),
      color=cut(delta,breaks=c(4,2,1,0,0.5,-1,-4,-8,-16,-32,-64))))+ 
  #geom_boxplot()+
  stat_summary(aes(group=fmethod),color="gray30",size=0.1,
               position = position_dodge(width = 0.7),geom="line") +
  stat_summary(position = position_dodge(width = 0.7),geom="point",size=2.75,color="black") +
  stat_summary(position = position_dodge(width = 0.7),geom="point",size=2) +
    #stat_summary(fun.y="mean", geom="line", aes(group=fmethod))+
  #scale_size_area(max_size = 3)+
  #scale_fill_binned(breaks=c(3,0,-1,-2,-4,-10,-50))+
  #scale_fill_gradient2(midpoint  = -10,high="#FF0000",mid="#4455EE", low="#00EE55")+
  #scale_x_continuous(trans="log", breaks=c(10, 40, 160, 640,2560))+
  scale_fill_brewer(name="Mean delta error", palette = "Spectral",direction=-1)+
  scale_color_brewer(name="Mean delta error", palette = "Spectral",direction=-1)+
  xlab("Number of queries")+   
  scale_y_continuous(name="CPU time (hours)",trans="identity")+
  scale_shape_manual(name="", values=c(21,23,24,22,25))+
  #facet_wrap(nrow=1, .~qsize) + 
  theme_classic()+
  theme(legend.position = c(0.22,.65))


ggsave("./figures/smallinsert_paths.pdf",width=7, height=7)

df4 <- read.table("./data/smallinsert_acc_qdrf.csv", header=F)
colnames(df4) = c("qd", "rf", "qsize", "mode","method")
head(df4)

df5 <- merge(df, df4, by=c("mode", "qsize", "method"))
df5$mode = as.factor(df5$mode)
levels(df5$mode) = c("auto", "fast")
df5$fmethod = paste(df5$method, df5$mode,sep="-")
df5 = df5[df5$fmethod != "APPLES2-fast",]

ggplot(data=df5, aes(y=rf,x=time/60/60, color=fmethod))+ #stat_summary(fun.y = mean, geom = "point") +
  #stat_summary(fun.y="mean", geom="line", aes(group=fmethod))+
  geom_point()+ geom_line()+
  scale_x_continuous(trans="log", breaks=c(10, 40, 160, 640,2560))+
  xlab("CPU hours")+
  #facet_wrap(nrow=1, .~qsize) + 
  theme_classic()

  #acc_sca$selection=as.factor(acc_sca$selection)
#levels(acc_sca$selection) <- c("best","random")

aggdata <- aggregate(df3, by=list(df3$fmethod, df3$qsize), FUN= function(c)sum(c==0)/sum(c>=0))
aggdata$Algorithm  = aggdata$Group.1
aggdata$size  = aggdata$Group.2
head(aggdata)

aggdata2 <- aggregate(df3, by=list(df3$fmethod, df3$qsize), FUN=mean)
aggdata2$Algorithm  = aggdata2$Group.1
aggdata2$size  = aggdata2$Group.2
head(aggdata2)

ggplot(data=df, aes(y=qsize, y=time/60/60/24)) + geom_path() +
  scale_x_continuous(trans="log", breaks=c(5,10,20,40,80,160,320,640,1280)) +
  ylab("CPU hours")+ xlab("Number of queries") + theme_classic()

ggsave("./figures/smallinsert.pdf",width=5, height=5)



##########################

df2 <-read.table("data/all_delta_errs.csv", header=F)
head(df2)
colnames(df2) = c("error", "rfpl", "rfbb","query","qsize","method")
ggplot(data=df2, aes(color=method, x=as.factor(qsize),y=rfpl)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
  theme_classic() #+ coord_cartesian(ylim= c(0,2))


########################

df <-read.table("data/16internalexternal.csv",header=F,sep=',')
#ggplot(data=df,aes(x=V1,fill=V2, y=V2)) + geom_violin(trim = FALSE,draw_quantiles = c(0.25, 0.5, 0.75)) + theme_classic() + xlab("BL")
df[df$V3=="sim",]$V1 = df[df$V3=="sim",]$V1/mean(df[df$V3=="sim",]$V1)
df[df$V3=="real",]$V1 = df[df$V3=="real",]$V1/mean(df[df$V3=="real",]$V1)
df[df$V3=="sim-inf",]$V1 = df[df$V3=="sim-inf",]$V1/mean(df[df$V3=="sim-inf",]$V1)
ggplot(data=df,aes(x=V1,color=V3, linetype=V2)) + stat_ecdf() + theme_classic()

ggsave("./figures/internal-external-ratio.pdf",width=5.5, height=5)


ggplot(data=df,aes(x=V1,fill=V2, y=V2)) +  facet_grid(V3~.) + 
  geom_violin(trim = FALSE)+#,draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_boxplot(width=.1,outlier.size = 0) +
  theme_bw() + xlab("BL")

ggsave("./figures/internal-external-ratio-violin.pdf",width=5.5, height=5)

############################
df <- read.table("./data/bratiosl_level.csv", header=F)
head(df)

colnames(df) <- c("level", "ratio", "denom", "datatype")
ggplot(df, aes(x = level, y = log10(ratio))) + 
  facet_grid(.~denom) +
  #geom_text(data=z[z$mc != "mc3-500",],aes(y=0.02, x=mc, label=con, color=variable),position = position_dodge(width = 0.7))+
  geom_point(alpha=0.5,aes(color=datatype), position = position_dodge(width = 0.7)) + #facet_wrap(.~size)+
  geom_hline(yintercept = 0, linetype="dotted") + 
  #geom_errorbar(aes(xmin = value - std.err, xmax = estimate + std.error), width = 0.3)+
  stat_summary( aes(group=datatype), position = position_dodge(width = 0.7), geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .5)+
  stat_summary( aes(group=datatype), position = position_dodge(width = 0.7), geom = "errorbar", width = .25, size=0.3)+
  coord_cartesian(ylim = c(-1,1))+
  theme_classic() +
  #scale_y_continuous(trans="sqrt", breaks = c(0.04,  0.06,  0.08,0.1 ))+
  theme(legend.position = "right", text = element_text(size=11), )+#axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(name="", values = cbPalette[c(2,3,4)]) + xlab("Branch Level")+ylab("Log Normalized BL") #+

ggsave("./figures/bratios_levels.pdf",width=6, height=4)

ggplot(df[df$denom == "all" & df$level <= 4,], aes( x = log10(ratio), color=datatype)) + 
  #facet_grid(.~denom) +
  #geom_text(data=z[z$mc != "mc3-500",],aes(y=0.02, x=mc, label=con, color=variable),position = position_dodge(width = 0.7))+
  #geom_point(alpha=0.5,aes(color=datatype), position = position_dodge(width = 0.7)) + #facet_wrap(.~size)+
  geom_density()+
  geom_vline(xintercept = 0, linetype="dotted") + 
  #geom_errorbar(aes(xmin = value - std.err, xmax = estimate + std.error), width = 0.3)+
  #stat_summary( aes(group=datatype), position = position_dodge(width = 0.7), geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .5)+
  #stat_summary( aes(group=datatype), position = position_dodge(width = 0.7), geom = "errorbar", width = .25, size=0.3)+
  coord_cartesian(xlim = c(-1,1))+
  theme_classic() +
  #scale_y_continuous(trans="sqrt", breaks = c(0.04,  0.06,  0.08,0.1 ))+
  theme(legend.position = c(0.15,0.85), text = element_text(size=11), )+#axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(name="", values = cbPalette[c(2,3,4)]) + ylab("Density")+xlab("Log Normalized BL") #+

ggsave("./figures/bratios_levels_new.pdf",width=4, height=4)

##############################

df <- read.table("./data/subs_collapse.csv", header=F)
head(df)
colnames(df) <- c("numspecies", "numinternal", "collapse",  "tree")

df$ratio = df$numinternal/(df$numspecies-1)
ggplot(df,aes(color=tree, y=ratio, x=collapse)) + geom_line(aes(group=tree)) + theme_classic() +
  xlab("Branch support collapse threshold") +
  ylab("Proportion of remaining branches") +
  scale_color_manual(name="", values = cbPalette[c(2,3,4)]) +
  theme(legend.position = c(0.15,0.15), text = element_text(size=11), )

ggsave("./figures/support_collapse.pdf",width=4, height=4) 

ggplot(df,aes(color=tree, y=numinternal, x=collapse)) + geom_line(aes(group=tree)) + theme_classic() +
  xlab("Branch support collapse threshold") +
  ylab("Number of remaining branches") +
  scale_color_manual(name="", values = cbPalette[c(2,3,4)]) +
  theme(legend.position = c(0.85,0.85), text = element_text(size=11), )

ggsave("./figures/support_collapse_absolute.pdf",width=4, height=4) 


###################33

df <- read.table("./data/discacross.csv", header=F)
head(df)
colnames(df) <- c("QD","nRF","gene","rep","dataset")

#df$dataset = as.factor(df$dataset)
#levels(df$dataset) = c("10k", "MD-500","HD-500")

ggplot(data=df, aes(color=dataset,x=nRF)) + geom_density() + theme_classic() + xlab("Error (nRF)") + 
  theme(legend.position = c(0.85,0.85), text = element_text(size=11), )
ggsave("./figures/discacross.pdf",width=4, height=4) 

###########################

df <- read.table("./data/discacross_hds.csv", header=F)
head(df)
colnames(df) <- c("QD","nRF","gene","rep","dataset")

df$dataset = as.factor(df$dataset)
levels(df$dataset) = c("HD-P1", "10k", "HD-P2","HD-P3","HD-P4","HD-P5")

ggplot(data=df, aes(color=dataset,x=nRF)) + geom_density() + theme_classic() + xlab("Error (nRF)") + 
  theme(legend.position = c(0.85,0.85), text = element_text(size=11), )
ggsave("./figures/discacross_hdps.pdf",width=4, height=4) 
