require(ggplot2)
require(reshape2)
require(scales)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

setwd("~/Workspace/btol/")



##########dissertation figures##########

sz=500
i=nRF
ggplot(dfm[dfm$size == sz & dfm$mc == "mc2"  & dfm$variable == "true",], aes(x = variable, y = value, fill=variable)) + 
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_violin(scale = "count", trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75) ) + #facet_wrap(.~size)+
  stat_summary()+theme_classic() +
  #scale_y_continuous(trans="log", breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = cbPalette[c(1,2,3)]) + xlab("")+ylab(sprintf("Gene Tree Discordance (%s)", i)) +
  #coord_cartesian(ylim=c(0.01, 0.8)) + 
  guides(fill=F)


ggsave("./figures/violin_gene_discordance_nRF_trueonly.pdf",width=2, height=4)


dfall <- read.delim("data/stee_full_all.csv",fill = T,na.strings = "", header=F)
head(dfall)

i="nRF"
sz="100"

colnames(dfall) <- c("rep",  "size", "uDanceq", "uDance", 
                     "FastTreeq", "FastTree2", "concatq", "concat","mc")
dfrf <- dfall[,c(1,2,4,6,8,9)]

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

ggplot(dfm[ dfm$mc == "mc2"  & dfm$variable == "uDance",], aes(x = variable, y = value)) + 
  facet_wrap(.~size)+
  #geom_text(data=z,aes(y=0.01, x=variable, label=con))+
  geom_point(aes(color=variable)) + #facet_wrap(.~size)+
  stat_summary( geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75)+ 
  stat_summary( geom = "errorbar", width = .25, size=0.3)+
  theme_classic() +
  #scale_y_continuous() + #, breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
  coord_cartesian(ylim=c(0.009, 0.10)) + 
  guides(color=F) #+ coord_cartesian(ylim = c(0,0.3))

ggsave(sprintf("./figures/violin_stee_full_%s_%s_udance.pdf", i,sz),width=2.5, height=4)


ggplot(dfm[ dfm$mc == "mc2"  & dfm$variable %in% c("uDance", "FT2-Astral"),], aes(x = variable, y = value)) + 
  geom_text(data=z[z$variable %in% c("uDance", "FT2-Astral" ) & z$mc == "mc2", ],aes(y=0.01, x=variable, label=con))+
  facet_wrap(.~size)+
  
  geom_point(aes(color=variable)) + #facet_wrap(.~size)+
  stat_summary( geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75)+ 
  stat_summary( geom = "errorbar", width = .25, size=0.3)+
  theme_classic() +
  #scale_y_continuous() + #, breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
  coord_cartesian(ylim=c(0.009, 0.10)) + 
  guides(color=F) #+ coord_cartesian(ylim = c(0,0.3))

ggsave(sprintf("./figures/violin_stee_full_%s_%s_udance_ft.pdf", i,sz),width=4, height=4)

ggplot(dfm[ dfm$mc == "mc2" ,], aes(x = variable, y = value)) + 
  geom_text(data=z[ z$mc == "mc2", ],aes(y=0.01, x=variable, label=con))+
  facet_wrap(.~size)+
  
  geom_point(aes(color=variable)) + #facet_wrap(.~size)+
  stat_summary( geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75)+ 
  stat_summary( geom = "errorbar", width = .25, size=0.3)+
  theme_classic() +
  #scale_y_continuous() + #, breaks = c(0.01, 0.1, 0.8))+
  theme(text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(values = cbPalette[c(2,3,4)]) + xlab("")+ylab(sprintf("Species Tree Estimation Error (%s)", i)) +
  coord_cartesian(ylim=c(0.009, 0.10)) + 
  guides(color=F) #+ coord_cartesian(ylim = c(0,0.3))

ggsave(sprintf("./figures/violin_stee_full_%s_%s_all.pdf", i,sz),width=5, height=4)

######################################
