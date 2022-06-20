#Libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(reshape2)
library(cowplot)
library(ggpmisc)

setwd("/home/pablo/ilmode/ilmode/morph")
file = "aF0xA0S.txt"
DA <- read.table(file)

DC <- subset(DA, V2>2)
DCS <- subset(DC, V9<=12)
DCL <- subset(DC, V9>=13)

#DCS <- subset(DC, V11<=12)
#DCL <- subset(DC, V11>=13)
r <- cor.test(DA$V4, DA$V5,  method = "pearson")
print (r)

r <- cor.test(DC$V4, DC$V5,  method = "pearson")
print (r)

r <- cor.test(DCS$V4, DCS$V5,  method = "pearson")
print (r)

r <- cor.test(DCL$V4, DCL$V5,  method = "pearson")
print (r)

cat(" DA inital ", mean(DA$V4), sd(DA$V4) )
cat(" DA final ", mean(DA$V5), sd(DA$V5) )

cat(" DC inital ", mean(DC$V4), sd(DC$V4) )
cat(" DC final ", mean(DC$V5), sd(DC$V5) )

cat(" DCS inital ", mean(DCS$V4), sd(DCS$V4) )
cat(" DCS final ", mean(DCS$V5), sd(DCS$V5) )

cat(" DCL inital ", mean(DCL$V4), sd(DCL$V4) )
cat(" DCL final ", mean(DCL$V5), sd(DCL$V5) )



bf<-ggplot() +
  #  geom_point(data = DA, aes(x = V2, y = V3), size=1, shape=1) +
  geom_point(aes(x = DCS$V4, y = DCS$V5, color= 'Shorter loops (10-12 residues)'), size=2, shape=2) +
  geom_point(aes(x = DCL$V4, y =DCL$V5, color = 'Longer loops(13-15 residues) '),   size =2,  shape=17 ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.8)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,12)) + 
  xlab("Initial backbone RMSD") + ylab("Final backbone RMSD") +
  theme_bw()+
  theme(plot.title =element_blank(), legend.position="top", legend.title = element_blank() )  

print(bf)


