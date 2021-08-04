#########################################################################################

#install packages
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if("BSDA" %in% rownames(installed.packages()) == FALSE){install.packages("BSDA")
} else {print (paste0("'BSDA' has already been installed in library"))}

#load packages
pacman::p_load("BSDA")

##########################################################################################
################################# INPUT FILES ############################################
# 1- cvs file with genomic diversity summary statistics
# 2- structure file with the genomic data 

##########################################################################################
####################################### PLOT #############################################

setwd("~/") # set working directory

data <- read.csv("sumstat.csv", header=TRUE, sep=",")

x1<- data$Exp_Het
x2<- data$Pi
x3<- data$Fis

#tiff("boxplot_sumstat_sp.tiff", width=20, height=15, unit="cm", res=300)
pdf("boxplot_sumstat_sp.pdf")
setEPS()
postscript("boxplot_sumstat.eps")
plot.new()    
par(oma=c(0,1,1,1));
par(mar=c(5,5,1,0));
plot.window(xlim=c(0,12), ylim=c(0, 0.5))

boxplot(x1, x2, x3, data=data)

points(x= 1, y=data[1,5], col = 'slategray4', bg = 'red', cex = 2, pch=21)
points(x= 1, y=data[2,5], col = 'slategray4', bg = 'orange', cex = 2, pch=21)
points(x= 1, y=data[3,5], col = 'slategray4', bg = 'darkgreen', cex = 2, pch=21)
points(x= 1, y=data[4,5], col = 'slategray4', bg = 'blue', cex = 2, pch=21)

points(x= 2, y=data[1,8], col = 'slategray4', bg = 'red', cex = 2, pch=21)
points(x= 2, y=data[2,8], col = 'slategray4', bg = 'orange', cex = 2, pch=21)
points(x= 2, y=data[3,8], col = 'slategray4', bg = 'darkgreen', cex = 2, pch=21)
points(x= 2, y=data[4,8], col = 'slategray4', bg = 'blue', cex = 2, pch=21)

points(x= 3, y=data[1,11], col = 'slategray4', bg = 'red', cex = 2, pch=21)
points(x= 3, y=data[2,11], col = 'slategray4', bg = 'orange', cex = 2, pch=21)
points(x= 3, y=data[3,11], col = 'slategray4', bg = 'darkgreen', cex = 2, pch=21)
points(x= 3, y=data[4,11], col = 'slategray4', bg = 'blue', cex = 2, pch=21)

#axis(2, at=seq(0, 0.5, by=0.1), cex.axis=1.15, las=1);
#legend<- c("East_hist", "East_exp", "West_hist",  "West_exp")
#legend("topleft", legend = legend, cex=0.6, bty="n", ncol=3,y.intersp=1, col = c('red', 'orange', 'darkgreen', 'blue'), pch= 16)

dev.off()

##########################################################################################
###################################### t TEST ############################################

### Obs_Het ### 

#EH X EE
p1<-tsum.test(mean.x=0.19305,   s.x=0.239645572, n.x=16,
          mean.y=0.15239, s.y=0.28167357, n.y=8)
p1

#EH X WE
p2<-tsum.test(mean.x=0.19305,   s.x=0.239645572, n.x=16,
              mean.y=0.11329, s.y=0.288738636, n.y=3)
p2

#EH X WH

p3<-tsum.test(mean.x=0.19305,   s.x=0.239645572, n.x=16,
              mean.y=0.12082, s.y=0.270905888, n.y=5)

p3

#EE X WE

p4<-tsum.test(mean.x=0.15239,   s.x=0.28167357, n.x=8,
              mean.y=0.11329, s.y=0.288738636, n.y=3)

p4

#EE X WH

p5<-tsum.test(mean.x=0.15239,   s.x=0.28167357, n.x=8,
              mean.y=0.12082, s.y=0.270905888, n.y=5)

p5

#WE X WH

p6<-tsum.test(mean.x=0.11329, s.x=0.288738636, n.x=3,
              mean.y=0.12082, s.y=0.270905888, n.y=5)

p6

### Pi ### 

#EH X EE
p7<-tsum.test(mean.x=0.27075,   s.x=0.199198394, n.x=16,
              mean.y=0.18635, s.y=0.27462702, n.y=8)
p7

#EH X WE
p8<-tsum.test(mean.x=0.27075,   s.x=0.199198394, n.x=16,
              mean.y=0.11871, s.y=0.283337255, n.y=3)
p8

#EH X WH

p9<-tsum.test(mean.x=0.27075,   s.x=0.199198394, n.x=16,
              mean.y=0.14107, s.y=0.267301328, n.y=5)

p9

#EE X WE

p10<-tsum.test(mean.x=0.18635,   s.x=0.27462702, n.x=8,
               mean.y=0.11871, s.y=0.283337255, n.y=3)

p10

#EE X WH

p11<-tsum.test(mean.x=0.18635,   s.x=0.27462702, n.x=8,
               mean.y=0.14107, s.y=0.267301328, n.y=5)

p11

#WE X WH

p12<-tsum.test(mean.x=0.11871, s.x=0.283337255, n.x=3,
               mean.y=0.14107, s.y=0.267301328, n.y=5)

p12


### Fis ### 

#EH X EE
p13<-tsum.test(mean.x=0.22863,   s.x=0.468785665, n.x=16,
              mean.y=0.06689, s.y=0.313177266, n.y=8)
p13

#EH X WE
p14<-tsum.test(mean.x=0.22863,   s.x=0.468785665, n.x=16,
              mean.y=0.00866, s.y=0.147444905, n.y=3)
p14

#EH X WH

p15<-tsum.test(mean.x=0.22863,   s.x=0.468785665, n.x=16,
              mean.y=0.03592, s.y=0.242631408, n.y=5)

p15

#EE X WE

p16<-tsum.test(mean.x=0.06689,   s.x=0.313177266, n.x=8,
               mean.y=0.00866, s.y=0.147444905, n.y=3)

p16

#EE X WH

p18<-tsum.test(mean.x=0.06689,   s.x=0.313177266, n.x=8,
               mean.y=0.03592, s.y=0.242631408, n.y=5)

p18

#WE X WH

p19<-tsum.test(mean.x=0.00866, s.x=0.147444905, n.x=3,
               mean.y=0.03592, s.y=0.242631408, n.y=5)

p19

Bonferroni<- p.adjust (c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19), method = "bonferroni", n = 19)
Bonferroni

#################################################################################################################
################################################# PCA ###########################################################

setwd("~/") #set working directory

