#################################################################################################################

#install packages
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("fossil" %in% rownames(installed.packages()) == FALSE){install.packages("fossil")
} else {print (paste0("'fossil' has already been installed in library"))}
if("vegan" %in% rownames(installed.packages()) == FALSE){install.packages("vegan")
} else {print (paste0("'vegan' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("adegenet", "dartR", "fossil", "vegan")


#################################################################################################################
############################################# INPUT FILES #######################################################
# 1- structure file with the genomic data 
# 2- cvs file with the geographic coordinates of each individual

##################################################################################################################
############################################ pairwise FST ########################################################

setwd("~/") #set working directory

input = read.structure('pop.str', n.ind = 32, 
                       n.loc = 7088, col.lab = 1, col.pop = 2, 
                       onerowperind = FALSE, row.marknames = 0, NA.char = 0, ask = FALSE)
input2 <- gi2gl (input) 

### fst based on Weir and Cocker- ham (1984) ###

FST_stat = gl.fst.pop(input2, nboots = 1000, percent = 95, nclusters = 1)

FST<-FST_stat$Fsts
write.table(FST,"populations_fst.csv", sep=",")

Pvalues<-FST_stat$Pvalues

Pvalues<-c(0.000,0.345,0.000, 0.064, 0.000,  0.000) 

#bonferroni correction

bonferroni<- p.adjust(Pvalues, method = "bonferroni", n = length(Pvalues))

###################################################################################################################
############################################ MANTEL TEST ##########################################################

coord <- read.csv("pop_coord.csv", header = T, sep=",")
coord2 <- coord[,2:3]
mDIST <- earth.dist(coord2, dist = TRUE)
mDIST <- as.dist(mDIST)

#FST<- read.csv("populations_fst.csv", header = T, row.names = 1, sep=",")

fst_cor <- FST/(1-FST)

mt <- mantel(mDIST, fst_cor, method="pearson", permutations = 1000)


#Perform a sequential population dropout procedure, in which the test was repeated 
#excluding one population at time, in order to confirm that the results were robust 

pop_counter= 4

for(i in 1:pop_counter) {
  temp_coord <- coord2[-i,]
  mDIST_temp <- earth.dist(temp_coord, dist = TRUE)
  mDIST_temp <- as.dist(mDIST_temp)
  fst_cor_temp <- fst_cor[-i,-i]
  fst_cor_temp <- as.dist(fst_cor_temp)
  mt <- mantel(mDIST_temp, fst_cor_temp, method="pearson")
  print(i)
  print(mt)
}

fst_cor <- as.dist(fst_cor)

## Plot Mantel test

tiff("fst_mantel.tiff", width=50, height=40, unit="cm", res=300)
setEPS()
postscript("fst_plot.eps")
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(100,250), ylim=c(0, 0.15));
points(mDIST[6],fst_cor[6], col = 'slategray4', bg = 'lightblue', cex = 2, pch=21)#chip_sch
points(mDIST[5],fst_cor[5], col = 'slategray4', bg = 'black', cex = 2, pch=21)#chip_men
points(mDIST[4],fst_cor[4], col = 'slategray4', bg = 'green', cex = 2, pch=21) #men_Sch
points(mDIST[3],fst_cor[3], col = 'slategray4', bg = 'purple', cex = 2, pch=21) #chey_Chip
points(mDIST[2],fst_cor[2], col = 'slategray4', bg = 'yellow', cex = 2, pch=21) #chey_Sch
points(mDIST[1],fst_cor[1], col = 'slategray4', bg = 'orange', cex = 2, pch=21) #chey_Men
abline(lm(fst_cor ~ mDIST), col = "gray30", lwd = 2)
axis(1, at=seq(50, 350, by=50), cex.axis=1);
axis(2, at=seq(0, 0.3, by=0.05), cex.axis=1, las=1);
mtext(side=1, text="Distance (Km)",line=2.5, cex=1)
mtext(side=2, text="FST/(1-FST)", line=2.8, cex=1)
legend('topleft', legend = "r = -0.614 , p=value = 0.125", bty = 'n', cex = 1)
dev.off()


