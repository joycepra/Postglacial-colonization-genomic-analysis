############################################### MODIFIED FROM Auteri and Knowles (2020) ##################################################################
##########################################################################################################################################################

#install packages
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("adegenet")

##########################################################################################################################################################
################################################################ INPUT FILES #############################################################################
# 1- structure file with the genomic data 

##########################################################################################################################################################
################################################################# Projected PCA ##########################################################################

setwd("~/") #set working directory

str1 <- read.structure('pop.str', n.ind = 32, 
                       n.loc = 7088, col.lab = 1, col.pop = 2, 
                       onerowperind = FALSE, row.marknames = 0, NA.char = 0, ask = FALSE)

gind <- scaleGen (str1, NA.method = "mean", scale = FALSE, center = FALSE)

#Split Data
rownames(gind)                                                              
dfgind <- as.data.frame(as.matrix(unlist(gind), nrow=(length(gind[,1])), byrow=T),stringsAsFactors=FALSE) #convert to data frame


hist <- dfgind[1:21,] # This will be the group which is projected onto the first PCA
exp <- dfgind[22:32,] #This will be whatever group the pca is first calculated with

#Check for zero variance columns 
zero_var <- which(apply(hist, 2, var)==0) #make a vector with all zero-variance columns
length(zero_var) #how many sites have no variance
histVar <- hist[,!names(hist) %in% names(zero_var)] #make a subset of the data with only variable sites. Do not do this if there are 0 non-variant sites. 

#Prepare second set of data for pca projection
#Remove sites with zero-variability in the first group from the second group. Do not do this if your first group (above) had 0 non-variant sites
exp <- exp[ , (names(exp) %in% names(histVar))] 

################################################################## Standard PCA ###########################################################################
#Make calculate PCA using first group
pcas <- prcomp(histVar, center = TRUE, scale = TRUE)
#Examine the PCA
summary(pcas) 
#Graph of proportion of data explained by each axis
plot(pcas, type = "l", main = "Proportion of data explained by each PCA axis for sources")

#Apply center and scale from first pca to the new group
pca2 <- scale(exp, pcas$center, pcas$scale) %*% pcas$rotation 
#Examine pca for second set of data. Figure out SDs then use these to calculate proportion of variance 
SDs <- rep(NA,length(colnames(pca2)))    #create vector for empty SD values
for (i in 1:length(colnames(pca2))){  # a loop to fill vector with SDs for each PCA axis
  SDs[i] <- sd(pca2[,i])
}
PropVars <- rep(NA,length(colnames(pca2)))  #Get the proportion of varinaces explained by each axis
for (i in 1:length(colnames(pca2))){     #A loop to fill vector with proportion explained by each axis
  PropVars[i] <- ((sd(pca2[,i]))^2)/sum(SDs^2)
}
sum(PropVars)       #Proportions of all variances should sum to 1
PropVars      #list proportion of variance explained by each axis for group 2
plot(PropVars, type = "l", main = "Proportion of data explained by each PCA axis for Descendents")

#Subset first three axes for group 1
grp1ax1_s1 <- pcas$x[1:16,1]
grp1ax1_s2 <- pcas$x[17:21,1]
grp1ax2_s1 <- pcas$x[1:16,2]
grp1ax2_s2 <- pcas$x[17:21,2]


#Subset first three axes for group 2
grp2ax1_s1 <- pca2[1:3,1]
grp2ax1_s2 <- pca2[4:11,1]
grp2ax2_s1 <- pca2[1:3,2]
grp2ax2_s2 <- pca2[4:11,2]


#Combine for scaling purposes
pc1 <- c(grp1ax1_s1,grp1ax1_s2, grp2ax1_s1, grp2ax1_s2)
pc2 <- c(grp1ax2_s1,grp1ax2_s2, grp2ax2_s1, grp2ax2_s2)


final_scores<- data.frame(cbind(pc1,pc2))

pops<- as.factor(c("EH", "EH","EH","EH","EH","EH","EH","EH",
                   "EH","EH","EH","EH","EH","EH","EH","EH",
                   "WH", "WH","WH","WH","WH", "WE","WE","WE","EE",
                   "EE","EE","EE","EE","EE","EE","EE"))

#####PCA 1 and 2
#tiff("pca_exp_projected.tiff", width=20, height=15, unit="cm", res=300)
setEPS()
postscript("pca_exp_projected.eps")
#pdf("pca_exp_projected_ELIP.pdf")
quartz.options(height=10, width=12, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(5,5,5,1));
plot.window(xlim=c(-100,120), ylim=c(-80, 140));
points(grp1ax1_s1,grp1ax2_s1, col = "#0000ff", bg = "#0000ff", cex = 2, pch=17)
points(grp1ax1_s2,grp1ax2_s2, col = "#fc6100", bg = "#fc6100", cex = 2, pch=16)
points(grp2ax1_s1,grp2ax2_s1, col = "#9f3d00", bg = "#9f3d00", cex = 2, pch=15)
points(grp2ax1_s2,grp2ax2_s2, col = "#fcbd9e", bg = "#fcbd9e", cex = 2, pch=18)

axis(1, at=seq(-100, 120, by=70), cex.axis=1.15);
axis(2, at=seq(-80, 100, by=60), cex.axis=1.15, las=1);

mtext(side=1, text='PC1(6.77%)',line=2.5, cex=1)
mtext(side=2, text='PC2(6.13%)', line=2.8, cex=1)
#legend<- c("East_hist", "East_exp", "West_hist",  "West_exp")
#legend("topleft", legend = legend, cex=1, bty="n", ncol=2, y.intersp=0.5, x.intersp=0.2, col = c('red', 'orange', 'darkgreen', 'blue'), pch= 16)

s.class(final_scores, fac=pops, pch= c(17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,16,16,16,16,15,
                                       15,15, 18,18,18,18,18,18,18,18), col=c("#0000ff", "#fc6100", "#9f3d00", "#fcbd9e"), 
                                       axesel=FALSE, addaxes =FALSE, cstar=0, cpoint=2, add.plot=TRUE)
dev.off()
