#### MODIFIED FROM Andrea Thomaz (2015), Sarp Kaya (2017), Giorgia Auteri (2018) Cecilia Fiorini (2019)
#######################################################################################################

#install packages
if("plyr" %in% rownames(installed.packages()) == FALSE){install.packages("plyr")
} else {print (paste0("'plyr' has already been installed in library"))}
if("pegas" %in% rownames(installed.packages()) == FALSE){install.packages("pegas")
} else {print (paste0("'pegas' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("plyr", "pegas")

#########################################################################################################
############################################# INPUT FILES ###############################################
# 1- .vcf file from stacks output 

#########################################################################################################
#READ VCF
setwd("~/")
data <- read.table('populations.snps.vcf', header = FALSE, sep = "\t")
head(data[1:10,1:10])


loci_num<-as.numeric(sub(':.*', '', data$V3))
pos2 <- as.numeric(gsub(".*[:]([^.]+)[:].*", "\\1", data$V3))


#creates dataframe with loci ID, the variable positions and the number of individuals in each loci
new_data <- data.frame(loci_ID = loci_num,
                       pos = pos2,  
                       ind = rowSums(data[,14:length(data)] != "./."))

#SEQUENCE LENGTH
seq_len <- as.numeric(max(new_data$pos))

head(new_data)
min(new_data$pos)#should always be position 5 (first positions are the adapters)
max(new_data$pos)#should always be position 139 !!!UNLESS!!! data is mapped to reference genome, may be longer if you allowed for some insertions
table(new_data$pos)
length(unique(new_data$loci_ID))#how many loci do I have

#pdf(file="segregating-sites.pdf", 
   # width=5, 
    #height=4, 
    #pointsize=12)
par(mar = rep(2, 4))
#saving graph with frequency of variable sites along the loci
#pdf("./SNPdistr_pos140bp.pdf")
hist(new_data$pos, breaks = c(seq(1, seq_len , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites'); #for unmapped loci
#play with first abline number to determine cutoff
abline(11500, 0, col = "red")#helps to find where starts to increase toward the end, last positions have strong increase
abline(v = 5, col = "red")#helps to figure out where to cut off before increase in bad calls
#move the lines around to visualize depending on my case

# #BASE ON THE GRAPH, CHOOSE HOW MANY POSITION TO DELETE FROM THE *end*
to_del <- 16 #how many sites to delete in the end of the sequence (10 is based on the 130 I chose for the ab line above)
seq_len_cut <- seq_len - to_del
# #create a whitelist to exclude those 10 (to_del) positions from the end
whitelist <- subset(new_data, pos < seq_len_cut & pos > 6)
# #pdf("./SNPdistr_pos_cutto125bp.pdf")
hist(whitelist$pos, xlim = c(0,seq_len_cut), breaks = c(seq(-1, seq_len_cut -1 , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');

#BASE ON THE GRAPH, CHOOSE HOW MANY POSITION TO DELETE FROM THE *beginning*
#to_del <- 6 #how many sites to delete in the beginning of the sequence
#seq_len_cut <- seq_len - to_del
#create a whitelist to exclude those 10 (to_del) positions from the beginning
#whitelist <- subset(whitelist, pos > seq_len_cut)

#whitelist <- new_data
table(whitelist$pos)
hist(whitelist$pos, breaks = c(seq(1, seq_len , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites after cut');

#calculating theta for all loci
var.sites <- count(whitelist, "loci_ID")
length(var.sites$loci_ID)
max(var.sites$freq) #max variable sites in one loci; 
theta_calc <- merge(unique(whitelist[,-2]), var.sites, by = "loci_ID")
theta_calc$theta <- 0
head(theta_calc)
for (i in 1:length(theta_calc$theta)){
   theta_calc[i,4] <- theta.s(theta_calc$freq[i], theta_calc$ind[i])/seq_len
 }

#calculating the 95% quantile to exclude loci extremely variable
quant <- quantile(theta_calc$theta, probs=0.95) #set the value to be variability threshold
quant 
#pdf("./theta125bp.pdf")
hist(theta_calc$theta)
abline(v = quant, col="red")
#dev.off()

#what is the maximum number of mutations in a loci
max(theta_calc$freq) # max theta (variable positions before ble positions in one loci)
x <- subset(theta_calc, theta < quant)
max(x$freq) # max theta after, make sure is realistic for a 140 bp sequence
#think about what mutation rate the spp might have

#saving whitelist for re-run populations in stacks (write blacklist and subtract it from whitelist)
blacklist <- subset(theta_calc, theta > quant)[,1]
#write.table(blacklist, file="blacklist.txt", sep = '\n', row.names = F, col.names = F)
#above blacklist would only be for highly variable loci

#removes the blacklist from the whitelist and write off white list
whitelist$blacklist <- match(whitelist$loci_ID, blacklist, nomatch = 0)
whitelist_final <- subset(whitelist, blacklist == 0)
length(unique(whitelist_final$loci_ID)) #final number of unique loci, this is the number I need to get out with "write random loci"
length(whitelist_final$loci_ID) #final number of snps
write.table(whitelist_final[,1:2], file="whitelist.txt", sep = '\t', row.names = F, col.names = F)

#re-enable scientific notation
options(scipen = 0)
dev.off()


