######################
#### Libraries #######
######################

library(adegenet)
library(diveRsity)
library(poppr)
library(Demerelate)
library(rworldmap)
library(geosphere)
library(raster)
library(sp)
library(sf)
library(rworldmap)
library(rgdal)
library(spdep)
library(rgeos)
library(dismo)
library(gbm)
library(AUC)
library(ggplot2)
library(plyr)
library(HH)
library(ltm)
library(geosphere)

###############################
#### Load in docs, set wd #####
###############################
my_dir <- "G:\\Shared drives\\Emily_Schumacher\\Oak_RAD\\arp_genind_files"
setwd(my_dir)

##convert files
arp_list <- list.files(path = my_dir, pattern = ".arp$")

for(i in 1:length(arp_list)){
  
  arp2gen(arp_list[[i]])
  
}

##list out genind files 
genind_list <- list.files(path = my_dir, pattern = ".gen$")

###read genind files
quac_qubo_genind <- list()

##load in genind files
for(j in 1:length(genind_list)){
  
  quac_qubo_genind[[j]] <- read.genepop(paste0(my_dir,"\\",genind_list[[j]]), ncode = 3)
  
}

##load in data frames 
setwd("G:\\Shared drives\\Emily_Schumacher\\Oak_RAD\\QUAC_QUBO_dataframes")

##create df list 
dataframe_list <- list.files(pattern = "_df.csv")

##save dfs 
QUAC_QUBO_df_list <- list()

for(k in 1:length(dataframe_list)){
  
  ##create a data frame list 
  QUAC_QUBO_df_list[[k]] <- read.csv(dataframe_list[[k]])
  
  ##name rownames of the genind to id individuals 
  rownames(quac_qubo_genind[[k]]@tab) <- QUAC_QUBO_df_list[[k]][1:length(rownames(quac_qubo_genind[[k]]@tab)),1]
  
}

#############################################
######### Clone Correct Code ################
#############################################
clone_test <- list()
clonetest_list <- list()
clone_index <- list()
genind_nocl_list <- list()
popr_nocl <- list()
genind_nocl_list <- list()

##loop to remove clones from genind 
for(i in 1:length(quac_qubo_genind)){
  ##try clone analyses for QUAC
  clone_test[[i]] <- as.genclone(quac_qubo_genind[[i]])
  clonetest_list[[i]] <- mlg.id(quac_qubo_genind[[i]])
  #Function to pull out individual indices where clone length greater than 1
  clone_index[[i]] <- which(sapply(clonetest_list[[i]],function(x) length(x)>1))
  
  ##remove clones from test doc
  popr_nocl[[i]] <- clonecorrect(clone_test[[i]])
  #genind2genalex(genclone2genind(popr_nocl),file="QH_clone_free.csv")
  #Create genpop and genind objects that now have no clones- GI_nocl, GP_nocl
  genind_nocl_list[[i]] <- genclone2genind(popr_nocl[[i]])
}

##Set output location
setwd("G:\\Shared drives\\Emily_Schumacher\\Oak_RAD")

##list for no clone data frame
nocl_df <- list()

##reduce data frames from clones 
for(i in 1:length(genind_nocl_list)){
  
  nocl_df[[i]] <- QUAC_QUBO_df_list[[i]][QUAC_QUBO_df_list[[i]][,1] %in% rownames(genind_nocl_list[[i]]@tab),]
}

###################################
######### Relatedness #############
###################################
QUAC_QUBO_rel <- list()
QUAC_QUBO_halfsib_names <- list()
QUAC_QUBO_halfsib_names_cleanfront <- list()
QUAC_QUBO_halfsib_names_cleanback <- list()
QUAC_QUBO_relate_ind <- list()
##create a data frame to average # of related inds 
rel_matrix <- matrix(nrow = 4, ncol = 2)

for(i in 1:length(nocl_df)){
 
  ##calculate relatedness
  QUAC_QUBO_rel[[i]] <- Demerelate(nocl_df[[i]], object = T, value = "loiselle")
  
  ##get the names of over 25% relatedness 
  QUAC_QUBO_halfsib_names[[i]] <- names(which(unlist(QUAC_QUBO_rel[[i]]$Empirical_Relatedness) > 0.25))
  
  ##clean the front of the names 
  QUAC_QUBO_halfsib_names_cleanfront[[i]] <- gsub("^.*\\.","", QUAC_QUBO_halfsib_names[[i]])
  
  ##clean the back - remove the second name 
  QUAC_QUBO_halfsib_names_cleanback[[i]] <- gsub("^.*\\_","", QUAC_QUBO_halfsib_names_cleanfront[[i]])
  
  ##now isolate to single individuals
  QUAC_QUBO_relate_ind[[i]] <- unique(QUAC_QUBO_halfsib_names_cleanback[[i]])
  
  ##save some results in a matrix
  rel_matrix[i,2] <- length(QUAC_QUBO_relate_ind[[i]])/length(rownames(genind_nocl_list[[i]]@tab))
  
  
}


