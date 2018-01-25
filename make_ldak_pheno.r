args <- commandArgs(trailingOnly = TRUE)
pca_file <- args[1] #PCA file
n_pcs <- args[2] #Number of PCs to retain
pheno_file <- args[3] #Phenotype file
pheno_name <- args[4] #Phenotype name
covar_name <- args[5] #Covariate names
outfile <- args[6] #Output filename

#This is a number, needs to be converted to one.
 n_pcs <- as.numeric(n_pcs)

#PCA results from LDAK
 pca <- read.table(pca_file, header=F,stringsAsFactors=F)

#Phenotype file. ASSUMES THAT -9 , NA, and #N/A MEAN MISSING!!
 pheno <- read.table(pheno_file,header=T,stringsAsFactors=F,na.strings=c("NA","#N/A","-9"))

 
#Name PCs simply
 names(pca)[1:2] <- c("FID","IID")
 names(pca)[3:dim(pca)[2]] <- paste("PC",1:(dim(pca)[2] - 2 ),sep="")
 
#Subset PCA file to just PCs you want

 pca <- pca[,c(1:(2+n_pcs))]
 
#List FID and IID of all subjects with PCs - This is the set of unrelated subjects
 usable <- subset(pca,select=c(FID,IID))
 
#Only take phenotyped subjects, i.e. those with phenotype values
 pheno <- subset(pheno,!is.na(pheno_name))

#Only take phenotyped subjects with PCs

 pheno <- merge(pheno,usable,by=c("FID","IID"))
 
#Write out the phenotype 
 write.table(subset(pheno,select=c("FID","IID",pheno_name)),paste(outfile,'.pheno',sep=''),row.names=F,quote=F)
  
#Write just the covars 
 if (covar_name != "xxxx")
 {
  
  covar1 <- merge(pca,pheno,by=c("FID","IID"))
  covar <- subset(covar1, select=c("FID","IID",covar_name))
 } else { covar <- pca }
 write.table(covar,paste(outfile,'.cov',sep=''),row.names=F,quote=F)
  
  
#Example of subsetting on gender
#  write.table(subset(pheno,Gender == 1, select=c(FID,IID,Pheno_fix)),"phenotypes/malptsd_nomega_20pct.pheno",row.names=F,quote=F) #  & included_m ==1 
#  write.table(subset(pheno,Gender == 2 ,  select=c(FID,IID,Pheno_fix)),"phenotypes/femptsd_nomega_20pct.pheno",row.names=F,quote=F) # & included_f == 1 
 
