args <- commandArgs(trailingOnly = TRUE)
pca_file <- args[1] #PCA file
exclusions <- args[2] #Subjects to remove
n_pcs <- args[3] #Number of PCs to retain
pheno_file <- args[4] #Phenotype file
pheno_name <- args[5] #Phenotype name
covar_name <- args[6] #Covariate names
studycov <- args[7] #Add study covar based on sample FIDs
outfile <- args[8] #Output filename

unlist_split2 <- function(x, colnum ,...)
{
    toret <- unlist(strsplit(x, ...) )[colnum]
    return(toret)
}
 


 
#This is a number, needs to be converted to one.
 n_pcs <- as.numeric(n_pcs)

#PCA results from LDAK
 pca <- read.table(pca_file, header=F,stringsAsFactors=F)


 
#Phenotype file. ASSUMES THAT -9 , NA, and #N/A MEAN MISSING!!
 pheno0 <- read.table(pheno_file,header=T,stringsAsFactors=F,na.strings=c("NA","#N/A","-9"))
#Exclude samples if called upon
if (exclusions != "xxxx")
{ 
 excluder <- read.table(exclusions,header=T,stringsAsFactors=F)
 pheno0_FID_IID <- paste(pheno0$FID,pheno0$IID,sep="_")
 excluder_FID_IID <- paste(excluder$FID,excluder$IID,sep="_")
 cutthem <- which(pheno0_FID_IID %in% excluder_FID_IID)
 pheno0 <- pheno0[-cutthem,]

}


 
#Name PCs simply
 names(pca)[1:2] <- c("FID","IID")
 names(pca)[3:dim(pca)[2]] <- paste("PC",1:(dim(pca)[2] - 2 ),sep="")
 
#Subset PCA file to just PCs you want
 pca <- pca[,c(1:(2+n_pcs))]
 
#List FID and IID of all subjects with PCs - This is the set of unrelated subjects
 usable <- subset(pca,select=c(FID,IID))
 
#Only take phenotyped subjects, i.e. those with phenotype values
 pheno1 <- subset(pheno0,!is.na(pheno_name))

#Only take phenotyped subjects with PCs
 pheno <- merge(pheno1,usable,by=c("FID","IID"))
 


 if (studycov == TRUE)
 {
  dummy <- pheno
  dummy$studycov_split <- sapply(pheno$FID,unlist_split2,colnum=3,split="_")
  recoded_site <-  model.matrix(~studycov_split, data=dummy)[,-1]  #Make a dummy covar, minus the intercept
  colnames(recoded_site) <- paste ("dummy", colnames(recoded_site),sep="_")
  pheno2 <- data.frame(pheno,recoded_site) #Put the dummy coded thing into the phenotype file itself
  
 }
 if (is.character(studycov)) #If study cov is a string
 {
  studycovfile <- read.table(studycov,header=T,stringsAsFactors=F)
  pheno1a <- merge(pheno,studycovfile,by=c("FID","IID"))
  
  recoded_site <-  model.matrix(~studycov, data=pheno1a)[,-1]  #Make a dummy covar, minus the intercept
  colnames(recoded_site) <- paste ("dummy", colnames(recoded_site),sep="_")
  pheno2 <- data.frame(pheno1a,recoded_site) #Put the dummy coded thing into the phenotype file itself
 }

#Write out the phenotype 
 write.table(subset(pheno2,select=c("FID","IID",pheno_name)),paste(outfile,'.pheno',sep=''),row.names=F,quote=F)
  
#Write just the covars 
 covar1 <- merge(pca,pheno2,by=c("FID","IID"))
 
 
 if (covar_name == "xxxx")
 {
  covarx <- covar1
   if (studycov != FALSE)
  {
   covar <- subset(covar1, select=c("FID","IID",names(pca)[-c(1:2)],colnames(recoded_site)))
  }
 } else { 
  covar <- subset(covar1, select=c("FID","IID",covar_name))
  if (studycov != FALSE)
  {
   covar <- subset(covar1, select=c("FID","IID",covar_name,colnames(recoded_site)))
  }
 }


 write.table(covar,paste(outfile,'.cov',sep=''),row.names=F,quote=F)
  
  
#Example of subsetting on gender
#  write.table(subset(pheno,Gender == 1, select=c(FID,IID,Pheno_fix)),"phenotypes/malptsd_nomega_20pct.pheno",row.names=F,quote=F) #  & included_m ==1 
#  write.table(subset(pheno,Gender == 2 ,  select=c(FID,IID,Pheno_fix)),"phenotypes/femptsd_nomega_20pct.pheno",row.names=F,quote=F) # & included_f == 1 
 
