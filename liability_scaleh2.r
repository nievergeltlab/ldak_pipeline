args <- commandArgs(trailingOnly = TRUE)
h2 <- args[1] #Observed scale h2
seh2 <-  args[2] # SE of h2
pop_prev <- args[3] #population prevalence
sample_prev <- args[4] #Sample prevalence (optional)
pheno_file <- args[5] #Phenotype (optional)

h2 <- as.numeric(h2)
seh2 <- as.numeric(seh2)
pop_prev <- as.numeric(pop_prev)

#If sample prevalence is given
if(sample_prev != "xxxx")
{
 sample_prev <- as.numeric(sample_prev)
}

#If there is a phenotype file supplied, assume that it has a header and that col 3 is the phenotype, coded as 1= control, 2=  case
if(pheno_file != "xxxx")
{
 pheno <- read.table(pheno_file,header=T,stringsAsFactors=F)
 sample_prev <- length(which(pheno[,3] == 2)) /  length(which(pheno[,3] %in% c(1,2)))
}

#From equation 23 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/ Estimating Missing Heritability for Disease from Genome-wide Association Studies
K <- pop_prev
P <- sample_prev
zv <- dnorm(qnorm(K))

h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
var_h2_liab <- ( seh2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2) ^2
p_liab <- pchisq(h2_liab^2/var_h2_liab,1,lower.tail=F)

print(c(h2_liab,sqrt(var_h2_liab),p_liab))