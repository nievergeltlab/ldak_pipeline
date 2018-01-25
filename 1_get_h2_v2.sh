### LDAK SNP based heritability estimates


#Set working directory 
wd=/home/cnieverg/gctam/ldak_test

#Call into it
 cd $wd

#Put location of LDAK binary here
 ldak="$wd"/ldak5.linux
#Location of multithread ldak
 ldakf="$wd"/ldak5.linux.fast
 
#Put location of PLINK1.9 binary here
 p_loc="$wd"/plink
 
#Name of your study (all outputs will have this name)
 sname=eurtest

#Set and make directories for outputs. Copy and paste this part, just leave as defaults. 


 #1: Location of merged genotypes . 
 outdir_mgeno="$wd"/merged_genotypes
 #2A1 and 2A2: Location of matrix for calculation of relatedness and pcs
 outdir_relatedness="$wd"/relatedness/bychr 
 #3A: Location of weight files
 outdir_weights="$wd"/weights
 #3B Kinship matrix location (by chromosome)
 kinship_dir="$wd"/kinships
 #3C:Merged GRM location
 combkinship_dir="$wd"/combined_kinships
 #4: LDAK results
 ldak_resultsdir="$wd"/results
 if [ ! -e $outdir_mgeno ]
 then
  mkdir $outdir_mgeno
 fi
 
 if [ ! -e $outdir_relatedness ]
 then
  mkdir -p $outdir_relatedness
 fi

 if [ ! -e $outdir_weights ]
 then
  mkdir  $outdir_weights
 fi

 if [ ! -e $kinship_dir ]
 then
  mkdir  $kinship_dir
 fi

 if [ ! -e $combkinship_dir ]
 then
  mkdir  $combkinship_dir
 fi
  if [ ! -e $ldak_resultsdir ]
 then
  mkdir  $ldak_resultsdir
 fi
 
### Step 1: Merge multiple datasets, per chromosome

##Assumes that studies have already been filtered by ancestry, and that most studies exist in their own separate files.
#Note: Generally it's going to be easier to analyze ancestries separately. If you want to analyze an ancestry group, it's easy to just put another copy of this code in a different location and change paths as necessary.

#Put the genotype files that you want to merge into a folder. They should be split by chr and have the filename end with _chromosome. e.g. pts1_22 is dataset pts1, chromosome 22
 starting_geno_dir="$wd"/genotypes/
 
##Set parameters for quality filters genotype data AFTER merge
#only retain markers genotyped in this % of people
 gtthresh=0.95
#maf: filter on this MAF
 mafthresh=0.01

 #PUT IN FILTERS IN CODE 
#Merge all datasets together (one job per chromosome)

#Note: By default this will export only SNPs with ACGT alleles (As of writing, LDAK does not like alleles with more than 1 character)
#If you don't want to do this, set -a to no. But to avoid errors I recommend just leaving this as is
 qsub -t1-22  -l walltime=00:15:00 scripts/1_process_data_merge.qsub -d $wd -e errandout/ -o errandout/ -F "-g $starting_geno_dir -p $p_loc -n $sname -o $outdir_mgeno -a yes -q $gtthresh -m $mafthresh" 
 
###Step 2: Identify related subjects

#We are doing this because SNP based heritability analysis requires related people to be removed

##2A1: Make kinship matrix for each LD pruned chromosome.

 #Program parameter settings
 ldakmodel=ldak #If value set to GCTA, no weights are used and power is set to -1. Otherwise power will be at user setting and weights will be used.
 power=-0.25 #Set power for LDAK GRM calculation. Speed recommends -0.25. GCTA is -1.
 pruning=yes #LD prune data prior to making GRM. Should be on for relatedness checking.
 windowkb=1000 #Window size for LD pruning.
 r2=0.2 #r2 for LD pruning.
 genodir_relatedness=$outdir_mgeno #Directory of genotypes used in GRM calculation.
 useweights=no #Use weights in computation of GRM. Set to 'no' if you don't want to use weights. For this step, don't use weights.

#If output folder doesnt exist yet, just make it
 qsub -t1-22 -l walltime=00:15:00 scripts/2a_identify_related.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldak -l ldakmodel -p $power -q $pruning -w $windowkb -r $r2 -n $sname -g $genodir_relatedness -u $useweights -o $outdir_relatedness "

##2A2: Merge chromosomal kinship matrices, identify related subjects to remove. Calculate PCs to include in analysis (probably redundant due to GRM, but OK to include)

 #Note: LDAK manual has a step to filter relatedness by dataset, then all data together. Not included here..
 #Note: LDAK manual suggested filtering out snps with p < 1e-20. I don't do this step, becausse I didn't find any. If you have markers like that, filter them out early on
 
 #Program parameter settings
 #Note: xxxx is what optional parameters should be set to if you want to ignore them
 rel_cutoff=xxxx #For relatedness, by default LDAK filters using the threshold c, where -c is smallest observed kinship. To set this value manually, set -r to a value from 0 to 1.
 pheno_file=xxxx #Phenotype file (phenotyped subjects are retained preferentially. If you specify a pheno file, when it finds a related pair, it will prefer to keep the phenotyped subject. 
 keepfile=xxxx #keep flag for retaining custom subset of people. If you specify a keep file, subjects NOT in this file will be removed from analysis
 kinship_matrix_location=$outdir_relatedness #Output of kinship analysis (PCAs, grm, etc)
 outdir_kinship_grm="$wd"/relatedness
 qsub -l walltime=00:15:00 scripts/2a2_merge_related_grm.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldak -p $pheno_file -k $keepfile -r $rel_cutoff -g  $kinship_matrix_location -n $sname -o $outdir_kinship_grm" 

 
##2B: Prepare phenotype and covariates data.
 
 #Program parameter settings
 pca_file="$outdir_kinship_grm"/"$sname"_relatedness_pruned.vect  #PCA file from last analysis
 phenofile=dummy.pheno #Phenotpye file should be blankspace delimited. It should have columns 1 and 2 be FID and IID. It must have a header. Put phenotype and covariates in here
 phenoname=ptsd #Name of phenotype in pheno file
 covar_names=xxxx #Write covariate names here, if any. Delimit using commas e.g. C1,C2,C3 
 npcs=5 #Number of PCs to include. If none are included, PCA file will just be used to determine unrelated people to include
 phenooutname="$sname"_ldak #output filename of pheno and covs
 Rscript scripts/make_ldak_pheno.r $pca_file $npcs $phenofile $phenoname $covar_names $phenooutname
  
  
###Step 3: Make GRM for analysis

#The previous step was an LD pruned matrix used for identifying relateds. This matrix will be based on the whole data and incorporate the LDAK or GCTA model.
#Furthermore this matrix will incorporate weights

#3A: Split genome into sections to weight. Calculate weightings for each section.

#Note: Manual suggests doing this by cohort to prevent genotyping errors from influencing weights. However,this is impractical because of our number of cohorts and their wildly different sizes, so I don't do it!
 
 weightings_keepflag="$outdir_relatedness"/"$sname"_relatedness_pruned.keep #Put a list of subjects if you want to calculate weights in just a subset of subjects, namely the pruned subject list

 qsub -t1-22 -l walltime=00:05:00 scripts/3_ldak_weights.qsub -d $wd -e errandout/ -o errandout/  -F "-e $ldak -g $outdir_mgeno -n $sname -k $weightings_keepflag -o $outdir_weights "

#Join weightings (Should be very quick,maybe not even needed to run job script)
 for chr in {1..22} ; do $ldak --join-weights  "$outdir_weights"/sections_$chr --bfile "$outdir_mgeno"/"$sname"_"$chr" ; done ; 
 
 #If it is too much of a computational burden, echo the above line to a file and submit a job.. e.g. qsub  -l walltime=00:10:00 -e errandout -o errandout -d $wd weights_join.qsub
 
#3B: Calculate kinship for each chromosome
 
 ldakmodel=ldak #If value set to gcta, no weights are used and power is set to -1. If set to 'unweighted', no weights will be used (power can be set with power). Otherwise, weights will be used
 power=-0.25 #Set power for LDAK GRM calculation. Speed recommends -0.25. GCTA is -1.

 qsub -t1-22 -l walltime=00:05:00 scripts/3B_ldak_kinship_v2_alt.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldak -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o $kinship_dir -p $power -w $outdir_weights"

#3C: Merge GRMs into entire autosome GRM (inputs: a GRM list and output file name)

#Program parameters:
 biascheck=yes #Set to yes if you want to calculate GRMs for 1/4 segments of the genome, to be used for estimating bias
 qsub -l walltime=00:05:00 scripts/3C_merge_kinships.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldak -k $kinship_dir -o $combkinship_dir -f kinships -b yes "

##4A: Heritablity estimates

#Program parameters:
 phenotype="$sname"_ldak.pheno  #Phenofile file location
 covar="$sname"_ldak.cov #Covariate file location
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned.keep  #If you want to prune extra subjects, beyond what was already done with the relatedness matrix and built into the phenotype file, do so here. E.g. for sex stratified analyses
 outname_append=sens #Extra output names to append
 prevalence=xxxx #Extra command to set prevalence, so h2 is on the right scale
 
 qsub -l walltime=00:05:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldak  -p $phenotype -c $covar -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a $outname_append -b yes -r $prevalence"

#Once this is done we can estimate the bias due to cryptic relatedness:
 grep Her_K "$ldak_resultsdir"/"$outname_append"quad{A,B,C,D,ALL}.reml | awk '(NR<=4){sum+=$2}(NR>4){sum2+=$2}END{I=(sum-sum2)/3;print "Inflation:", I, "=", I/sum2*100,"%"}'

#Now we can do the actual estimate of h2 snp
 qsub -l walltime=00:05:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p $phenotype -c $covar -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a $outname_append -b no -r $prevalence"
 grep Her_K "$ldak_resultsdir"/"$outname_append"full.reml

  

 
