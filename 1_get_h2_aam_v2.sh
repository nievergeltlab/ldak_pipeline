### LDAK SNP based heritability estimates
#Note fix mrs data, aam files are in there twice, but once saved as 'eur'
#I've pulled out the AAM datasets into a folder
# for study in mrsc safr dnhs gsdc fscd cogb nss1 nss2 pts1 stro psy3 meg2 gtpc comc grac wrby
# do
 # tar xvf /archive/maihofer/dac/"$study"_aam_cogbgfile_v1.tar --strip-components=5
# done

# for study in mrsc safr dnhs gsdc fscd cogb nss1 nss2 pts1 stro psy3 meg2 gtpc comc grac wrby
# do
 # tar zxvf /archive/maihofer/dac/"$study"_qcbfile_v1.tgz --wildcards "*.mds_cov" --strip-components=2

# done

# tar xvzf /archive/maihofer/gsdc_qc_v2_mar12_2017.tgz   --strip-components=3 --wildcards "*.mds_cov" 
# tar xvzf /archive/maihofer/mrsc_qc_v2_mar12_2017.tgz   --strip-components=3 --wildcards "*.mds_cov" 
# tar xvzf /archive/maihofer/pts1_qc_v2_mar12_2017.tgz   --strip-components=3 --wildcards "*.mds_cov" 

#cat mds_cov/* | awk '{print $1,$2}' > aam_analyzed.subjects

#Set working directory 
wd=/home/maihofer/ldak

#Call into it
 cd $wd
 
#Are you working on LISA? Write yes if you are, no if you are not
lisa=yes


#If not, make a dir for scratch files
if [ $lisa != "yes" ]
then
 if [ ! -e "$wd"/scratch ]
 then
  mkdir "$wd"/scratch
 fi
 lisa="$wd"/scratch
fi
 

#Put location of LDAK binary here
 ldak="$wd"/ldak5.linux
#Location of multithread ldak
 ldakf="$wd"/ldak5.linux.fast
 
#Put location of PLINK1.9 binary here
 p_loc="$wd"/plink
 
#Name of your study (all outputs will have this name)
 sname=aamptsd

#Set and make directories for outputs. Copy and paste this part, just leave as defaults. 

 if [ ! -e "$wd"/errandout ]
 then
  mkdir "$wd"/errandout
 fi
 
 #1: Location of merged genotypes . 
 outdir_mgeno="$wd"/merged_genotypes
 #2A1 and 2A2: Location of matrix for calculation of relatedness and pcs
 outdir_relatednessa="$wd"/relatedness
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
 if [ ! -e snplist ]
 then
  mkdir  snplist
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
 qsub -t1-22  -l walltime=00:15:00 scripts/1_process_data_merge.qsub -d $wd -e errandout/ -o errandout/ -F "-g $starting_geno_dir -p $p_loc -n $sname -o $outdir_mgeno -a yes -q $gtthresh -m $mafthresh -x $lisa" 

 #Calculate for snplist at different mafs, cutoffs
 
 for chr in {1..22}
 do
  $p_loc --bfile "$outdir_mgeno"/"$sname"_"$chr"  --maf 0.01 --geno 0.05 --write-snplist --out snplist/"$sname"_maf01geno95_"$chr" --memory 2000
  $p_loc --bfile "$outdir_mgeno"/"$sname"_"$chr"  --maf 0.10 --geno 0.05 --write-snplist --out snplist/"$sname"_maf10geno95_"$chr" --memory 2000
 done
 
 
###Step 2: Identify related subjects

#We are doing this because SNP based heritability analysis requires related people to be removed

##2A1: Make kinship matrix for each LD pruned chromosome.

 #Program parameter settings
 ldakmodel=gcta #ldak #If value set to GCTA, no weights are used and power is set to -1. Otherwise power will be at user setting and weights will be used.
 power=-1 #-0.25 #Set power for LDAK GRM calculation. Speed recommends -0.25. GCTA is -1.
 pruning=yes #LD prune data prior to making GRM. Should be on for relatedness checking.
 windowkb=1000 #Window size for LD pruning.
 r2=0.2 #r2 for LD pruning.
 genodir_relatedness=$outdir_mgeno #Directory of genotypes used in GRM calculation.
 useweights=no #Use weights in computation of GRM. Set to 'no' if you don't want to use weights. For this step, don't use weights.
 snplist="$wd"/snplist/"$sname"_maf01geno95
 
#If output folder doesnt exist yet, just make it
 qsub  -t1-8 -l walltime=01:05:00 scripts/2a_identify_related.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l ldakmodel -p $power -q $pruning -w $windowkb -r $r2 -n $sname -g $genodir_relatedness -u $useweights -o $outdir_relatedness -x $lisa -z $snplist "
 qsub -t9-22 -l walltime=00:35:00 scripts/2a_identify_related.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l ldakmodel -p $power -q $pruning -w $windowkb -r $r2 -n $sname -g $genodir_relatedness -u $useweights -o $outdir_relatedness -x $lisa -z $snplist "
 
##2A2: Merge chromosomal kinship matrices, identify related subjects to remove. Calculate PCs to include in analysis (probably redundant due to GRM, but OK to include)

 #Note: LDAK manual has a step to filter relatedness by dataset, then all data together. Not included here..
 #Note: LDAK manual suggested filtering out snps with p < 1e-20. I don't do this step, becausse I didn't find any. If you have markers like that, filter them out early on
 
 #Program parameter settings
 #Note: xxxx is what optional parameters should be set to if you want to ignore them
 rel_cutoff=xxxx #For relatedness, by default LDAK filters using the threshold c, where -c is smallest observed kinship. To set this value manually, set -r to a value from 0 to 1.
 pheno_file=xxxx #Phenotype file (phenotyped subjects are retained preferentially. If you specify a pheno file, when it finds a related pair, it will prefer to keep the phenotyped subject. 
 keepfile="$wd"/aam_analyzed.subjects #keep flag for retaining custom subset of people. If you specify a keep file, subjects NOT in this file will be removed from analysis
 kinship_matrix_location=$outdir_relatedness #Output of kinship analysis (PCAs, grm, etc)
 outdir_kinship_grm="$wd"/relatedness
 qsub -l walltime=00:15:00 scripts/2a2_merge_related_grm.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -p $pheno_file -k $keepfile -r $rel_cutoff -g  $kinship_matrix_location -n $sname -o $outdir_kinship_grm -x $lisa" 

 #Recalcualte relatedness PCs in just a set pool of subjects
  max_rel=.2
  append_name=_allgwas
  qsub -l walltime=00:20:00 scripts/2a3_calcualte_pcs_only.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -p $pheno_file -k $keepfile -r $max_rel -g  $kinship_matrix_location -n $sname -o $outdir_kinship_grm -x $lisa -z $append_name" 

  
##2B: Prepare phenotype and covariates data.
 #NOte: add study covariate as you did in the Raj analysis
 #Program parameter settings
 echo "FID IID PTSD" > header_pheno.txt
cat genotypes/*_chr22.fam | awk '{print $1,$2,$6}' | cat header_pheno.txt - > "$sname".pheno #Must have FID, IID, and pheno name in header
cat genotypes/*_chr22.fam | awk '{if($5 == "1") print $1,$2,$6}' | cat header_pheno.txt -   > "$sname"_males.pheno
cat genotypes/*_chr22.fam | awk '{if($5 == "2") print $1,$2,$6}' | cat header_pheno.txt -  > "$sname"_females.pheno


 cat mds_cov/*cogb* mds_cov/*comc*  mds_cov/*meg2* mds_cov/*grac* mds_cov/*pts1* mds_cov/*safr* | awk '{print $1,$2}' > "$sname"_males.exclude
 cat mds_cov/*cogb* mds_cov/*comc*  mds_cov/*meg2* mds_cov/*grac* mds_cov/*pts1* mds_cov/*safr*  mds_cov/*nss1* mds_cov/*ppds*  | awk '{print $1,$2}' > "$sname"_males_nonss.exclude
 cat mds_cov/*cogb* mds_cov/*comc* mds_cov/*grac* mds_cov/pts_meg2_mix_am-qc-af1_pca.menv.mds_cov  mds_cov/*mrsc*   mds_cov/*mrsc* mds_cov/*ppds* mds_cov/*pts1*  mds_cov/*stro*  mds_cov/*wrby* | awk '{print $1,$2}' > "$sname"_females.exclude
 
 rm "$sname".studyindicator
 
 #Make study indicator variable by trolling the mds cov folder
 for files in $(ls mds_cov | grep -v unused)
 do
  fname=$(echo $files | awk -F"_|-" '{print $2_$6}' )
  echo $fname
  awk -v fname=$fname '{print $1,$2,fname}' mds_cov/$files >> "$sname".studyindicator
 done
  awk '{if(NR == 1) $3="studycov"; if(NR == 1 || $1 != "FID") print}' "$sname".studyindicator > "$sname".studyindicator_use
  


 pca_file="$outdir_kinship_grm"/"$sname"_relatedness_pruned_allgwas.vect  #PCA file from last analysis
 exclusions="$sname"_females.exclude #xxxx # pgc_aam.exclude
 phenofile="$sname"_females.pheno #Phenotpye file should be blankspace delimited. It should have columns 1 and 2 be FID and IID. It must have a header. Put phenotype and covariates in here
 phenoname=PTSD #Name of phenotype in pheno file
 covar_names=xxxx #Write covariate names here, if any. Delimit using commas e.g. C1,C2,C3 
 npcs=5 #Number of PCs to include. If none are included, PCA file will just be used to determine unrelated people to include
 studycov="$sname".studyindicator_use #INclude study covariate. Assumes samples in phenotpye file are named by ricopili ID
 phenooutname="$sname"_20rel_females #output filename of pheno and covs
 Rscript scripts/make_ldak_pheno.r $pca_file $exclusions $npcs $phenofile $phenoname $covar_names $studycov $phenooutname
  
  
###Step 3: Make GRM for analysis

#The previous step was an LD pruned matrix used for identifying relateds. This matrix will be based on the whole data and incorporate the LDAK or GCTA model.
#Furthermore this matrix will incorporate weights

#3A: Split genome into sections to weight. Calculate weightings for each section.

#Note: Manual suggests doing this by cohort to prevent genotyping errors from influencing weights. However,this is impractical because of our number of cohorts and their wildly different sizes, so I don't do it!
 
 weightings_keepflag="$outdir_relatednessa"/"$sname"_relatedness_pruned.keep #Put a list of subjects if you want to calculate weights in just a subset of subjects, namely the pruned subject list

 qsub -t1-22 -l walltime=00:30:00 scripts/3_ldak_weights.qsub -d $wd -e errandout/ -o errandout/  -F "-e $ldak -g $outdir_mgeno -n $sname -k $weightings_keepflag -o $outdir_weights -x $lisa"

#Join weightings (Should be very quick,maybe not even needed to run job script)
 for chr in {1..22} ; do $ldak --join-weights  "$outdir_weights"/sections_$chr --bfile "$outdir_mgeno"/"$sname"_"$chr" ; done ; 
 
 #If it is too much of a computational burden, echo the above line to a file and submit a job.. e.g. qsub  -l walltime=00:10:00 -e errandout -o errandout -d $wd weights_join.qsub
 
#3B: Calculate kinship for each chromosome
 
 ldakmodel=gcta #If value set to gcta, no weights are used and power is set to -1. If set to 'unweighted', no weights will be used (power can be set with power). Otherwise, weights will be used
 power=-1 #Set power for LDAK GRM calculation. Speed recommends -0.25. GCTA is -1.

 #40 min for chr 8-22, 1:30 for 1-7, for 12k subjects
 snplist="$wd"/snplist/"$sname"_maf01geno95
 qsub -t1-7 -l walltime=02:00:00 scripts/3B_ldak_kinship_v2_alt_maf.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o $kinship_dir -p $power -w $outdir_weights -x $lisa -z $snplist"
 qsub -t8-22 -l walltime=01:00:00 scripts/3B_ldak_kinship_v2_alt_maf.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o $kinship_dir -p $power -w $outdir_weights -x $lisa -z $snplist"
 qsub -t22 -l walltime=00:05:00 scripts/3B_ldak_kinship_v2_alt_maf.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o $kinship_dir -p $power -w $outdir_weights -x $lisa -z $snplist"

 #Make kinships for other MAF settings

  mkdir "$kinship_dir"/"maf10"
  snplist="$wd"/snplist/"$sname"_maf10geno95
 qsub -t1-7 -l walltime=02:00:00 scripts/3B_ldak_kinship_v2_alt_maf.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o "$kinship_dir"/"maf10" -p $power -w $outdir_weights -x $lisa -z $snplist"
 qsub -t8-22 -l walltime=01:00:00 scripts/3B_ldak_kinship_v2_alt_maf.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o "$kinship_dir"/"maf10" -p $power -w $outdir_weights -x $lisa -z $snplist"

 mkdir "$kinship_dir"/"maf10_a25"
 snplist="$wd"/snplist/"$sname"_maf10geno95
 ldakmodel=unweighted #If value set to gcta, no weights are used and power is set to -1. If set to 'unweighted', no weights will be used (power can be set with power). Otherwise, weights will be used
 power=-0.25 #Set power for LDAK GRM calculation. Speed recommends -0.25. GCTA is -1.
 qsub -t1-7 -l walltime=02:00:00 scripts/3B_ldak_kinship_v2_alt_maf.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o "$kinship_dir"/"maf10_a25" -p $power -w $outdir_weights -x $lisa -z $snplist"
 qsub -t8-22 -l walltime=01:00:00 scripts/3B_ldak_kinship_v2_alt_maf.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -l $ldakmodel -p $power -n $sname -g $outdir_mgeno -o "$kinship_dir"/"maf10_a25" -p $power -w $outdir_weights -x $lisa -z $snplist"

 
#3C: Merge GRMs into entire autosome GRM (inputs: a GRM list and output file name)

#Program parameters:
 biascheck=yes #Set to yes if you want to calculate GRMs for 1/4 segments of the genome, to be used for estimating bias
 qsub -l walltime=00:25:00 scripts/3C_merge_kinships.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -k $kinship_dir -o $combkinship_dir -f kinships -b $biascheck -x $lisa"

 #Remove biascheck, just merge them...
 mkdir "$combkinship_dir"/maf10
 qsub -l walltime=00:25:00 scripts/3C_merge_kinships.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -k $kinship_dir -o $combkinship_dir -f kinships -b no -x $lisa"
 qsub -l walltime=00:25:00 scripts/3C_merge_kinships.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -k "$kinship_dir"/maf10 -o "$combkinship_dir"/maf10 -f kinships -b no -x $lisa"
mkdir "$combkinship_dir"/maf10_a25
  qsub -l walltime=00:25:00 scripts/3C_merge_kinships.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf -k "$kinship_dir"/maf10_a25 -o "$combkinship_dir"/maf10_a25 -f kinships -b no -x $lisa"

  
##4A: Heritablity estimates

#Program parameters:
 phenotype="$sname"_ldak.pheno  #Phenofile file location
 covar="$sname"_ldak.cov #Covariate file location
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned.keep  #If you want to prune extra subjects, beyond what was already done with the relatedness matrix and built into the phenotype file, do so here. E.g. for sex stratified analyses
 outname_append=sens #Extra output names to append

 qsub -l walltime=01:05:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p $phenotype -c $covar -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a $outname_append -b yes -r $prevalence -x $lisa"

#Once this is done we can estimate the bias due to cryptic relatedness:
 grep Her_K "$ldak_resultsdir"/"$outname_append"quad{A,B,C,D,ALL}.reml | awk '(NR<=4){sum+=$2}(NR>4){sum2+=$2}END{I=(sum-sum2)/3;print "Inflation:", I, "=", I/sum2*100,"%"}'

#Now we can do the actual estimate of h2 snp

#Now we can do the actual estimate of h2 snp
 outname_append=defrel_maf01
 prevalence=0.1 #Extra command to set prevalence, so h2 is on the right scale
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned.keep 
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel.pheno -c "$sname"_defrel.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append" -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males.pheno -c "$sname"_defrel_males.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_females.pheno -c "$sname"_defrel_females.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females -b no -r $prevalence -x $lisa"

 
 outname_append=20rel_maf01
 prevalence=0.1 #Extra command to set prevalence, so h2 is on the right scale
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned_allgwas.keep 
 qsub -l walltime=00:05:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_20rel.pheno -c "$sname"_20rel.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append" -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:05:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_20rel_males.pheno -c "$sname"_20rel_males.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:05:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_20rel_females.pheno -c "$sname"_20rel_females.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females -b no -r $prevalence -x $lisa"

 
 outname_append=defrel_maf10
 prevalence=0.1 #Extra command to set prevalence, so h2 is on the right scale
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned.keep 
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel.pheno -c "$sname"_defrel.cov -g "$combkinship_dir"/maf10 -k $keeplist -o $ldak_resultsdir -a "$outname_append" -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males.pheno -c "$sname"_defrel_males.cov -g "$combkinship_dir"/maf10 -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_females.pheno -c "$sname"_defrel_females.cov -g "$combkinship_dir"/maf10 -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females -b no -r $prevalence -x $lisa"

 
 
 outname_append=20rel_maf10
 prevalence=0.1 #Extra command to set prevalence, so h2 is on the right scale
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned_allgwas.keep 
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_20rel.pheno -c "$sname"_20rel.cov -g "$combkinship_dir"/maf10 -k $keeplist -o $ldak_resultsdir -a "$outname_append" -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_20rel_males.pheno -c "$sname"_20rel_males.cov -g "$combkinship_dir"/maf10 -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_20rel_females.pheno -c "$sname"_20rel_females.cov -g "$combkinship_dir"/maf10 -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females -b no -r $prevalence -x $lisa"

  
 outname_append=defrel_maf10_a25
 prevalence=0.1 #Extra command to set prevalence, so h2 is on the right scale
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned.keep 
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel.pheno -c "$sname"_defrel.cov -g "$combkinship_dir"/maf10_a25 -k $keeplist -o $ldak_resultsdir -a "$outname_append" -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males.pheno -c "$sname"_defrel_males.cov -g "$combkinship_dir"/maf10_a25 -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_females.pheno -c "$sname"_defrel_females.cov -g "$combkinship_dir"/maf10_a25 -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females -b no -r $prevalence -x $lisa"

 

 outname_append=defrel
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned.keep 
 qsub -l walltime=00:25:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel.pheno -c "$sname"_defrel.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a $outname_append -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males.pheno -c "$sname"_defrel_males.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_females.pheno -c "$sname"_defrel_females.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females -b no -r $prevalence -x $lisa"

 qsub -l walltime=00:10:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males_nonss.pheno -c "$sname"_defrel_males_nonss.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males_nonss -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:10:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males_nonmiss.pheno -c "$sname"_defrel_males_nonmiss.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males_nonmiss -b no -r $prevalence -x $lisa -j xxxx"
 
  qsub -l walltime=00:10:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males.pheno -c "$sname"_defrel_males.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males_constrained -b no -r $prevalence -x $lisa -j YES"
 
 #I put teh new maf10 grm here for htis
  qsub -l walltime=00:10:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males.pheno -c "$sname"_defrel_males.cov -g "$combkinship_dir"/maf10 -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males_maf10 -b no -r $prevalence -x $lisa -j xxxx"
 #and high genotyping : name may need fixing
  qsub -l walltime=00:10:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_males.pheno -c "$sname"_defrel_males.cov -g "$combkinship_dir"/maf10genohigh -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males_maf10 -b no -r $prevalence -x $lisa -j xxxx"
   qsub -l walltime=00:10:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_defrel_females.pheno -c "$sname"_defrel_females.cov -g "$combkinship_dir"/maf10genohigh -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females_maf10geno -b no -r $prevalence -x $lisa -j xxxx"
 
 
 qsub -l walltime=00:40:00  plink_assoc_test.sh -e errandout/ -o errandout/ 
 
 outname_append=analysis2_allgwasp2
 keeplist="$outdir_kinship_grm"/"$sname"_relatedness_pruned_allgwasp2.keep 
 qsub -l walltime=00:25:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_lowrel.pheno -c "$sname"_lowrel.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a $outname_append -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_lowrel_males.pheno -c "$sname"_lowrel_males.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_males -b no -r $prevalence -x $lisa"
 qsub -l walltime=00:15:00 scripts/4_h2.qsub -d $wd -e errandout/ -o errandout/ -F "-e $ldakf  -p "$sname"_lowrel_females.pheno -c "$sname"_lowrel_females.cov -g $combkinship_dir -k $keeplist -o $ldak_resultsdir -a "$outname_append"_females -b no -r $prevalence -x $lisa"

  for pop_prev in 0.1 0.5
  do
   for results1 in 20rel defrel 
   do
    for results2 in maf01 maf10
    do
    #Overall
    res=$(grep Her_K "$ldak_resultsdir"/"$results1"_"$results2"full.reml | awk '{print $2,$3}')
    h2=$(echo $res | awk '{print $1}')
    se=$(echo $res | awk '{print $2}')
    lrt=$(grep LRT_P "$ldak_resultsdir"/"$results1"_"$results2"full.reml | awk '{print $2}')
    ncases=$(awk '{print $3}' "$sname"_"$results1".pheno | grep 2 | wc -l | awk '{print $1}')
    ncontrols=$(awk '{print $3}' "$sname"_"$results1".pheno | grep 1 | wc -l | awk '{print $1}')
        
    rj=$(Rscript scripts/liability_scaleh2.r $h2 $se $pop_prev xxxx "$sname"_"$results1".pheno)
    echo "Overall "$results1"_"$results2"_"$pop_prev" $ncases $ncontrols $rj $lrt "
    #Males
    res=$(grep Her_K "$ldak_resultsdir"/"$results1"_"$results2"_malesfull.reml | awk '{print $2,$3}')
    h2=$(echo $res | awk '{print $1}')
    se=$(echo $res | awk '{print $2}')
    lrt=$(grep LRT_P "$ldak_resultsdir"/"$results1"_"$results2"_malesfull.reml | awk '{print $2}')
    ncases=$(awk '{print $3}' "$sname"_"$results1"_males.pheno | grep 2 | wc -l | awk '{print $1}')
    ncontrols=$(awk '{print $3}' "$sname"_"$results1"_males.pheno | grep 1 | wc -l | awk '{print $1}')
     
    rj=$(Rscript scripts/liability_scaleh2.r $h2 $se $pop_prev xxxx "$sname"_"$results1"_males.pheno)
   
    echo "Males "$results1"_"$results2"_"$pop_prev" $ncases $ncontrols $rj $lrt "
    #Females
    res=$(grep Her_K "$ldak_resultsdir"/"$results1"_"$results2"_femalesfull.reml | awk '{print $2,$3}')
    h2=$(echo $res | awk '{print $1}')
    se=$(echo $res | awk '{print $2}')
        lrt=$(grep LRT_P "$ldak_resultsdir"/"$results1"_"$results2"_femalesfull.reml | awk '{print $2}')
    ncases=$(awk '{print $3}' "$sname"_"$results1"_females.pheno | grep 2 | wc -l | awk '{print $1}')
    ncontrols=$(awk '{print $3}' "$sname"_"$results1"_females.pheno | grep 1 | wc -l | awk '{print $1}')
     
    rj=$(Rscript scripts/liability_scaleh2.r $h2 $se $pop_prev xxxx "$sname"_"$results1"_females.pheno)
    echo "Females "$results1"_"$results2"_"$pop_prev" $ncases $ncontrols $rj $lrt "
    done
  done
done
 

 
