# ldak_pipeline
Operation: Make a working directory. Inside of this folder, put all scripts in a folder called 'scripts'. Follow instructions in 1_get_h2_v2.sh to perform LDAK analysis.  


to extract data to run an h2 analysis,  
tar xvf /archive/maihofer/freeze2gwas_ldakeur.tar --wildcards "ldakeur/scripts/*" --wildcards "ldakeur/combined_kinships/*" --wildcards "*.pheno" --wildcards "*.cov" --wildcards "*.sh" --wildcards "*.vect"
