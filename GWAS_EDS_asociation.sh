# GWAS analysis

plink2 --bfile genome_3523 --pheno phen.txt --pheno-name EDS204,EDS81 --covar phen.txt --covar-name age,sex,SNPPC1,SNPPC2,SNPPC3,SNPPC4,SNPPC5 --covar-variance-standardize --glm hide-covar --maf 0.01 --out GWAS_EDS
