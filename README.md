# EpigeneticDriftScore
DNA methylation; Epigenetic Drift Score; Aging

![Build](https://github.com/fan-7/EpigeneticDriftScore)

# Introduction

Epigenetic drift refers to the accumulation of stochastic changes in DNA methylation over time that is particularly noticeable at specific cytosine-phosphate-guanine (CpG) dinucleotides. Unlike epigenetic clock sites, where changes in DNA methylation linearly associate with age, epigenetic drift sites are characterized by an increase in variance in the methylated positions with age. Biologically, epigenetic 
drift may reflect the heightened heterogeneity in health status among the elderly population. This epigenetic drift stochasticity may set the pace of aging and provide opportunities for therapies that decelerate the aging process. Conventional methods for identifying drifting CpGs have shown limitations in handling non-linear drifts and lack a robust quantitative assessment score. In response, we introduce a robust white heteroscedasticity method for detecting drift-CpGs. We executed an Epigenome-Wide Drift Study (EWDS) in a discovery cohort, followed by replication studies in independent cross-ethnicity cohorts.  Cross-ethnicity analysis confirmed the universal nature of epigenetic drift, with 99% of drift-CpGs showing increased variance and a small subset exhibiting decreased variance with age. We developed two epigenetic drift scores based on 204 positive and 81 negative drift-CpGs, providing robust measures of epigenome-wide drift. The positive drift score was strongly associated with age, other epigenetic measures, and lipoprotein metabolites. GWAS identified 9q33.1 and 2p21 as susceptibility loci for the positive and negative drift scores, respectively. These findings enhance our understanding of epigenetic drift's role in aging.


This repo contains code needed to generate figures for the paper 

    A robust epigenetic drift score reveals genetic and metabolic insights into aging

The code is under development

## The analysis is performed in 6 steps:

1. Benchmarking of Epigenetic Drift Statistical Method
2. Development of a Robust Drift-CpG List
3. Creation of a Genome-Wide Epigenetic Drift Score
4. Association with Epigenetic Clocks
4. Association with Lipid Metabolites and Biomarkers
5. Genetic Basis

R Markdown Notebooks Used for Analysis and Figure Generation



## Step 1. 

Each CpG site is processed separately with heteroscedasticity test:

    source("FUN_EpigeneticDrift_WhiteP_robustcoef.r")
    source("FUN_EpigeneticDrift_BP_rawcoef_test.r")
    source("FUN_EpigeneticDrift_BP_LMR_test.r")
    source("FUN_EpigeneticDrift_DGLM_bacon_test.r")
    getWhitePdriftCpG(method="white", alpha, nSamples, nIterations, age, cpgdata, covdata)
    getBPdriftCpG_rawcoef(method="BP", cpg, nIterations, alpha)
    getBPdriftCpG_fit_age_dispersion(y,fm, mf)
    getBPdriftCpG_runDGLM(x,cpgnames)


## Step 2.

Each individual is predicted separately with EDS:

    source("FUN_EDS_POS.r")
    source("FUN_EDS_NEG.r")
    input <- readRDS("input.rds")
    eds_pos <- EDS_POS(input, ref_stat, ref_drift, ref_coef)
    eds_neg <- EDS_NEG(input, ref_stat, ref_drift, ref_coef)
    

cd into folder where you want to do the analyis

mkdir EpigeneticDriftScore
cd EpigeneticDriftScore
Clone the git repo

git clone https://github.com/fan-7/EpigeneticDriftScore.git
Create conda environment with all necessary tools installed by:

conda env create -f EpigeneticDrift/envs/evmake.yaml
Install extra R libraries that are not present in conda or have trouble working:

TODO 


## More

More info about Epigenetic drift can be found on our [website](https://github.com/fan-7/EpigeneticDriftScore.git ).  
