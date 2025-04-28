# EpigeneticDriftScore(EDS)
DNA methylation; Epigenetic Drift Score; Aging

![Schematic representation]([https://github.com/fan-7/EpigeneticDriftScore](https://github.com/fan-7/EpigeneticDriftScore/blob/main/Figures/Schematic%20representation.png))

# Introduction

Epigenetic drift refers to the progressive, stochastic accumulation of molecular alterations across the epigenome during aging. These include changes in DNA methylation, histone modifications, chromatin remodeling, and non-coding RNAs. Collectively, such alterations disrupt gene regulatory networks, leading to transcriptional dysregulation, loss of cellular homeostasis, and increased vulnerability to age-related diseases. Among these, DNA methylation drift has emerged as the most extensively characterized component of epigenetic aging. It is typified by increased variance in methylation levels at specific CpG sites-termed drift-CpGs-over chronological age, distinguishing it from the linear age-associated methylation changes at epigenetic clock sites. Importantly, epigenetic drift likely unfolds more gradually across wider temporal windows, capturing interindividual and cellular heterogeneity in a way that reflects biological aging beyond the linear progression captured by clock-based models. This stochastic accumulation is hypothesized to reflect rising intercellular and interindividual heterogeneity with age, potentially capturing biological aging more dynamically than static methylation clocks (Meyer and Schumacher 2024). Understanding and quantifying epigenetic drift may offer novel biomarkers for aging trajectories, disease susceptibility, and therapeutic interventions aimed at mitigating age-associated decline.

We first performed a comprehensive evaluation of four commonly used statistical methods for identifying drift-CpGs, using both computer simulations and empirical population data. We then applied this method to 735,267 CpG sites measured in 3,538 Chinese individuals, systematically identifying a high-confidence set of drift-CpGs. Our replication analysis was conducted in 2,423 individuals from two additional Chinese cohorts and two European cohorts. Both positive and negative drift-CpGs were characterized to reveal complementary biological patterns. We further investigated whether drift-CpGs exhibit cell-type specificity or represent shared signatures of hematopoietic aging by integrating single cell RNA-seq analyses. To quantify individual-level epigenetic drift, we developed an EDS based on a defined subset of drift-CpGs, offering a standardized alternative to entropy-based metrics. Finally, we evaluated the functional relevance of EDS through association analyses with age, lipidomic profiles, and genome-wide genetic variation, uncovering potential mechanistic links between stochastic epigenetic changes and age-related complex traits. 

This repo contains code needed to generate figures for the paper 

    A robust epigenetic drift score reveals genetic and metabolic insights into aging

The code is under development

## The analysis consists of 7 steps:
1. Benchmark the epigenetic drift statistical method to ensure accuracy.
2. Develop a robust drift-CpG list for epigenetic drift analysis.
3. Analyze the cellular heterogeneity and transcriptional associations of drift.
4. Create a genome-wide epigenetic drift score (positive/negative).
5. Explore the association with epigenetic clocks related to aging.
6. Examine the links with lipid metabolites and biomarkers.
7. Investigate the genetic basis influencing epigenetic drift.

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

Install the R package EDS to calculate epigenetic drift score capturing directional methylation variability:

    download EDS_0.1.0.tar.gz
    remotes::install_local("EDS_0.1.0.tar.gz",upgrade = F,dependencies = T)
    library(EDS)
    input <- read.csv("input.csv", header = T)
    eds_pos <- EDS_POS(input, ref_stat, ref_drift, ref_coef)
    eds_neg <- EDS_NEG(input, ref_stat, ref_drift, ref_coef)
    

## Step 3.

Or run the customized script manually:

    source("FUN_EDS_POS.r")
    source("FUN_EDS_NEG.r")
    input <- read.csv("input.csv", header = T)
    eds_pos <- EDS_POS(input, ref_stat, ref_drift, ref_coef)
    eds_neg <- EDS_NEG(input, ref_stat, ref_drift, ref_coef)
    

cd into folder where you want to do the analyis

mkdir EpigeneticDriftScore
cd EpigeneticDriftScore
Clone the git repo
git clone https://github.com/fan-7/EpigeneticDriftScore.git


TODO 


## More

More info about Epigenetic drift can be found on our [website](https://github.com/fan-7/EpigeneticDriftScore.git ).  
