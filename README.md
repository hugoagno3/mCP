# mCP

![imagen](https://user-images.githubusercontent.com/82643524/179713559-a0f27ad1-07af-4d69-ac72-380d4eba304b.png)

mini-Complexome Profiler


---
title: "mCP_tutorial_Humans_samples"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{mCP_tutorial_Humans_samples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(mCP)
```
# Overview 

The `mCP` package provides a set of functions for targeted protein complex detection in mass spectrometry-based proteomics co-fractionation experiments. It is designed to work with experimental data in the form of protein abundance matrices, and to compare the data to the CORUM database of known protein complexes to identify complexes that are present in the experimental data.

This vignette provides a step-by-step guide to using the `mCP` package to analyze proteomics data and identify protein complexes of interest.
We provide a processing with "hsapiens" dataset Hek293 processed by BNE page of 35 fractions. We explain in the last section how to run an extra function that will be very useful for getting the input matrix as the average of the replicates in the experiment. All functions have different options, therefore we provide a detailed explanation about the functions at the end of this tutorial.

# Installation

To install the mCP package from GitHub, run the following code:
```{r pressure, echo=FALSE}
library(devtools)
devtools::install_github("hugoagno3/mCP")
```


# Data input (two files Corum database and Experimental wide format results)
   mCP R package is focused on detection of protein complexes and it accepts only protein data as an input matrix. The mass spectrometry data aquisition can be done by Data Dependent Acquisition mode or Data independent Aquisition.  The following steps are recommended for pre-processing the data before getting into mCP:

1. Normalize the data between fractions in experiments.  
2. Replace missing values by 0.
2. For comparison the normalization of the data to plot complexome profiling plots can be done inside each protein complexe by activateion the relativization fuction relative= TRUE in ccp_plotter function.
3. Normalize per Biologica replicates if fractionation is reproducible. 


To perform targeted protein complexes detection, two input files are needed:
           1- A CORUM database of protein complexes (or a targeted list of interest) dataframe with 3 columns:´  col1= "complex_id", col2= "complex_name" and col3= "protein_id").
           2- An experiment file *data.frame* with your experiment results in wide-format- first column called "protein_id" and the next columns contains protein intensity and belong to the fractions in the co-fractionation experiment: can be numerics names from 1 to the last number of fractions, for example 1,2,3,..,35 for 35 fractions, means 35 columns.
           In this example, we will use the Corum_Humans_Database file as the protein complex database and the Hek293_P2_1.csv and Hek293_P2_2.csv files as the experiment files.

```{r pressure, echo=FALSE}
library(readr)
library(dplyr)

# Read the Corum protein complex database file
      data(Corum_Humans_Database)

# Read the experiment files
      data(Hek293_P2_1)
```

# Data Processing (Input of matrix with long names wide format)

To process the input data, we need to run *Option 1*- the mCP function or *Option 2* 3 functions provided by our mCP package.
# Option 1: Running mCP function
   This function is an integrated function of mCP package, that needs as input an experimental data and returns a list of plots, binary total hits, id of proteins of binary hits and heatmaps_seaborn of know protein complexes detected in CORUM database. In addition, it plots  6 files as  outputs: 
 1- pdf file with All candidates detected as protein complexes profiles from Corum database.
 2- pdf with with all candidates heatmaps of the detected protein complexes.
 3- pdf file with the detected as protein complexes profiles with a controlled FDR.
 4- pdf with heatmaps of the detected protein complexes with a controlled FDR.
 5- txt file with numbers about general false positive when at least 1 hit is consider as filter.
 6- CVS file containing all protein complexes detected, hits of binary interactions inside the protein complexes, FDR detected by MonteCarloSimulation.
   An example can be found here:
 
```{r pressure, echo=FALSE}
library(mCP)
# Read the Corum protein complex database file
      data(Corum_Humans_Database)
# Read the experiment files
      data(Hek293_P2_1)
##### Example #####
out_Hek_P2_1_teste <- mCP(corum_database = Corum_Humans_Database,
                          experiment_data = Hek293_P2_1, 
                          N_fractions = 35, 
                          specie = "hsapiens",
                          method_cor = "pearson",
                          heatmap_seaborn = TRUE,
                          format = "pdf", 
                          output_name = "m_CP_analysis_2",
                          filter = 0.93,
                          heat_map = TRUE,
                          relative = FALSE,
                          n_simulations= 9,
                          display_weights = TRUE,
                          standard_weights = list(list(x =11, label= "1049KDa"), 
                                             list(x = 13, label ="720 KDa")))
```
 
# Option 2: Runn mCP function individually

```{r pressure, echo=FALSE}
library(dplyr)
library(mCP)
# Read the Corum protein complex database file
      data(Corum_Humans_Database)
# Read the experiment files
      data(Hek293_P2_1)
#### RUN mCP list that creates a list of potential protein complexes
CL_hek_P2_1<- mcp_list(corum_database =  Corum_Humans_Database,
                       experiment_data = Hek293_P2_1, 
                       N_fractions = 35, 
                       specie = "hsapiens",
                       method_cor = "pearson",
                       heatmap_seaborn = TRUE)
##### Run the output of mcp_list into the cpp_ploter function. 
out_Hek_P2_1 <- cpp_plotter(complex_list = CL_hek_P2_1,
                            format = "pdf", 
                            output_name = "m_CP_analysis",
                            filter = 0.93,
                            N_fractions = 35,
                            heat_map = TRUE,
                            relative = FALSE,
                            display_weights = TRUE,
                            standard_weights = list(list(x =11, label= "1049KDa"), 
                                                    list(x = 13, label ="720 KDa")))
###### Run the Experimental matrix and the output of cpp_ploter into the fdr_mcp function. Keep the same filter value!
FDR_DIANN_dDIA_Hek_P2_1_<- fdr_mCP(corum_database= Corum_Humans_Database,
                               Output_cpp_plotter = out_Hek_P2_1, 
                              experiment_data=Hek293_P2_1,
                              file_name = "m_CP_analysis",
                              N_fractions = 35,
                              specie = "hsapiens",
                              fdr_limit = 0.05,
                              filter=0.93,
                              n_simulations= 185)
##Note this function will perform a FDR filter and simultation it is important to do the simulation with the same filter then the previous function. if filter value is different the simulation do not makes sense. 
```




Protein/protein complexes (PPC) are key components of cellular and physiological processes. Several cardiac diseases, such as Atrial Fibrillation1, are associated with disruptions of these complexes. Due to their role in health and disease, PPCs have been the focus of much research2–4. Here, we introduce a new workflow named mini-Complexome Profiler (m-CP) that allows sensitive detection of PPCs from HEK293 and mouse primary myocytes. Unlike previously established approaches such as CCProfiler2, m-CP is suitable for PPC identification in low fraction number data, making it a useful tool for screening or analyzing difficult-to-obtain samples. We were able to reproduce previous observations showing that phospholamban (PLN) is part of a cardiac Ca+2 handling complex, and that its deletion results in alteration of this structure in isolated mouse ventricular cardiomyocytes knockout conditions (KO-PLN-/-)1. 
![imagen](https://user-images.githubusercontent.com/82643524/179712587-c01491a4-9a5c-49d1-8d2d-e90077970e21.png)
m-CP is a workflow for identification and comparison of PPC stablished for HEK and ventricular cardiomyocytes cells
m-CP is suitable for rapid screening of PPC to detect variations between two or more different conditions as shown in figure 5
It is ideal to study difficult-to-obtain or scarce samples. We showed increased sensitivity in low-fraction-number samples compared to a previous state of art approach2
The workflow requires 34 FR, 35-71µg protein per sample and 14 days lab & bioinformatic work
![imagen](https://user-images.githubusercontent.com/82643524/179712766-6099285a-7512-4b1d-a6a3-d9eb7df89fa3.png)
Bibliography
1.PMID: 31185731
2.PMID: 30642884
3. PMID: 22939629
4. PMID: 30357367
5.PMID: 17010801 


