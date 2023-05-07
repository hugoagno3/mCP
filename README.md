# mCP

![imagen](https://user-images.githubusercontent.com/82643524/179713559-a0f27ad1-07af-4d69-ac72-380d4eba304b.png)

mini-Complexome Profiler


 "mCP_tutorial_Humans_samples"

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
## Last step to get the protein complexes filtered by FDR is to get back to the list and run cpp_plotter again.
Note that if you wish to have the plots of this last FDR protein complexes  you could filter the names of the protein complexes detected into the first list (from mcp_list) and run cpp_ploter again. 


```{r}
#### list to plot
   CL_final_output<- CL_hek_P2_1[names(out_Hek_P2_1)]
#### Then run again this list on the cpp_plotter function as follows. 
   out_Hek_P2_1_final_output <- cpp_plotter(complex_list = CL_final_output,
                            format = "pdf", 
                            output_name = "m_CP_analysis",
                            filter = 0.93,
                            N_fractions = 35,
                            heat_map = TRUE,
                            relative = FALSE,
                            display_weights = TRUE,
                            standard_weights = list(list(x =11, label= "1049KDa"), 
                                                    list(x = 13, label ="720 KDa")))

```

### PLots
The output list out_Hek_P2_1_final_output has a list of protein complexes each protein complex is an element of that list with 4 objects
1) The first object is a profile of proteins, fractions vs intensities. 




![13S condensin Complex](https://user-images.githubusercontent.com/82643524/236702562-3e656a6d-3272-422f-9419-2c698a1d0a09.png)




```{r my-plot, fig.cap = "My Plot Caption"}
my_plot<- out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[2]]
  
png("my_plot.png")
 ![13S condensin Complex](https://user-images.githubusercontent.com/82643524/236702514-c9356bfe-dd80-4c9d-8c2a-5504c10c7ba8.png)

  plot(my_plot)
 dev.off()
 knitr::include_graphics("my_plot.png")
```

2) The second object is the number of Hits detected in the experiment as the number of significant binary interactions. So binare pearson correlations values in the experiment higher than a filter.  

```{r}
 out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[2]]

[[2]]
  Hits
1   10
```


## The output list out_Hek_P2_1_final_output has 4 objects
3) The name of the proteins involved in the Binary interaction detected (Experimental Protein-binary Pearson correlatio nhigher than the filter) and its Pearson correlation of thes significant binary interactions.

```{r}
out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[3]]

[[3]]
[[3]]$interactions
 [1] "NCAPH:NCAPD2" "SMC4:NCAPD2"  "NCAPH:NCAPG"  "SMC4:NCAPG"   "NCAPG:NCAPH"  "SMC4:NCAPH"  
 [7] "NCAPG:SMC2"   "SMC4:SMC2"    "NCAPG:SMC4"   "SMC2:SMC4"   

[[3]]$Pearson
 [1] 0.9939944 0.9552057 0.9912144 0.9463247 0.9912144 0.9585936 0.9776026 0.9662482 0.9463247 0.9662482

```

4) Network heatmap constructed from the correlation matrix. 

```{r my-plot, fig.cap = "My Plot Caption"}
my_plot<- out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[4]]
  
png("my_plot.png")
 
  plot(my_plot)
 dev.off()
 knitr::include_graphics("my_plot.png")
```


# Data Processing (Input of matrix with protein_id as first column and Factions)
  Most of the time we measure samples per duplicate in co-fractionation experiemnt, mCP provide a function to deal with these extensive number of columns. The function input can be a data.frame imported by the following function read.table
```{r}
NAmatrix_P2_1 <- read.table("Hek293_P2_1.csv",sep =",",dec = ".", header= T)
```
  once imported you will get this file, we provide it in the R package, just write data(NAmatrix_P2_1)
```{r}

```
##  After that you can run the funtion calc_mean_matrix
  The function calc mean matrix will recognize the replicates by a tag like _A_ and _B_ as part of the names for the MS file and the *Frac_index*. It is the position in which you can find the number of fraction in the column name separated by underscores. for example the frac_index=5 is date_Surname_Measurement_A_01 because the number for the fraction is 5 spaces between underscores, like a_b_c_d_FractionNumber. 
```{r}
  ###  Note: Here the frac_index= is 17. and the pattern_group= _A_ and pattern_group= _B_  
     names(Hek293_P2_1[,2])
    [1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_A_01_GA1_1_6588.d"
    names(Hek293_P2_1[,37])
    [1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_B_01_GA1_1_6626.d"

```
  
  So we can run the function as follow. 

```{r}

Hek_293_MA_P2_1 <- Calc_mean_matrix(NAmatrix = NAmatrix_P2_1,
                            pattern_group_A = "_A_",
                            pattern_group_B = "_B_",
                            frac_index=17,
                            Protein_ID_column = 1,
                            save_file = TRUE,
                            save_name = "Hek293_P2_1.csv")
```

Note: The user can also use its own scrip to generate the average of the replicates. This function is limited to 2 replicates. It is our traditional measurement. mCP do not acepts NA values. All NA values in our case are replaced by 0. 

# How does mCP works deteil fo the functions and explanation step by step

The `mCP` package provides the following functions:

- `Calc_mean_matrix()`: pre-processes the data before running `mcp_list()`.

- `mcp_list()`: Extracts Potential protein complexes in the experiment from CORUM database.
    It extract intensities of all protein complexes from the CORUM database and generates a list of all potential protein complexes presents. Each element of the list is identified as a detected potential protein complex and named with its corresponding complex name. Each protein complex element contains 2 elements: 

 *my_protein_complex_mcp_list[1][1]* is a dataframe composed of the Id proteins of the proteins, names, and fractions and its correspondig intensities per fraction. Names are annotated to the uniprot accesion Ids by a repository R package called gprofiler2. 

 *my_protein_complex_mcp_list[1][2]* is a correlation matrix of these elements,  this correlation matrix is done by Pearson correlations (but kendall and spearman algoriths are also posible in that step). and calculates a matrix correlation between all detected component of each present protein complex.

- `cpp_plotter()`: The cpp_plotter() function takes a list of potential protein complexes as input. The function first filters out any complexes composed of less than 2 components with non-zero intensities. Then, it applies a filter based on the minimum definition of a protein complex, keeping only complexes with at least one significant co-elution hit. The function also filters the list of complexes based on the presence of co-eluting proteins, requiring at least one binary interaction higher than the filter for proteins within the candidate complex. Finally, cpp_plotter() generates complexome profiling plots, heatmaps, and network plots of the proteins within the selected complexes.

The list of potential protein complexes is the input for the next function called cpp_plotter. This is the first filter for  Then it filters protein complexes composed by at least 2 component with intensities different from cero. Then the first filter is runned: The filter is based in a minimun definition of protien complex. So in this step it keep only protein complexes with at least one significant coelution hit.
filters the list of protein complexes based on the presence of co-eluting proteins calculates as a minimun of 1 binary interactions higher that the filter of protien within the candidate protein complex and creates complexome profiling plots, heatmaps and networks plots of proteins with in protein complexes of the selected complexes.

*my_protein_complex_mpp_ploter[1][1]* The first object is a profile of proteins, fractions vs intensities.
*my_protein_complex_mpp_ploter[1][2]* The second object is the number of Hits detected in the experiment as the number of significant binary interactions. So binare pearson correlations values in the experiment higher than a filter.  
*my_protein_complex_mpp_ploter[1][3]* The name of the proteins involved in the Binary interaction detected.
*my_protein_complex_mpp_ploter[1][3]* Pearson correlation detected higher than the filter of the detected binaries interactions. 
*my_protein_complex_mpp_ploter[1][4]* Network heatmap constructed from the correlation matrix.


