# mCP

![imagen](https://user-images.githubusercontent.com/82643524/179713559-a0f27ad1-07af-4d69-ac72-380d4eba304b.png)

mini-Complexome Profiler (mCP) 


 "mCP_tutorial_Humans_samples"

# Overview 

The `mCP` package provides a set of functions for targeted protein complex detection in mass spectrometry-based proteomics co-fractionation experiments. It is designed to work with experimental data in the form of protein abundance matrices, and to compare the data to the CORUM database of known protein complexes to identify complexes that are present in the experimental data. The package contains two datasets as example: one example for "hsapiens" and on efor "mmusculus" Co fractionation experiment. 

This vignette provides a step-by-step guide to using the `mCP` package to analyze proteomics data and identify protein complexes of interest.
We provide an Experimental dataset Hek293 called "hsapiens" fractionated by BNE-PAGE of about 35 fractions.
In the last section,the user can find how to run an extra function called `Calc_mean_matrix()` that will be very useful for getting the input matrix as the average of the replicates in the experiment. All functions have different options, therefore 
we provide a detailed explanation about the functions at the end of this tutorial.



![heat eIF3 complex](https://user-images.githubusercontent.com/82643524/236703288-735a5d6e-25aa-49b9-8520-071a3df7ea4c.png)
![Eif3j heatmap](https://user-images.githubusercontent.com/82643524/236703265-8cd3be1e-4252-4878-8f91-58481916cb84.png)

# Installation

To install the mCP package from GitHub, run the following code:
```{r pressure, echo=FALSE}
library(devtools)
devtools::install_github("hugoagno3/mCP", biuld_vignette= TRUE)
```

# Data input 
 In this tutorial, we will use the Corum_Humans_Database file as the protein complex database and the Hek293_P2_1 file as experiment files. To perform targeted protein complexes detection, only these two input files are needed:
```{r pressure, echo=FALSE}
# Open the Corum protein complex database file
      data(Corum_Humans_Database)
# open the experiment files
      data(Hek293_P2_1)
```
* A CORUM database of protein complexes (or a targeted list of interest) dataframe with 3 columns:´  col1= "complex_id", col2= "complex_name" and col3= "protein_id").
* An experiment file *data.frame* with your experiment results in wide-format- first column called "protein_id" with a single uniprot accession per row as elements and the next columns contain protein intensity detected in the fractions in the co-fractionation experiment: can be numerics names from 1 to the last number of fractions, for example 1,2,3,..,35 for 35 fractions, means 35 columns.

Note: mCP R package is focused on detection of protein complexes and it accepts only protein data Level as an input matrix. The mass spectrometry data aquisition can be done by Data Dependent Acquisition mode or Data independent Aquisition.  

# Data Processing

To process the input data, we need to run *Option 1*- the mCP function or *Option 2*- 3 functions provided by our mCP package.
# Option 1: Running mCP function
   This function is an integrated function of mCP package, that needs as input an experimental data and returns a list of plots, binary total hits, id of proteins of binary hits and heatmaps_seaborn of know protein complexes detected in CORUM database. In addition, it plots  6 files as  outputs: 
 1- pdf file with All candidates detected as protein complexes profiles from Corum database.
 2- pdf with with all candidates heatmaps of the detected protein complexes.
 3- pdf file with the detected as protein complexes profiles with a controlled FDR.
 4- pdf with heatmaps of the detected protein complexes with a controlled FDR.
 5- txt file with numbers about general false positive when at least 1 hit is consider as filter.
 6- CVS file containing all protein complexes detected, hits of binary interactions inside the protein complexes, FDR detected by MonteCarloSimulation.
   An example can be found here:
 
```{r pressure3, eval=FALSE, message=FALSE, include=TRUE}
library(mCP)
```

```{r pressure3A, eval=FALSE, include=TRUE}
# Read the Corum protein complex database file
      data(Corum_Humans_Database)
```

```{r pressure3B, eval=FALSE, include=TRUE}
# Read the experiment files
      data(Hek293_P2_1)
```

```{r pressure3C, eval=FALSE, include=TRUE}
      ##### Example #####
mCP_Hek_P2_1 <- mCP(corum_database = Corum_Humans_Database,
                          experiment_data = Hek293_P2_1, 
                          N_fractions = 35, 
                          specie = "hsapiens",
                          method_cor = "pearson",
                          heatmap_seaborn = TRUE,
                          format = "pdf", 
                          output_name = "m_CP_analysis",
                          filter = 0.81,
                          heat_map = TRUE,
                          relative = FALSE,
                          fdr_limit = 0.05,
                          n_simulations= 9,
                          monomeric_filter = FALSE,
                          Risk_fraction = 31,
                          set_seed = TRUE,
                          display_weights = TRUE,
                          standard_weights = list(list(x =6, label= "2700 KDa"), 
                                                          list(x = 11, label ="950 KDa"),
                                                          list(x = 14, label = "750 KDa"), 
                                                          list(x =27, label ="146 KDa"),
                                                          list(x =30, label ="60 KDa")))
```
 # Outputs plots per Protein complex detected:

 mCP_TEND_out_Hek_P2_1_teste_11$`Respiratory chain complex I (intermediate VII/650kD), mitochondrial`

 [[1]]
 ![respiratory mitochondrial](https://user-images.githubusercontent.com/82643524/236703671-d1cdcda6-73c0-4ede-8258-679c28fe234e.png)
[[2]]
  Hits
1    9

[[3]]
[[3]]$interactions
[1] "NDUFS5:NDUFA6" "NDUFA6:NDUFA9" "NDUFS3:NDUFA9" "NDUFA9:NDUFS2" "NDUFS7:NDUFS2" "NDUFS2:NDUFS3"
[7] "NDUFA6:NDUFS5" "NDUFA6:NDUFV2" "NDUFS3:NDUFV2"

[[3]]$Pearson
[1] 0.9434851 0.9346499 0.9487094 0.9306028 0.9495674 0.9727912 0.9434851 0.9427735 0.9544865
 
 [[4]]
 ![seahet mitochondrial respiratory](https://user-images.githubusercontent.com/82643524/236703676-5ac98435-02b4-48d8-a1b6-e4ceeeae5737.png)

# Option 2: Runn mCP function individually
The following code guide you to get the output step by step

```{r pressure4, eval=FALSE, include=TRUE}
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
                            filter = 0.81,
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
                              filter=0.81,
                              Risk_fraction = 31,
                              monomeric_filter = FALSE,
                              set_seed = TRUE,
                              n_simulations= 185)
##Note this function will perform a FDR filter and simultation it is important to do the simulation with the same filter then the previous function. if filter value is different the simulation do not makes sense. 
```
## Last step to get the protein complexes filtered by FDR is to get back to the list and run cpp_plotter again.
Note that if you wish to have the plots of this last FDR protein complexes  you could filter the names of the protein complexes detected into the first list (from mcp_list) and run cpp_ploter again. 


```{r5}
#### list to plot
   CL_final_output<- CL_hek_P2_1[names(FDR_DIANN_dDIA_Hek_P2_1_)]
#### Then run again this list on the cpp_plotter function as follows. 
   out_Hek_P2_1_final_output <- cpp_plotter(complex_list = CL_final_output,
                            format = "pdf", 
                            output_name = "m_CP_analysis",
                            filter = 0.81,
                            N_fractions = 35,
                            heat_map = TRUE,
                            relative = FALSE,
                            display_weights = TRUE,
                            standard_weights = list(list(x =6, label= "2700 KDa"), 
                                                          list(x = 11, label ="950 KDa"),
                                                          list(x = 14, label = "750 KDa"), 
                                                          list(x =27, label ="146 KDa"),
                                                          list(x =30, label ="60 KDa")))
```

### PLots
The output list out_Hek_P2_1_final_output has a list of protein complexes each protein complex is an element of that list with 4 objects
1) The first object is a profile of proteins, fractions vs intensities. 
```{r6 my-plot1, fig.cap = "13S condensin complex;Condensin I complex"}
out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[1]]
 knitr::include_graphics("my_plot.png")
```
![13S condensin Complex](https://user-images.githubusercontent.com/82643524/236702562-3e656a6d-3272-422f-9419-2c698a1d0a09.png)

2) The second object is the number of Hits detected in the experiment as the number of significant binary interactions. So binare pearson correlations values in the experiment higher than a filter.  

```{r7}
 out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[2]]

[[2]]
  Hits
1   10
```


## The output list out_Hek_P2_1_final_output has 4 objects
3) The name of the proteins involved in the Binary interaction detected (Experimental Protein-binary Pearson correlatio nhigher than the filter) and its Pearson correlation of thes significant binary interactions.

```{r8}
out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[3]]

[[3]]
[[3]]$interactions
 [1] "NCAPH:NCAPD2" "SMC4:NCAPD2"  "NCAPH:NCAPG"  "SMC4:NCAPG"   "NCAPG:NCAPH"  "SMC4:NCAPH"  
 [7] "NCAPG:SMC2"   "SMC4:SMC2"    "NCAPG:SMC4"   "SMC2:SMC4"   

[[3]]$Pearson
 [1] 0.9939944 0.9552057 0.9912144 0.9463247 0.9912144 0.9585936 0.9776026 0.9662482 0.9463247 0.9662482

```

4) Network heatmap constructed from the correlation matrix. 
```{r9 my-plot2, fig.cap = "13S condensin complex;Condensin I complex"}
out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[4]] 
 knitr::include_graphics("my_plot.png")
```
![Network](https://user-images.githubusercontent.com/82643524/236702686-3e03d641-9f98-4d14-a25d-cc96895a664d.png)


# Data Processing (Input of matrix with protein_id as first column and Factions)
  Most of the time we measure samples per duplicate in co-fractionation experiemnt, mCP provide a function to deal with these extensive number of columns. The function input can be a data.frame imported by the following function read.table
```{r10}
NAmatrix_P2_1 <- read.table("Hek293_P2_1.csv",sep =",",dec = ".", header= T)
```
  once imported you will get this file, we provide it in the R package, just write:
```{r11}
 data(NAmatrix_P2_1)
```
##  After that you can run the funtion calc_mean_matrix
  The function calc mean matrix will recognize the replicates by a tag like _A_ and _B_ as part of the names for the MS file and the *Frac_index*. It is the position in which you can find the number of fraction in the column name separated by underscores. for example the frac_index=5 is date_Surname_Measurement_A_01 because the number for the fraction is 5 spaces between underscores, like a_b_c_d_FractionNumber. 
```{r12}
  ###  Note: Here the frac_index= is 17. and the pattern_group= _A_ and pattern_group= _B_  
     names(Hek293_P2_1[,2])
    [1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_A_01_GA1_1_6588.d"
    names(Hek293_P2_1[,37])
    [1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_B_01_GA1_1_6626.d"

```
  
  So we can run the function as follow. 

```{r13}

Hek_293_MA_P2_1 <- Calc_mean_matrix(NAmatrix = NAmatrix_P2_1,
                            pattern_group_A = "_A_",
                            pattern_group_B = "_B_",
                            frac_index=17,
                            Protein_ID_column = 1,
                            save_file = TRUE,
                            save_name = "Hek293_P2_1.csv")
```

Note: The user can also use its own scrip to generate the average of the replicates. This function is limited to 2 replicates. It is our traditional measurement. mCP do not acepts NA values. All NA values in our case are replaced by 0. 

# How does mCP works 

The `mCP` package provides the following functions:

- `Calc_mean_matrix()`: pre-processes the data before running `mcp_list()`.

- `mcp_list()`: Extracts Potential protein complexes in the experiment from CORUM database.
    It extract intensities of all protein complexes from the CORUM database and generates a list of all potential protein complexes presents. Each element of the list is identified as a detected potential protein complex and named with its corresponding complex name. Each protein complex element contains 2 elements: 

 *my_protein_complex_mcp_list[1][1]* is a dataframe composed of the Id proteins of the proteins, names, and fractions and its correspondig intensities per fraction. Names are annotated to the uniprot accesion Ids by a repository R package called gprofiler2. 

 *my_protein_complex_mcp_list[1][2]* is a correlation matrix of these elements,  this correlation matrix is done by Pearson correlations (but kendall and spearman algoriths are also posible in that step). and calculates a matrix correlation between all detected component of each present protein complex.
 
The list of potential protein complexes is the input for the next function called cpp_plotter. This is the first filter for  Then it filters protein complexes composed by at least 2 component with intensities different from cero. Then the first filter is runned: The filter is based in a minimun definition of protien complex. So in this step it keep only protein complexes with at least one significant coelution hit.


- `cpp_plotter()`: The cpp_plotter() function takes a list of potential protein complexes as input. The function first filters out any complexes composed of less than 2 components with non-zero intensities. Then, it applies a filter based on the minimum definition of a protein complex, keeping only complexes with at least one significant co-elution hit. The function also filters the list of complexes based on the presence of co-eluting proteins, requiring at least one binary interaction higher than the filter for proteins within the candidate complex. Finally, cpp_plotter() generates complexome profiling plots, heatmaps, and network plots of the proteins within the selected complexes.

- `cpp_plt_sim()`: Is an intermediate function used in Monte-Carlo-Simulations fdr_mCP(). The function is the same function that cpp_plotter(), but it is prepared to use a vector of different filters according to the average of pearson correllation detected to each protein complex (higher than the first filter). The function requires a list of potential protein complexes as input. The function first filters out any complexes composed of less than 2 components with non-zero intensities. Then, it applies a filter based on the minimum definition of a protein complex, keeping only complexes with at least one significant co-elution hit. The function also filters the list of complexes based on the presence of co-eluting proteins, requiring at least one binary interaction higher than the filter for proteins within the candidate complex. Finally, cpp_plt_sim() generates complexome profiling plots, heatmaps, and network plots of the proteins within the selected complexes.


*my_protein_complex_mpp_ploter[1][1]* The first object is a profile of proteins, fractions vs intensities.
*my_protein_complex_mpp_ploter[1][2]* The second object is the number of Hits detected in the experiment as the number of significant binary interactions. So binare pearson correlations values in the experiment higher than a filter.  
*my_protein_complex_mpp_ploter[1][3]* The name of the proteins involved in the Binary interaction detected.
*my_protein_complex_mpp_ploter[1][3]* Pearson correlation detected higher than the filter of the detected binaries interactions. 
*my_protein_complex_mpp_ploter[1][4]* Network heatmap constructed from the correlation matrix.


- `mcp_fdr()`: This function performs different matrix from your experimental matrix and evaluates this "fake matrix" into  mCP workflow by using the selected filter then calculates the FDR based in montecarlo simulations and the real result of protein complexes detected in your experiment.


