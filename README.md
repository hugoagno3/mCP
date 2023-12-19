# mCP


![imagen](https://user-images.githubusercontent.com/82643524/179713559-a0f27ad1-07af-4d69-ac72-380d4eba304b.png)

mini-Complexome Profiler (mCP) 


 "mCP_tutorial_Humans_samples"

# Overview 

The `mCP` package provides a set of functions for targeted protein complex detection from data generated from co-fractionation mass spectrometry-based proteomics experiments. It is designed to work with experimental data in the form of protein abundance matrices. The program, identifies protein complexes present in the experimental data, using an external protein complexes database (e.g. CORUM or Complex Portal). The R package contains co-fractionation experiment datasets as examples of two different species: Homo sapiens (Hek293T cells) and Mus musculus (cardiomyocytes). 

This vignette provides a step-by-step guide to using the `mCP` package to analyze proteomics data and identify protein complexes of interest.
In this tutorial, we will analyze a single co-fractionation experiment. An experimental dataset "Homo sapiens", derived from Hek293T cells, fractionated by BNE-PAGE (35 fractions).
In addition, we provide a detailed explanation of all of the R package functions at the end of this tutorial.

https://private-user-images.githubusercontent.com/82643524/291553448-c0924ec2-5110-4f5d-8ccd-45d57b53a14f.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDI5ODExODUsIm5iZiI6MTcwMjk4MDg4NSwicGF0aCI6Ii84MjY0MzUyNC8yOTE1NTM0NDgtYzA5MjRlYzItNTExMC00ZjVkLThjY2QtNDVkNTdiNTNhMTRmLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFJV05KWUFYNENTVkVINTNBJTJGMjAyMzEyMTklMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjMxMjE5VDEwMTQ0NVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTA0N2QwNmIzOTIxYmQ5OGJkY2I0NTdlOWEzMmE4Y2ZmMGM0NzVkY2Q0MjU5ODY2ZDFmMGU2MWZiMjcwZDFkMzEmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.hlkuJjIf2jRzjuzGI6XTPO7azImljU39Ijlv8bGfuts

![Rplot02](https://github.com/hugoagno3/mCP/assets/82643524/0f62593f-a702-4330-98f9-c569c9f46ca8)

![heat eIF3 complex](https://user-images.githubusercontent.com/82643524/236703288-735a5d6e-25aa-49b9-8520-071a3df7ea4c.png)
![Eif3j heatmap](https://user-images.githubusercontent.com/82643524/236703265-8cd3be1e-4252-4878-8f91-58481916cb84.png)

# Installation

To install the mCP package from GitHub, run the following code:
```{r pressure, echo=FALSE}
library(devtools)
devtools::install_github("hugoagno3/mCP", biuld_vignette= TRUE)
```

# Data input 
 In this tutorial, we will use the Corum_Humans_Database file as the external protein complex reference database, and the Hek293_P2_1 file as experimental input. To perform targeted protein complex detection, only these two input files are needed:
```{r pressure, echo=FALSE}
# Open the Corum protein complex database file
      data(Corum_Humans_Database)
# open the experiment files
      data(Hek293_P2_1)
```
* A CORUM database of protein complexes (or a targeted list of interest) data frame with 3 columns:´  col1= "complex_id", col2= "complex_name" and col3= "protein_id").
* An experiment file *data.frame* with your experimental results in wide format. The first column is called "protein_id", with a single UniProt accession per row as elements. The rest of the columns contain protein intensities detected in each fraction of the co-fractionation experiment. The column's names of the fractions. It can be numeric names from 1 to the last number of fractions. For example, 1,2,3,..,35 for 35 fractions, means 35 columns.

Note: The `mCP` R package is focused on the detection of protein complexes and it accepts only protein data Level as an input matrix. The mass spectrometry data acquisition can be done by Data Dependent Acquisition mode or Data independent Aquisition.  

# Data Processing

To process the input data, we need to run *Option 1*- the mCP function or *Option 2*- 3 functions provided by our `mCP` package.
# Option 1: Running mCP function
   mCP() is an integrated function of `mCP` R-package. It requires an experimental protein abundance matrix dataset as input, and returns the following: a curated list of identified protein complexes, accompanied by a controlled False Discovery Rate filter (FDR). Each detected protein complex is an output element organized into the resulting list output of mCP. The output list is structured with four distinct elements.
   1- One fraction profile plot (absolute or relative abundance of each protein vs. fraction).
   2- The number of intra binary interactions (protein-protein) with a Pearson correlation value higher than the filter (number of total binary hits) within component of each protein complex in the experiment.
   3- The id of proteins of binary hits.
   4- A heatmaps_seaborn of known protein complexes detected in the protein complexes database, in this tutorial we use the CORUM database (http://mips.helmholtz-muenchen.de/corum/).
   The above data is outputted in the following format: 
 1- pdf file with fraction-profile plots (absolute or relative abundance of each protein vs. fraction) of each potential protein complex to be detected. Before the FDR analysis.
 2- pdf with heatmaps of potential candidates from 1 (before the FDR analysis).
 3- pdf file with the detected protein complexes profiles with a controlled FDR (this is analogous to 1, but without protein complexes, excluded based on FDR evaluation).
 4- pdf with heatmaps of the detected protein complexes with a controlled FDR (this is analogous to 2, but without protein complexes excluded based on FDR evaluation).
 5- txt file with numbers about general false positives.
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
![Respirationry chain complex I for github](https://github.com/hugoagno3/mCP/assets/82643524/c7b4d75f-3e35-4a49-9f20-ea85d887c79b)

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
The following code guides you to get the output step-by-step

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
##Note this function will perform an FDR filter and simulation it is important to do the simulation with the same filter then the previous function. If the filter value varies, the simulation loses its meaningfulness. 
```
## The last step to get the protein complexes filtered by FDR is to get back to the list and run cpp_plotter again.
Note that if you wish to have the plots of this last FDR protein complexes you could filter the names of the protein complexes detected into the first list (from mcp_list) and run cpp_ploter again. 


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
![13s condensin complex for hithub](https://github.com/hugoagno3/mCP/assets/82643524/9bfe668f-a8d2-4ad3-b335-e2dc8243dc5a)

2) The second object is the number of hits detected in the experiment as the number of significant binary interactions. So, binary Pearson's correlation values in the experiment that are higher than a filter.  

```{r7}
 out_Hek_P2_1_final_output[["13S condensin complex;Condensin I complex"]][[2]]

[[2]]
  Hits
1   10
```


## The output list out_Hek_P2_1_final_output has 4 objects
3) The name of the proteins involved in the Binary interaction detected (Experimental Protein-binary Pearson correlation higher than the filter) and its Pearson correlation of the significant binary interactions.

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
  In a typical co-fractionation experiment, samples are measured per duplicate. mCP R package provides a function to deal with this extensive number of columns. The function input can be a data.frame imported by the following function read.table
```{r10}
NAmatrix_P2_1 <- read.table("Hek293_P2_1.csv",sep =",",dec = ".", header= T)
```
  once imported you will get this file, we provide it in the R package, just write:
```{r11}
 data(NAmatrix_P2_1)
```
##  After that you can run the function calc_mean_matrix
  The function calc mean matrix will recognize the replicates by a tag like _A_ and _B_ as part of the names for the MS file and the *Frac_index*. It is the position in which you can find the number of fractions in the column name separated by underscores. For example, the frac_index=5 is date_Surname_Measurement_A_01 because the number for the fraction is 5 spaces between underscores, like a_b_c_d_FractionNumber. 
```{r12}
  ###  Note: Here the frac_index= is 17. and the pattern_group= _A_ and pattern_group= _B_  
     names(Hek293_P2_1[,2])
    [1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_A_01_GA1_1_6588.d"
    names(Hek293_P2_1[,37])
    [1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_B_01_GA1_1_6626.d"

```
  
  So we can run the function as follows. 

```{r13}

Hek_293_MA_P2_1 <- Calc_mean_matrix(NAmatrix = NAmatrix_P2_1,
                            pattern_group_A = "_A_",
                            pattern_group_B = "_B_",
                            frac_index=17,
                            Protein_ID_column = 1,
                            save_file = TRUE,
                            save_name = "Hek293_P2_1.csv")
```

Note: The user can also use their own script to generate the average of replicate experiments. At the moment, this function is limited to 2 replicates. mCP function does not accept NA values. All NA values in our case are replaced by 0. 

# How mCP works 

The `mCP` package provides the following functions:

- `Calc_mean_matrix()`: pre-processes the data before running `mcp_list()`.

- `mcp_list()`: Extracts potential protein complexes from the experimental dataset, using the protein complexes database as a reference.
    This function extracts intensities of all protein complexes from the protein complexes database (CORUM or complex portal) and generates a list of all potential protein complexes present in the experimental dataset. Each element of the list is identified as a detected potential protein complex and named with its corresponding complex name. Each protein complex element contains 2 elements: 

 *my_protein_complex_mcp_list[1][1]* is a dataframe composed of the Id proteins of the proteins, names, and fractions and their corresponding intensities per fraction. Names are annotated to the UniProt accesion ids by a repository R package called gprofiler2. 

 *my_protein_complex_mcp_list[1][2]* is a correlation matrix of these elements,  this correlation matrix is done by Pearson correlations (but Kendall and Spearman algorithms are also possible in that step). The matrix was calculated as a correlation matrix of all detected component of each detected protein complex.
 
The list of potential protein complexes is the input for the next function called `cpp_plotter()`. This is the presearch stage. Where protein complexes composed by at least 2 component with intensities different from cero are selected. Then the first filter is runned: The filter is based in a minimun definition of protien complex. So in this step it keep only protein complexes with at least one significant coelution hit (a binary interaction with a Pearson correlation higher than the filter).

- `cpp_plotter()`: The cpp_plotter() function takes a list of potential protein complexes as input. The function first filters out any complexes composed of less than 2 components with non-zero intensities. Then, it applies a filter based on the minimum definition of a protein complex, keeping only complexes with at least one significant co-elution hit. The function also filters the list of complexes based on the presence of co-eluting proteins, requiring at least one binary interaction higher than the filter for proteins within the candidate complex. Finally, cpp_plotter() generates complexome profiling plots, heatmaps, and network plots of the proteins within the selected complexes.

- `cpp_plt_sim()`: Is an intermediate function used in Monte-Carlo-Simulations fdr_mCP(). The function is the same function that cpp_plotter(), but it is prepared to use a vector of different filters according to the average of Pearson correlation detected to each protein complex (higher than the first filter). The function requires a list of potential protein complexes as input. The function first filters out any complexes composed of less than 2 components with non-zero intensities. Then, it applies a filter based on the minimum definition of a protein complex, keeping only complexes with at least one significant co-elution hit. The function also filters the list of complexes based on the presence of co-eluting proteins, requiring at least one binary interaction higher than the filter for proteins within the candidate complex. Finally, cpp_plt_sim() generates complexome profiling plots, heatmaps, and network plots of the proteins within the selected complexes.


*my_protein_complex_mpp_ploter[1][1]* The first object is a profile of proteins, fractions vs intensities.
*my_protein_complex_mpp_ploter[1][2]* The second object is the number of Hits detected in the experiment as the number of significant binary interactions. So, binary Pearson correlation values in the experiment that are higher than a filter.  
*my_protein_complex_mpp_ploter[1][3]* The name of the proteins involved in the Binary interaction detected.
*my_protein_complex_mpp_ploter[1][3]* Pearson correlation detected higher than the filter of the detected binaries interactions. 
*my_protein_complex_mpp_ploter[1][4]* Network heatmap constructed from the correlation matrix.

- `mcp_fdr()`: This function creates different decoys matrixes from your experimental matrix and evaluates this "decoy matrix" into  mCP workflow by using the selected filter then calculates the FDR based in montecarlo simulations. The characteristics of the decoys are compared agains the candidates in each simulations. A decoy matrix is consider as a potential false positive when the number of hits of Pearson correlation are higher that the candidates. Then the false discovery rate is calculates as False positive/number of simulations for each protein complex. 


