#' mCP function: Targeted analysis of protein complexes extended to a low number of fractions. 
#'
#' @description This function is an integrated function of mCP package, that needs as input an experimental data and returns a list of plots, binary total hits, id of proteins of binary hits and heatmaps_seaborn of know protein complexes detected in CORUM database. In addition, it plots 4 files as  outputs: 
#' 1- pdf file with detected protein complexes profiles from Corum database.
#' 2- pdf with heatmaps of the detected protein complexes.
#' 3- txt file with numbers about general false positive when just 1 hit is consider as filter.
#' 4- CVS file containing all protein complexes detected, hits of binary interactions inside the protein complexes, FDR detected by MonteCarloSimulation.
#'   
#' @param corum_database a data.frame with 3 columns: first "complex_id", "complex_name" and "protein_id"
#' @param experiment_data A *data.frame* with your experiment results in wide-format- first column called "protein_id" and the next columns numerics names from 1 to number of fractions.  
#' @param N_fractions number of protein fractions obtained in the co-fractionation experiment. 
#' @param specie = could be "mmusculus", "hsapiens" , check gconvert vignette for more options (gprofiler2 package) https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
#' @param method_cor = can be "kendall", "spearman" or "pearson" (default "pearson")
#' @param heatmap_seaborn always TRUE perform a correlation matrix compatible with the heat map seaborn form corrr package. Not active only in the FDR function.
#' @param format should be "pdf" to plots pdf. There is an additional possible value is ".": this will omit pdf. In both cases the user will have a list of plot for the detected protein complexes in R. 
#' @param output_name main core name for the files to be printed, it is a string: see example. 
#' @param filter targeted value of minimum accepted binary interaction hits between 2 proteins within a protein complex (we recomend a 0.93 value). 
#' @param heat_map a TRUE or FALSE statement to plot heatmaps. 
#' @param relativea TRUE or FALSE statement to plot relative complexome profiling plots, useful to compare between conditions.
#' @param display_weights TRUE or FALSE statement to display molecular weight markers in the complexome profiling plots. 
#' @param standard_weights following core= list(list(x =11, label= "1049KDa"), list(x = 13, label ="720 KDa"))). It is possible to add many markers. you have to extend the code, for example to add a thrid marker= list(list(x =11, label= "1049KDa"), list(x=12, label="900 kDa"), list(x = 13, label ="720 KDa"))). Display_weights muss be TRUE 
#' @param fdr_limit report only protein complexes detected with a false discovery rate lower than a numeric value, for example = 0.05. 
#' @param n_simulations This is the number of simulations we recomend 185 simulations (this part could take 12 hrs for human dataset).
#' @param Risk_fraction  is a number indicating the fraction where all monomeric components of the protein complex should be present. It is set to 85% of the n_fraction. 
#' @param monomeric_filter a TRUE or FALSE setting that takes out potential protein complexes sie a hight porcentaje of monomeric conformation. This is related to the previous filter by defoult is off so the user can decide. 
#' @param set_seed a TRUE or FALSE setting to standarise simulations.
#' @param mw TRUE or FALSE setting plot moleculare weight (requires std weights in kDa). 
#'
#' @return a list of protein complexes plots, the significant interactions, binary pairs proteins detected, and a network heatmap based.
#' @export
#' @import gprofiler2
#' @import gdata
#' @import corrr
#' @import dplyr
#' @import assertthat
#' @examples 
#' 
#' Co-fractionation experiments with NP-40 mild-detergent_Humans(Hek293 cells)
#' 
#'data("Corum_Humans_Database")
#'data("Hek293_P2_1")
#'mCP_Hek_P2_1 <- mCP(corum_database = Corum_Humans_Database,
#'                                      experiment_data = Hek293_P2_1, 
#'                                      N_fractions = 35, 
#'                                      specie = "hsapiens",
#'                                      method_cor = "pearson",
#'                                      heatmap_seaborn = TRUE,
#'                                      format = "pdf", 
#'                                      output_name = "mCP_Hek293_P2_1",
#'                                      filter = 0.81,
#'                                      heat_map = TRUE,
#'                                      relative = FALSE,
#'                                      fdr_limit = 0.05,
#'                                      n_simulations= 7,
#'                                      monomeric_filter = FALSE,
#'                                      mw = TRUE,
#'                                      Risk_fraction = 31,
#'                                      set_seed = FALSE,
#'                                      display_weights = TRUE,
#'                                      standard_weights = list(list(x =6, label= "2700 KDa"), 
#'                                                         list(x = 11, label ="950 KDa"),
#'                                                         list(x = 14, label = "750 KDa"), 
#'                                                         list(x =27, label ="146 KDa"),
#'                                                         list(x =30, label ="60 KDa")))
#'                                                              
#'Co-fractionation experiments with C12E8 mild-detergent_Mouse (89FR)
#'
#'data("Corum_Mouse_Database")
#'data("CM_LV_1")
#'mCP_TEND_out_WT_Mouse_Cardiac <- mCP(corum_database = Corum_Mouse_Database,
#'                                     experiment_data = CM_LV_1, 
#'                                     N_fractions = 35, 
#'                                     specie = "mmusculus",
#'                                     method_cor = "pearson",
#'                                     heatmap_seaborn = TRUE,
#'                                     format = "pdf", 
#'                                     output_name = "INTO_m_CP_mouse_analysis_8",
#'                                     filter = 0.81,
#'                                     heat_map = TRUE,
#'                                     relative = FALSE,
#'                                     fdr_limit = 0.05,
#'                                     n_simulations= 14,
#'                                     monomeric_filter = FALSE,
#'                                     Risk_fraction = 31,
#'                                     set_seed = TRUE,
#'                                     mw = TRUE
#'                                     display_weights = TRUE,
#'                                     standard_weights =  list( 
#'                                     list(x = 9, label = "1048 KDa"), 
#'                                     list(x = 11.5, label = "720 KDa"),
#'                                     list(x = 16, label = "480 KDa"),
#'                                     list(x = 21, label = "242 KDa"),
#'                                     list(x = 25, label = "242 KDa"),
#'                                     list(x = 29, label = "166 KDa")))


mCP <- function(corum_database, experiment_data, N_fractions=35, specie= "hsapiens",
                method_cor="pearson", heatmap_seaborn= TRUE, format="pdf", output_name= mCP_analysis,
                filter=0.81, heat_map= TRUE, relative= FALSE, display_weights=TRUE, 
                standard_weights=TRUE, fdr_limit=0.05 , n_simulations=185,
                Risk_fraction=floor(N_fractions*0.85) , monomeric_filter= FALSE, set_seed= TRUE, mw= FALSE){
  
  # initialiye progress bar
  
  pb <- txtProgressBar(min= 0, max= 10, style =  3)
  
  # Record the start time
  start_time <- Sys.time()
  
  # step1: Run mcp_list function
  #pb$tick(msg = "Running mcp_list function...")
  CL_hek_P2_1<- mcp_list(corum_database = corum_database,
                         experiment_data = experiment_data, 
                         N_fractions = N_fractions, 
                         specie = specie,
                         method_cor = method_cor, 
                         heatmap_seaborn = heatmap_seaborn)
  # Print message and elapsed time
  cat("mcp_list function completed in", difftime(Sys.time(), start_time), "seconds.\n")
  
  # Update progress bar
  t<-10/(1.5+5.4+3.75*n_simulations+5.4)
  setTxtProgressBar(pb, 1.5*t)
  
  # Step 2:  Run cpp_plotter function
  #pb$tick(msg = "Running cpp_plotter function...")
  
  out_Hek_P2_1 <- cpp_plotter(complex_list = CL_hek_P2_1,
                              format = format, 
                              output_name = paste0("all_Candidates_", output_name),
                              filter = filter, N_fractions = N_fractions,
                              heat_map = heat_map,
                              heatmap_seaborn= heatmap_seaborn,
                              relative = relative,
                              display_weights = display_weights,
                              standard_weights = standard_weights, mw=mw)
  # Print message and elapsed time
  cat("cpp_plotter function completed in", difftime(Sys.time(), start_time), "seconds.\n")
  # Step 3: Run fdr_mCP function
  #pb$tick(msg = "Running fdr_mCP function...")
  
  # Update progress bar
  setTxtProgressBar(pb, 5.4*t)
  
  FDR_DIANN_dDIA_P2_1_<- fdr_mCP(corum_database = corum_database,
                                 Output_cpp_plotter = out_Hek_P2_1, 
                                 experiment_data = experiment_data,
                                 file_name = output_name,
                                 N_fractions = N_fractions,
                                 specie = specie,
                                 filter = filter,
                                 fdr_limit = fdr_limit,
                                 n_simulations = n_simulations,
                                 Risk_fraction= Risk_fraction,
                                 monomeric_filter= monomeric_filter,
                                 set_seed= set_seed)
                                 
  
  filist<- names(FDR_DIANN_dDIA_P2_1_)
  CL_list_fdr<-CL_hek_P2_1[filist]
  # Update progress bar
  setTxtProgressBar(pb, 5.4*t)
  out_Hek_P2_1 <- cpp_plotter(complex_list = CL_list_fdr,
                              format = format, 
                              output_name = output_name,
                              filter = filter, N_fractions = N_fractions,
                              heat_map = heat_map,
                              heatmap_seaborn= heatmap_seaborn,
                              relative = relative,
                              display_weights = display_weights,
                              standard_weights = standard_weights,mw=mw)
  
  # # Calculate the elapsed time and print it to the console
  # end_time <- Sys.time()
  # elapsed_time <- end_time - start_time
  # cat("Elapsed time:", round(elapsed_time, 2), "mins\n")
  
  # Print message and elapsed time
  cat("fdr_mCP function completed in", difftime(Sys.time(), start_time), "seconds.\n")
  # Update progress bar
  setTxtProgressBar(pb, 10) 
  # End progress bar
  close(pb)
  return(out_Hek_P2_1)
}