#' fdr_mCP_connections
#' 
#' @description This function in an internal function of fdr_mCP_standard_modified(). This function needs as inputs a list of complexes outputs generates by mCP() or cpp_ploter, the experimental matrix and a Target list of protein complexes (CORUM or complex portal adapted). It performs the first simulation to detect number of binary interactions that are higher than a clustering filter (pearson). It builds different matrices from your experimental matrix and evaluates this "fake matrix" into  mCP workflow by using the selected filter. The output is numeric matrix with rows named as protein complexes that indicates detected false positives (1)  or false negative (0). 
#' @param corum_database The corum data base. 
#' @param experiment_data A matrix that has protein_id in the first column and the detected intensities in the other colums.
#' @param N_fractions Number of fractions in your experimet
#' @param filter Pearson Correlation filter for the analysis
#' @param specie Write: "hsapiens" if you work with homosapiens. 
#' @param n_simulations number of simulations 10 by default
#' @param Output_cpp_plotter output list of detected protein complexes function cpp_plotter
#' @param file_name this is the file name of the output, is a string, It should be writen within "". For example "outputname"
#' @param save_file could be TRUE or FALSE it generate an csv file with results of Fold Discovery Rate per protein complexes
#' @param Risk_fraction  is a number indicating the fraction where all monomeric components of the protein complex should be present. It is set to 85% of the n_fraction. 
#' @param monomeric_filter a TRUE or FALSE setting that takes out potential protein complexes sie a hight porcentaje of monomeric conformation. This is related to the previous filter by defoult is off so the user can decide. 
#' @param set_seed a TRUE or FALSE setting to standarise simulations. 
#'
#' @return a list of controlled FDr Portein compleses. 
#' @import gprofiler2
#' @import gdata
#' @import corrr
#' @import dplyr
#' @import assertthat
#' @import ggplot2
#' @export
#'
#' @examples
#' 
#' data("Corum_Humans_Database")
#' data("Hek293_P2_1")
#' mCP_Hek_P2_1 <- mCP(corum_database = Corum_Humans_Database,
#'                                      experiment_data = Hek293_P2_1, 
#'                                      N_fractions = 35, 
#'                                      specie = "hsapiens",
#'                                      method_cor = "pearson",
#'                                      network = TRUE,
#'                                      format = "pdf", 
#'                                      output_name = "mCP_Hek293_P2_1",
#'                                      filter = 0.81,
#'                                      heat_map = TRUE,
#'                                      relative = FALSE,
#'                                      fdr_limit = 0.05,
#'                                      n_simulations= 7,
#'                                      monomeric_filter = FALSE,
#'                                      mw = TRUE,
#'                                      dynamic = TRUE,
#'                                      Risk_fraction = 31,
#'                                      set_seed = FALSE,
#'                                      display_weights = TRUE,
#'                                      standard_weights = list(list(x =6, label= "2700 KDa"), 
#'                                                         list(x = 11, label ="950 KDa"),
#'                                                         list(x = 14, label = "750 KDa"), 
#'                                                         list(x =27, label ="146 KDa"),
#'                                                         list(x =30, label ="60 KDa")))
#'
#'  #### Run the function
#'FDR_DIANN_dDIA_P2_12_<- fdr_mCP_connections (corum_database= Corum_Humans_Database,
#'                                             Output_cpp_plotter = mCP_Hek_P2_1, 
#'                                             experiment_data= Hek293_P2_1,
#'                                             N_fractions = 35,
#'                                             specie = "mmusculus",
#'                                             filter=0.81,
#'                                             n_simulations= 4,
#'                                             fdr_limit = 0.05,
#'                                             set_seed= TRUE)
#'
#'


fdr_mCP_connections<-    function (corum_database, experiment_data, N_fractions = 35, 
                       specie = "hsapiens", filter = 0.81, n_simulations = 10, 
                       Output_cpp_plotter= Output_cpp_plotter, fdr_limit= 0.05,
                        set_seed= TRUE) 
{
  ifelse(length(names(Output_cpp_plotter))>0, standar_Experiment<- extract_mcp(Output_cpp_plotter), return(print("0 Protein Complexes detected")))
  complex_names <-  names(Output_cpp_plotter)
  ######## simulation##########
  if (set_seed){
        set.seed(123)
  }
  #####################Specific filter for FRD according to the average of Person binary interaction of each protein complex. 
  e <- function(z) {
    means <- lapply(z, function(x) {
      min(x[[3]][["Cor"]], na.rm = TRUE)
    })
    return(means)
  }
  ee<-e (Output_cpp_plotter)
  
  
  X <- replicate(n_simulations, {
    Output_cpp_plotter<- Output_cpp_plotter
    experiment_data <- experiment_data
    Indeces_in_Corum <- experiment_data$protein_id %in% 
      corum_database$protein_id
    fake_Ids <- sample(experiment_data$protein_id[!Indeces_in_Corum], 
                       sum(Indeces_in_Corum), replace = F)
    experiment_data$protein_id[experiment_data$protein_id %in% 
                                 fake_Ids] <- experiment_data$protein_id[sample(which(Indeces_in_Corum), 
                                                                                sum(Indeces_in_Corum), replace = F)]
    experiment_data$protein_id[Indeces_in_Corum] <- fake_Ids
    
##### Make an mCP list function, that will put out the complexes if the number of hits is lower than the arbitrary filter in comparison to the complexes if that is not the case it will put the average if not it will put 0, that let the complex out.
    CL_hek_P1_2_a <- mcp_list(corum_database = corum_database, 
                              experiment_data = experiment_data, N_fractions = N_fractions, 
                              specie = specie, network = FALSE)
    CL_hek_P1_2_a <- CL_hek_P1_2_a [names(Output_cpp_plotter)]
    out_Hek_p1_2_a <-cpp_plt_sim (complex_list = CL_hek_P1_2_a, Output_cpp_plotter= Output_cpp_plotter, #lapply(ee,filti)
                                  output_name = paste0(format(Sys.time(), "%H_%M_%OS3"),
                                                       "fake.pdf"), format = ".", 
                                  N_fractions = N_fractions, display_weights = FALSE)                    
 
  })
  FDs <- sapply(complex_names, function(comp_name) {
    return(sum(sapply(X, function(Simulation) {
      return(comp_name %in% names(Simulation) && Simulation[[comp_name]][[2]][1,1] >= standar_Experiment[comp_name,"hits"])
    })))
  })
  FDs[is.na(FDs)] <- 0
  res_matrix <- matrix(0,nrow = length(complex_names),ncol = n_simulations,dimnames = list(complex_names,paste0("Simulation_hits_",1:n_simulations)))
  
  for (i in 1:n_simulations) {
    
    subset_rows <- names(X[[i]])[names(X[[i]]) %in% rownames(res_matrix)]
    subset_values <- sapply(sapply(X[[i]], "[[", 2), "[[", 1)[names(X[[i]]) %in% rownames(res_matrix)]
    
    res_matrix[subset_rows, i] <- ifelse(length(subset_values) > 0, subset_values, 0)
  }
  return(res_matrix)
}
  