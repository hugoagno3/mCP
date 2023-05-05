#' fdr_mCP
#' 
#' @description This function performs different matrix from your experimental matrix and evaluetes this "fake matrix" into  mCP workflow by using the selected filter then calculates the FDR based in montecarlo simulations and the real result of protein complexes detected in your experiment.
#' @param corum_database The corum data base. 
#' @param experiment_data A matrix that has protein_id in the first column and the detected intensities in the other colums.
#' @param N_fractions Number of fractions in your experimet
#' @param filter Pearson Correlation filter for the analysis
#' @param specie Write: "hsapiens" if you work with homosapiens. 
#' @param n_simulations number of simulations 10 by default
#' @param Output_cpp_plotter output list of detected protein complexes function cpp_plotter
#' @param file_name this is the file name of the output, is a string, It should be writen within "". For example "outputname"
#' @param save_file could be TRUE or FALSE it generate an csv file with results of Fold Discovery Rate per protein complexes
#'
#' @return
#' @export
#'
#' @examples

fdr_mCP<-   function (corum_database, experiment_data, N_fractions = 35, 
                      specie = "hsapiens", filter = 0.93, n_simulations = 10, 
                      Output_cpp_plotter, file_name = "exp_id", save_file = TRUE, fdr_limit= 0.05) 
{
  standar_Experiment<- extract_mcp(Output_cpp_plotter)
  ######## simulation##########
  X <- replicate(n_simulations, {
    experiment_data <- experiment_data
    Indeces_in_Corum <- experiment_data$protein_id %in% 
      corum_database$protein_id
    fake_Ids <- sample(experiment_data$protein_id[!Indeces_in_Corum], 
                       sum(Indeces_in_Corum), replace = F)
    experiment_data$protein_id[experiment_data$protein_id %in% 
                                 fake_Ids] <- experiment_data$protein_id[sample(which(Indeces_in_Corum), 
                                                                                sum(Indeces_in_Corum), replace = F)]
    experiment_data$protein_id[Indeces_in_Corum] <- fake_Ids
    CL_hek_P1_2_a <- mcp_list(corum_database = corum_database, 
                              experiment_data = experiment_data, N_fractions = N_fractions, 
                              specie = specie, heatmap_seaborn = FALSE)
    out_Hek_p1_2_a <- cpp_plotter(complex_list = CL_hek_P1_2_a, 
                                  output_name = paste0(format(Sys.time(), "%H_%M_%OS3"), 
                                                       "fake.pdf"), format = ".", filter = filter, 
                                  N_fractions = N_fractions, display_weights = FALSE)
    #standar_simulation<- extract_mcp(out_Hek_p1_2_a)
  })
  
  complex_names <-  names(Output_cpp_plotter)
  FDs <- sapply(complex_names, function(comp_name) {
    return(sum(sapply(X, function(Simulation) {
      return(comp_name %in% names(Simulation) && Simulation[[comp_name]][[2]][1,1] >= standar_Experiment[comp_name,"hits"])
    })))
  })
  FDs[is.na(FDs)] <- 0
  
  res_matrix <- matrix(0,nrow = length(complex_names),ncol = n_simulations,dimnames = list(complex_names,paste0("Simulation_hits_",1:n_simulations)))
  for(i in 1:n_simulations){
    res_matrix[names(X[[i]])[names(X[[i]]) %in% rownames(res_matrix)],i] <- sapply(sapply(X[[i]],"[[",2),"[[",1)[names(X[[i]]) %in% rownames(res_matrix)]
  }
  ############################################################## 
  res_DF_1 <- data.frame(complex_names,Hits = standar_Experiment$hits, Discovered = FDs, rel.discoveries =FDs/n_simulations, 
                         N_subunits = as.numeric(table(corum_database$complex_name)[complex_names]),rbind(res_matrix)) #,rbind(res_matrix) to include the hits in each simulation
  ################################################################
  ######################################################
  #FDR_P11_today_50023 <- FDR_P11_2_801_tests_100
  VS<-res_DF_1 %>% dplyr::select(6:(5+n_simulations))
  VES<-res_DF_1$Hits
  papa<- as.data.frame(ifelse(VS>=VES,1,0))
  disco_11<- data.frame(Discovered= rowSums(papa))
  ##############################################################
  res_DF <- data.frame(complex_names,Hits_in_experiment = standar_Experiment$hits, False_positives = disco_11$Discovered, FDR =disco_11$Discovered/n_simulations, 
                       N_subunits = as.numeric(table(corum_database$complex_name)[complex_names]),rbind(res_matrix)) #,rbind(res_matrix) to include the hits in each simulation
  res_DF<- res_DF %>% dplyr::filter(FDR<= fdr_limit)
  filti<- res_DF$complex_names
  out_Hek_P2_1 <-  Output_cpp_plotter[filti]
  ############################################################
  if (save_file) {
    write.csv(res_DF, file = paste0(file_name,".csv"), row.names = FALSE)
     VS_f<- res_DF %>% dplyr::select(1,3,4,6:(5+n_simulations))
     write.csv(VS_f, file = paste0(file_name,"MontecarloSimulationFDR_results.csv"), row.names = FALSE)
  }
  sink(paste0(file_name, "_FDR.txt"))
  print("Fold discovery rate results")
  print(paste0("FDR (McS)= ", mean(sapply(X, "length")/length(Output_cpp_plotter))))
  print(paste0("PPC_Detected= ", length(Output_cpp_plotter)))
  print(paste0("SD_Fdr= ", sd(sapply(X, "length")/length(Output_cpp_plotter))))
  sink()
  return(out_Hek_P2_1)
}