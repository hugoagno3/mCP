#' fdr_mCP
#' 
#' @description This function performs different matrix from your experimental matrix and evaluates this "fake matrix" into  mCP workflow by using the selected filter then calculates the FDR based in montecarlo simulations and the real result of protein complexes detected in your experiment.
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
#' data(Corum_Humans_Database)
#' 
#' FDR_DIANN_dDIA_P2_1_<- fdr_mCP(corum_database= Corum_Humans_Database,
#'Output_cpp_plotter = out_Hek_P2_1, 
#'experiment_data=Hek293_P2_1,
#'file_name = "m_CP_analysis",
#'N_fractions = 35,
#'specie = "hsapiens",
#'filter=0.81,
#'n_simulations= 2)
#' 

fdr_mCP<-    function (corum_database, experiment_data, N_fractions = 35, 
                       specie = "hsapiens", filter = 0.81, n_simulations = 10, 
                       Output_cpp_plotter= Output_cpp_plotter, file_name = "exp_id", 
                       save_file = TRUE, fdr_limit= 0.05, Risk_fraction = floor(N_fractions*0.85), 
                       monomeric_filter= FALSE, set_seed= TRUE) 
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
      mean(x[[3]][["Cor"]], na.rm = TRUE)
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
  
  #### estimate Monomeric risk #################################################
  Max_positions <- sapply(Output_cpp_plotter,function(x){
    Intensity_profile <- x[[1]][["data"]] %>% group_by(SEC_FR) %>% summarise(AV_profile= mean(Intensity))
    return(which(Intensity_profile$AV_profile == max(Intensity_profile$AV_profile)))
  })
  
  Monomeric_risk_vec <- Max_positions > Risk_fraction
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
  res_DF <- data.frame(complex_names, Monomeric_risk= Monomeric_risk_vec,
                       Hits_in_experiment = standar_Experiment$hits,
                       False_positives = disco_11$Discovered,
                       FDR =disco_11$Discovered/n_simulations, 
                       N_subunits = as.numeric(table(corum_database$complex_name)[complex_names]),
                       rbind(res_matrix)) #,rbind(res_matrix) to include the hits in each simulation
  res_DF<- res_DF %>% dplyr::filter(FDR<= fdr_limit)
  filti<- res_DF$complex_names
  out_Hek_P2_1 <-  Output_cpp_plotter[filti]
  N_detected <- as.data.frame(sapply(out_Hek_P2_1, function(x){return(nrow(unique(x[[1]][["data"]][,"protein_id"])))}))
  colnames(N_detected)<-"N_detected"
  ################################################################################################################
  UniprotIDs <- as.data.frame(sapply(out_Hek_P2_1, function(x){
    return(paste(unique(x[[1]][["data"]][,"protein_id"]), collapse = ";", sep = ""))
  }))
  df_collapsed <- apply(UniprotIDs, 1, function(x) {paste(unlist(eval(parse(text=x))), collapse = ";")})
  UniprotIDs_1<- as.data.frame(df_collapsed)
  colnames(UniprotIDs_1)<-"uni_ids"
  ####################################################################################################
  interx <- as.data.frame(sapply(out_Hek_P2_1, function(x){
    return(paste(unique(x[[3]][["interactions"]][]), collapse = ";", sep = ""))
  }))
  colnames(interx)<-"interac"
  ###################################################################################################
  Pearsx<- as.data.frame(sapply(out_Hek_P2_1, function(x){
    return(paste(unique(x[[3]][["Cor"]][]), collapse = ";", sep = ""))
  }))
  colnames(Pearsx)<-"pearson"
  ################################################################################################################
  Protnames <- as.data.frame(sapply(out_Hek_P2_1, function(x){
    return(paste(unique(x[[1]][["data"]][,"prot_name"]), collapse = ";", sep = ""))
  }))
  df_collapsed_prot <- apply(Protnames, 1, function(x) {paste(unlist(eval(parse(text=x))), collapse = ";")})
  protnames_1<- as.data.frame(df_collapsed_prot)
  colnames(protnames_1)<-"prot_name"
  ################################################################################################################
  porc_detected<- 100*(N_detected/res_DF[,"N_subunits"])
  colnames(porc_detected)[1]<-"porc_detected"                     
  ############################################################
  if (save_file) {
    write.csv(res_DF, file = paste0(file_name, "MontecarloSimulationFDR_results.csv"), row.names = FALSE)
    VS_f<- data.frame(complex_names= res_DF[,"complex_names"],
                      Hits_in_experiment = res_DF[,"Hits_in_experiment"],
                      res_DF[,"FDR"], 
                      Monomeric_risk= res_DF[,"Monomeric_risk"],
                      N_detected = N_detected,
                      N_subunits = res_DF[,"N_subunits"],
                      UniprotID = UniprotIDs_1[,"uni_ids"],
                      Gen_Symbol = protnames_1[,"prot_name"],
                      Ids_binary_Interactions= interx[,"interac"],
                      Correlation_Binary_interactions= Pearsx[,"pearson"],
                      Percentage_detected =porc_detected[,"porc_detected"])
    write.csv(VS_f, file = paste0(file_name,"main.csv"), row.names = FALSE) 
  }
  if (monomeric_filter) {
   
    mono_filter<- VS_f$complex_names[VS_f$Monomeric_risk ==TRUE]
    out_Hek_P2_2<- out_Hek_P2_1 [! (names(out_Hek_P2_1) %in% mono_filter)] 
    
    sink(paste0(file_name, "_FDR.txt"))
    print("Fold discovery rate results")
    print(paste0("FDR (McS)= ", mean(sapply(X, "length")/length(Output_cpp_plotter))))
    print(paste0("PPC_candidate_Detected= ", length(Output_cpp_plotter)))
    print(paste0("PPC_Detected= ", length(out_Hek_P2_2)))
    print(paste0("SD_Fdr= ", sd(sapply(X, "length")/length(Output_cpp_plotter))))
    print(paste0("filter=" , filter))
    print(paste0("number of simulations= ", n_simulations))
    print(paste0("fdr_limit", fdr_limit))
    print(paste0("Monomeric_risk_vec= ", length(VS_f$Monomeric_risk[VS_f$Monomeric_risk==TRUE])))
    sink()
    
      return(out_Hek_P2_2)
  }
  
  sink(paste0(file_name, "_FDR.txt"))
  print("Fold discovery rate results")
  print(paste0("FDR (McS)= ", mean(sapply(X, "length")/length(Output_cpp_plotter))))
  print(paste0("PPC_candidate_Detected= ", length(Output_cpp_plotter)))
  print(paste0("PPC_Detected= ", length(out_Hek_P2_1)))
  print(paste0("SD_Fdr= ", sd(sapply(X, "length")/length(Output_cpp_plotter))))
  print(paste0("filter= ", filter))
  print(paste0("number of simulations= ", n_simulations))
  print(paste0("fdr_limit= ", fdr_limit))
  print(paste0("Monomeric_risk_vec= ", length(VS_f$Monomeric_risk[VS_f$Monomeric_risk==TRUE])))
  sink()

  return(out_Hek_P2_1)
}
