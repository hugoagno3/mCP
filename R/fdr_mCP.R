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

fdr_mCP<-  function (corum_database, experiment_data, N_fractions = 35, 
                     specie = "hsapiens", filter=0.93, n_simulations= 10, Output_cpp_plotter, file_name= "exp_id", save_file=TRUE) 
{
  #sink(paste0(file_name, "_FDR.txt"))
  #print("starting Monte Carlo Simulation")
  X<- replicate(n_simulations,{
    experiment_data <- experiment_data
    # Identifie those Protein Ids that are in the data set and in the Corum DB
    Indeces_in_Corum <- experiment_data$protein_id %in% corum_database$protein_id
    # choose random Protein Ids from Proteins, that are not in the Corum DB
    fake_Ids <- sample(experiment_data$protein_id[!Indeces_in_Corum],sum(Indeces_in_Corum),replace = F)
    # replace choosen Prot ids, from proteins that are not part of a complex with Prot Ids from Prot that are in the Corum DB
    experiment_data$protein_id[experiment_data$protein_id %in% fake_Ids] <- experiment_data$protein_id[sample(which(Indeces_in_Corum),sum(Indeces_in_Corum),replace = F)]
    # replace Prot Ids, of proteins that are in the Corum DB with random Protein Ids from Proteins that are not in the Corum DB
    experiment_data$protein_id[Indeces_in_Corum] <- fake_Ids
    
    # this part we don't need anymore
    #experiment_data[,2:ncol(experiment_data)]<- experiment_data[sample(1:nrow(experiment_data),nrow(experiment_data),replace = F),2:ncol(experiment_data)]
    
    CL_hek_P1_2_a<- mcp_list(corum_database =  corum_database,
                             experiment_data = experiment_data, N_fractions = N_fractions, specie = specie)
    
    out_Hek_p1_2_a <- cpp_plotter(CL_hek_P1_2_a,output_name = paste0(format(Sys.time(), "%H_%M_%OS3"),"fake.pdf"), format = ".", filter = filter, N_fractions = N_fractions,
                                  display_weights = FALSE)
  }
  )
  
  complex_names <- unique(corum_database$complex_name[corum_database$protein_id %in% experiment_data$protein_id])
  FDs <- sapply(complex_names,function(comp_name){
    return(sum(sapply(X,function(Simulation){
      return(comp_name %in% names(Simulation))
    })))
  })
  
  res_DF <- data.frame(complex_names,
                       Discovered = FDs,
                       rel.discoveries = FDs/n_simulations,
                       N_subunits = as.numeric(table(corum_database$complex_name)[complex_names]))
  
  if(save_file){
    write.csv(res_DF[names(Output_cpp_plotter),],file = paste0(file_name,".csv"))
  }
  
  
  return(res_DF[names(Output_cpp_plotter),])
  
  #print("Fold discovery rate results")
  #print(paste0("FDR (McS)= ", mean(sapply(X, "length")/length(Output_cpp_plotter))))
  #print(paste0("PPC_Detected= ", length(Output_cpp_plotter)))
  #print(paste0("SD_Fdr= ", sd(sapply(X, "length")/length(Output_cpp_plotter))))
  #sink()
  #return(X)
}
