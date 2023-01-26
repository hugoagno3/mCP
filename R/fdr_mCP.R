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
#'
#' @return
#' @export
#'
#' @examples

fdr_mCP<- function (corum_database, experiment_data, N_fractions = 35, 
                    specie = "hsapiens", filter=0.93, n_simulations= 10, Output_cpp_plotter, file_name= "exp_id") 
{
  sink(paste0(file_name, "_FDR.txt"))
  print("starting Monte Carlo Simulation")
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
    
    CL_hek_P1_2_a<- mcp_list(corum_database =  Corum_Humans_Database,
                             experiment_data = experiment_data, N_fractions = N_fractions, specie = specie)
    
    out_Hek_p1_2_a <- potter(CL_hek_P1_2_a,output_name = paste0(format(Sys.time(), "%H_%M_%OS3"),"fake.pdf"), format = ".", filter = filter, N_fractions = N_fractions,
                                  display_weights = TRUE, standard_weights =list(list(x = 11, label= "1049KDa"), 
                                                                                 list(x = 13, label ="720 KDa"),
                                                                                 list(x = 17, label = "480 KDa"), 
                                                                                 list(x =22, label ="146 KDa"),
                                                                                 list(x =29, label ="60 KDa")))
  }
  )
  print("Fold discovery rate results")
  print(paste0("FDR (McS)= ", mean(sapply(X, "length")/length(Output_cpp_plotter))))
  print(paste0("PPC_Detected= ", length(Output_cpp_plotter)))
  print(paste0("SD_Fdr= ", sd(sapply(X, "length")/length(Output_cpp_plotter))))
  sink()
}
