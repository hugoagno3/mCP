#' FDR_mCP
#'
#' @param corum_database
#' @param experiment_data 
#' @param N_fractions 
#' @param filter 
#' @param specie 
#' @param n_simulations
#'
#' @return
#' @export
#'
#' @examples

FRD_mCP<- function (corum_database, experiment_data, N_fractions = 35, 
                    specie = "hsapiens", filter=0.93, n_simulations= 10) 
{
  X<- replicate(n_simulations,{
    Hek_293_fake <- experiment_data
    Hek_293_fake[,2:ncol(experiment_data)]<- experiment_data[sample(1:nrow(experiment_data),nrow(experiment_data),replace = F),2:ncol(experiment_data)]
    
    CL_hek_P1_2_a<- mcp_list(corum_database =  Corum_Humans_Database,
                             experiment_data = Hek_293_fake, N_fractions = N_fractions, specie = specie)
    
    out_Hek_p1_2_a <- cpp_plotter(CL_hek_P1_2_a,output_name = paste0(format(Sys.time(), "%H_%M_%OS3"),"fake.pdf"), format = ".", filter = 0.93, N_fractions = 35,
                                  display_weights = TRUE, standard_weights =list(list(x = 11, label= "1049KDa"), 
                                                                                 list(x = 13, label ="720 KDa"),
                                                                                 list(x = 17, label = "480 KDa"), 
                                                                                 list(x =22, label ="146 KDa"),
                                                                                 list(x =29, label ="60 KDa")))
  }
  )
  Analysis_result  <- mcp_list(corum_database =  Corum_Humans_Database,
                               experiment_data = experiment_data, N_fractions = N_fractions, specie = specie)
  n_Complexes <- sum(sapply(Analysis_result, function(Comp){
    return(sum(Comp$corMat >= filter)>nrow(Comp$corMat))
  }))
  
  fdrs <- sapply(X, "length")/n_Complexes
  print(paste0("FDR (McS)= ", mean(fdrs)))
  FDRsimulations<- X
}
