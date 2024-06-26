#' extract_mcp 
#' 
#' @description This function Extract data from the output list of mCP and generates a data.frame file 
#' @param z The output list of mCP 
#'
#' @return a dataframe with all de information in the list of plot generated by cpp_plotter function
#' @export
#' @import gprofiler2
#' @import gdata
#' @import corrr
#' @import dplyr
#' @import assertthat
#' @examples extract_mcp(out_Hek_P2_1)
#'#### for  example RUN mCP list that creates a list of potential protein complexes
#'data("Hek293_P2_1")
#'
#'CL_hek_P2_1<- mcp_list(corum_database =  Corum_Humans_Database,
#'                       experiment_data = Hek293_P2_1, 
#'                       N_fractions = 35, 
#'                       specie = "hsapiens",
#'                       method_cor = "pearson",
#'                       heatmap_seaborn = TRUE)
#'out_Hek_P2_1 <- cpp_plotter(complex_list = CL_hek_P2_1,
#'                            format = "pdf", 
#'                            output_name = "m_CP_analysis",
#'                            filter = 0.93,
#'                            N_fractions = 35,
#'                            heat_map = TRUE,
#'                            relative = FALSE,
#'                            display_weights = TRUE,
#'                            standard_weights = list(list(x =11, label= "1049KDa"), 
#'                                                    list(x = 13, label ="720 KDa")))                      
#'extraction<- extract_mcp(out_Hek_P2_1)                       
#'                  
extract_mcp<- function(z){
  printa<- function(z){
    as.data.frame(z[[2]])
  }
  pavor<- sapply(z,printa)
  pavor_2<- as.data.frame(unlist(t(pavor)), row.names=names(pavor))
  colnames(pavor_2)[1]<-"hits"
  #################################
  printi<- function(z){
    as.data.frame(z[[3]])
  }
  lavor<- sapply(z,printi)
  lavor_1<- as.data.frame(unlist(t(lavor)), row.names=names(lavor))
  VERSR<-as.data.frame(t(lavor))
  sds<-as.data.frame(pavor_2$hits)
  colnames(sds)[1]<-"lab"
  vidi<- sds$lab
  VERSR$hits <-vidi  
  VERSR$names<- rownames(VERSR)
  Results_1<- VERSR %>% select(4,3,1,2)
}