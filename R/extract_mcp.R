#' List of Protein Complexes
#' 
#' @description This function Extract data from the output list of mCP and generates a data.frame file 
#' @param z the output list of mCP 
#'
#' @return
#' @export
#' @import gprofiler2
#' @import assertthat
#' @examples extratract_mcp(output_mcp_list)
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