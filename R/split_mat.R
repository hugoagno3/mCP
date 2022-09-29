#' split_mat
#'
#' @param Mat 
#' @param keep 
#' @param pattern 
#' @param Column 
#'
#' @return
#' @export
#'
#' @examples

split_mat <- function(Mat,keep = NULL,pattern,Column = T){
  if(Column){
    return(Mat[,unique(c(keep,grep(pattern = pattern,colnames(Mat))))])
  }
  else{
    return(Mat[unique(c(keep,grep(pattern = pattern,rownames(Mat)))),])
  }
}
