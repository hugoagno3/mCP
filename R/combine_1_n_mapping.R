#' combine_1_n_mapping
#'
#' @param vec_1 
#' @param Mat
#'
#' @return
#' @export
#'
#' @examples

# reduce petide level to protein level
combine_1_n_mapping <- function(vec_1,Mat){
  if(length(vec_1) != nrow(Mat)){
    stop("vec_1 must have same length as number of rows in Mat")
  }
  res_mat <- matrix(numeric(length(unique(vec_1))*ncol(Mat)),nrow = length(unique(vec_1)),
                    ncol = ncol(Mat),dimnames = list(unique(vec_1),colnames(Mat)))
  for (i in 1:nrow(res_mat)) {
    if(sum(vec_1==rownames(res_mat)[i]) == 1){res_mat[i,] <- Mat[vec_1==rownames(res_mat)[i],]}
    else{res_mat[i,] <- colSums(Mat[vec_1==rownames(res_mat)[i],])/sum(vec_1==rownames(res_mat)[i])}
  }
  return(res_mat)
}
