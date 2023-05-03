#' List of Protein Complexes
#' 
#' @description This function detects a list of potential list of complexes... 
#' @param corum_database 
#' @param experiment_data A *data.frame* with your experiment results 
#' @param N_fractions 
#' @param organism = could be "mmusculus", "hsapiens" , check gconvert vignette for more options https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
#' @param method_cor = can be "kendall", "spearman" or "pearson" (default "pearson")
#' @param heatmap_seaborn always TRUE perform a correlation matrix compatible with the heat map seaborn form corrr package. Not active only in the FDR function.
#'
#' @return
#' @export
#' @import gprofiler2
#' @import assertthat
#' @examples Co-fractionation experiments with Digitonin detergent
mcp_list <- function(corum_database, experiment_data, N_fractions = 34, specie = "mmusculus", method_cor= "pearson", heatmap_seaborn= TRUE) {
  
  # datacleaning
  # - check user input
  # - NA, remove duplicates
  # join dataframes by protein_id
  
  #assertthat::assert_that(class(corum_database) == "data.frame")
  # ... 
  Matrix_Clast1 <-
    experiment_data %>%
    dplyr::distinct(protein_id, .keep_all = TRUE) %>% 
    dplyr::select(1:(N_fractions + 1))
  colnames(Matrix_Clast1) <- c("protein_id", 1:N_fractions)
  df <- dplyr::inner_join(corum_database, Matrix_Clast1, by = 'protein_id')
  
  
  
  prot_names <- gprofiler2::gconvert(
    query = df$protein_id,
    organism =specie ,
    target = "ENSG",
    mthreshold = 1,
    filter_na = FALSE
  ) %>% 
    dplyr::select(name, input) %>% #TODO: include description?
    dplyr::mutate(name = ifelse(is.na(name), yes = input, no = name)) %>% # replace NA
    dplyr::rename(protein_id = input)
  
  assertthat::are_equal(prot_names$protein_id, df$protein_id)
  df1 <- df %>% tibble::add_column(prot_name = prot_names$name, .before = 3)
  
  df_long <- df1 %>% 
    tidyr::gather("SEC_FR", "Intensity", as.character(1):as.character(N_fractions)) %>% 
    dplyr::arrange(complex_id)
  df_long$SEC_FR <- as.numeric(df_long$SEC_FR)
  
  output <- unique(df_long$complex_id) %>%
    lapply(function(id) {
      subset <- df_long %>% dplyr::filter(complex_id == id)
      
      corMat <- subset %>% 
        dplyr::select(!protein_id) %>%  distinct(prot_name, SEC_FR,.keep_all = TRUE) %>% 
        tidyr::spread(key = "prot_name", value = "Intensity") %>%
        dplyr::select(!c(complex_id, complex_name, SEC_FR)) %>%
        cor(method= method_cor)
      if (heatmap_seaborn & method_cor=="pearson") {
        corMat_rrr <- subset %>% 
          dplyr::select(!protein_id) %>%  distinct(prot_name, SEC_FR,.keep_all = TRUE) %>% 
          tidyr::spread(key = "prot_name", value = "Intensity") %>%
          dplyr::select(!c(complex_id, complex_name, SEC_FR)) %>%
          corrr::correlate(method= method_cor, diagonal = 1)
      }
      ifelse(heatmap_seaborn & method_cor=="pearson",return(list(data = subset, corMat = corMat, CorMat_rrr= corMat_rrr)),return(list(data = subset, corMat = corMat)))
      #  return(list(data = subset, corMat = corMat, CorMat_rrr= corMat_rrr))
    })
  extri<- function(x){
    a<-unique(x[][["data"]]$complex_name)
    return(a)
  }
  VER<-sapply(output, extri)
  names(output)<- VER
  return(output)
}

