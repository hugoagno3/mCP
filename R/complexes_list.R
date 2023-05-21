#' List of Protein Complexes
#' 
#' @description This function produces a list of potential protein complexes. Its input is an experimental matrix and the CORUM database.
#' @param corum_database The CORUM database of protein complexes (or a targeted list of interest). The data frame should contain three columns: col1= "complex_id", col2= "complex_name", and col3= "protein_id".
#' @param experiment_data A *data.frame* with your experiment results. 
#' @param N_fractions The number of fractions obtained in the co-fractionation experiment.
#' @param organism A string that indicates the specie. It could be "mmusculus", "hsapiens". Check gconvert vignette for more options https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
#' @param method_cor A string that indicates the correlation algorithm to make the matrix correlation. It can be "kendall", "spearman" or "pearson" (default "pearson")
#' @param heatmap_seaborn a logical value indicating whether to perform a correlation matrix compatible with a network plot from the corrr package (default TRUE).
#'
#' @return a list of potential protein complexes in the experiment which has two elements inside the data of the protein complex detected and its correlation matrix. If  heatmap_seaborn= TRUE it also includes a correlataion matrix made by the corrr package funcion correlate() 
#' @export
#' @import gprofiler2
#' @import gdata
#' @import corrr
#' @import dplyr
#' @import assertthat
#' @examples 
#' # Read the experiment files
#'data(Hek293_P2_1)
#'#### RUN mCP list that creates a list of potential protein complexes
#'CL_hek_P2_1<- mcp_list(corum_database =  Corum_Humans_Database,
#'                       experiment_data = Hek293_P2_1, 
#'                       N_fractions = 35, 
#'                       specie = "hsapiens",
#'                       method_cor = "pearson",
#'                       heatmap_seaborn = TRUE)
#' 
#' 
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
  
  
  max_attempts <- 5  # Number of maximum attempts
  current_attempt <- 1  # Current attempt count
  
  while (current_attempt <= max_attempts) {
    prot_names <- try({gprofiler2::gconvert(
      query = df$protein_id,
      organism =specie ,
      target = "ENSG",
      mthreshold = 1,
      filter_na = FALSE
    )}, silent = TRUE)  
    
    # Check if the function call was successful
    if (!inherits(prot_names, "try-error")) {
      # Function call was successful, break out of the loop
      break
    }
    
    # Function call failed, increment the attempt count
    current_attempt <- current_attempt + 1
    
    # Wait for some time before the next attempt
    Sys.sleep(5)  # Adjust the sleep time as per your requirement
  }
  # Check if the maximum number of attempts was reached
  if (current_attempt > max_attempts) {
    # Handle the case when the function call consistently failed
    print("Function call failed after multiple attempts.")
  } else {
    # Function call was successful
    print("Function call succeeded.")
  } 
  prot_names<- prot_names %>%
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