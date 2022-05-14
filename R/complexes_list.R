# corum_database <- New_Corum_Mouse_09032022
# experiment_data <- pepQuant_Dig_C12E8_Triton_NP_40_Corrected
# View(corum_database)
# View(experiment_data)
#' List of Protein Complexes
#' 
#' @description This function detects a list of potential list of complexes... 
#' @param corum_database 
#' @param experiment_data 
#' @param N_fractions 
#' @param condition 
#'
#' @return
#' @export
#'
#' @examples
mcp_list <- function(corum_database, experiment_data, N_fractions = 34, condition) {
  
  # datacleaning
  # - check user input
  # - NA, remove duplicates
  # join dataframes by protein_id
  
  #assertthat::assert_that(class(corum_database) == "data.frame")
  # ... 
  Matrix_Clast1 <-
    experiment_data %>% distinct(protein_id, .keep_all = TRUE)
  Viaexp <- Matrix_Clast1 %>% select(1, 'Pi1...2':'Pi34...35')
  df <- dplyr::inner_join(corum_database, Viaexp, by = 'protein_id')
  
  prot_names <- gprofiler2::gconvert(
    query = df$protein_id,
    organism = "mmusculus",
    target = "ENSG",
    mthreshold = 1,
    filter_na = FALSE
  ) %>% 
    select(name, input) %>% #TODO: include description?
    mutate(name = ifelse(is.na(name), yes = input, no = name)) %>% # replace NA
    rename(protein_id = input)
  
  assertthat::are_equal(prot_names$protein_id, df$protein_id)
  
  df1 <- df %>% tibble::add_column(prot_name = prot_names$name, .before = 3)
  df_long <- df1 %>%
    tidyr::gather("SEC_FR", "Intensity", 'Pi1...2':'Pi34...35') %>% 
    dplyr::arrange(complex_id)
  
  return(df_long)
}


