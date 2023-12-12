#' quant_mCP function: Analysis of total area of protein complexes extended to a low number of fractions. 
#'
#' @description This function is an integrated function of mCP package, that needs as input a list of mCP detected protein complexes (output of mCP function). It returns total area of the whole chromatogram of the average detected protein complex (With and Without correction), a list of average profile plots, and a plot with baseline correction and without baseline correction (by using Baseline R package). In addition, it plots 2 files as  outputs: 
#' 1- a .pdf file with average protein complexes profiles detected by mCP from a database.
#' 2- a .csv file containing 3 colums, complex name, area detected with corrections, and area detected without correction. 
#'   
#' @param list_ppc a List of protein complexes output of mCP function (generated with mCP function). 
#' @param pdf_name main core name for the pdf file to be printed, it is a string: see example (coud be "name.pdf")
#' @param area_name_file main core name for the csv file to be printed, it is a string: see example (coud be "name.csv") it is important to add .csv at the end.
#' @param method so it will always be 'rollingBall'. We recommend to use this function together with Baseline R package in chase to try other integration you will need to use a function that requires 2 parameters. If not, just create a customized function by using this as starting point is possible.  
#' @param wm   is width of local window for smoothing. We recommend consulting the baseline R package. The recommended range of this value is between 1-20. In particular as starting 1/10*number of fractions could be used. 
#' @param ws  is width of local window for minimization/maximization. We recommend consulting the baseline R package. The recommended range of this value is between 1-20. In particular as starting 1/10*number of fractions could be used.   
#'
#' @return a list of protein complexes plots, with 4 elements, The area corrected, the area non corrected, a plot of the Profile average of the protein complex and a plot of the Baseline R package with and without correction.  
#' @export
#' @import baseline
#' @import ggplot2
#' @import dplyr
#' @import pracma
#' @import assertthat
#' @examples 
#' 
#' 
#' 
#' Co-fractionation experiments with NP-40 mild-detergent_Humans(Hek293 cells)
#' 
#'data("Corum_Humans_Database")
#'data("Hek293_P2_1")
#'mCP_Hek_P2_1 <- mCP(corum_database = Corum_Humans_Database,
#'                                      experiment_data = Hek293_P2_1, 
#'                                      N_fractions = 35, 
#'                                      specie = "hsapiens",
#'                                      method_cor = "pearson",
#'                                      network = TRUE,
#'                                      format = "pdf", 
#'                                      output_name = "mCP_Hek293_P2_1",
#'                                      filter = 0.81,
#'                                      heat_map = TRUE,
#'                                      relative = FALSE,
#'                                      fdr_limit = 0.05,
#'                                      n_simulations= 7,
#'                                      monomeric_filter = FALSE,
#'                                      mw = TRUE,
#'                                      Risk_fraction = 31,
#'                                      set_seed = FALSE,
#'                                      display_weights = TRUE,
#'                                      standard_weights = list(list(x =6, label= "2700 KDa"), 
#'                                                         list(x = 11, label ="950 KDa"),
#'                                                         list(x = 14, label = "750 KDa"), 
#'                                                         list(x =27, label ="146 KDa"),
#'                                                         list(x =30, label ="60 KDa")))
#' Relative_area_PPC_mCP_Hek293_P2_1 <- Quant_mCP(list_ppc = mCP_Hek_P2_1,
#'                                        pdf_name = "relative_quant.pdf", 
#'                        area_name_file = "Hek293_P2_1_Relative_quant.csv", 
#'                                     wm= 4,
#'                                     ws= 4)                                                             
#'
#'data("Corum_Mouse_Database") ### Example with Mouse
#'data("CM_LV_1")
#'#### Important, mCP parameter relative = FALSE if that is not false, the function will not work. We recommend WS and WM should be 1/10 of the number of fractions
#'mCP_TEND_out_WT_Mouse_Cardiac <- mCP(corum_database = Corum_Mouse_Database,
#'                                   experiment_data = CM_LV_1, 
#'                                     N_fractions = 35, 
#'                                     specie = "mmusculus",
#'                                     method_cor = "pearson",
#'                                     network = TRUE,
#'                                     format = "pdf", 
#'                                     output_name = "INTO_m_CP_mouse_analysis_8",
#'                                     filter = 0.81,
#'                                     heat_map = TRUE,
#'                                     relative = FALSE,
#'                                     fdr_limit = 0.05,
#'                                     n_simulations= 14,
#'                                     monomeric_filter = FALSE,
#'                                     Risk_fraction = 31,
#'                                     set_seed = TRUE,
#'                                     mw = TRUE
#'                                     display_weights = TRUE,
#'                                     standard_weights =  list( 
#'                                     list(x = 9, label = "1048 KDa"), 
#'                                     list(x = 11.5, label = "720 KDa"),
#'                                     list(x = 16, label = "480 KDa"),
#'                                     list(x = 21, label = "242 KDa"),
#'                                     list(x = 25, label = "242 KDa"),
#'                                     list(x = 29, label = "166 KDa")))
#'                                     
#' Relative_area_PPC_mCP_cardiomyocytes_LV_1 <- Quant_mCP(list_ppc = mCP_TEND_out_WT_Mouse_Cardiac,
#'                                        pdf_name = "relative_quant_LV_WT_MOUSE.pdf", 
#'                        area_name_file = "LV_WT_MOUSE_Relative_quant.csv", 
#'                                     wm= 4,
#'                                     ws= 4)                                                             
#'
#'                                     
#'                                     
#'                                     

quant_mCP <- function(list_ppc= mCP_Hek_P2_1, pdf_name = "output_plot.pdf", area_name_file = "Relative_total_abundance.csv",
                      method= 'rollingBall', wm=4, ws=4) {
  result_list <- list()
  
  # Open PDF file in the working directory
  pdf(file = pdf_name, width = 10, height = 6)
  
  for (name in names(list_ppc)) {
    # Check if data is NULL or empty
    data <- list_ppc[[name]][[1]][["data"]]
    if (is.null(data) || nrow(data) == 0) {
      warning(paste("Skipping", name, "as data is NULL or empty"))
      next
    }
    
    # Check if data frame has columns 5 and 6
    if (!(5 %in% seq_along(data)) || !(6 %in% seq_along(data))) {
      warning(paste("Skipping", name, "as data frame does not have columns 5 and 6"))
      next
    }
    
    Complex539 <- data %>% dplyr::select(5, 6)
    colnames(Complex539)[1] <- "x"
    colnames(Complex539)[2] <- "y"
    Complex540 <- Complex539 %>% group_by(x) %>% summarise(y = mean(y))
    Complex540$x <- as.numeric(Complex540$x)
    Complex540$y <- as.numeric(Complex540$y)
    
    PL <- Complex540 %>% 
      ggplot(aes(x = x, y = y)) +
      geom_line() +
      ggtitle(name) +
      geom_point(color = "red", size = 1.5, show.legend = TRUE) +
      theme_minimal() +
      xlab("Fractions") +
      ylab("Intensity")
    # Convert y column to a matrix
    spectra_matrix <- t(matrix(Complex540$y, ncol = 1))
    
    # Baseline correction using the baseline package with "irls" method
    baseline_corrected <- baseline(spectra_matrix, wm=wm, ws=ws,
                                   method=method)
    Complex540 <- data.frame(x = Complex539$x, y_corrected = baseline_corrected@corrected[, 1])
    
    # Group by x and summarize y
    Complex541 <- Complex540 %>% group_by(x) %>% summarise(y = mean(y_corrected))
    
    # Create ggplot with baseline
    # PL <- ggplot() +
    #   geom_line(data = Complex541, aes(x = x, y = y), color = "blue", linetype = "solid", size = 1) +
    #   geom_line(data = Complex540, aes(x = x, y = y_corrected), color = "black", linetype = "solid", size = 1) +
    #   ggtitle(name) +
    #   geom_point(data = Complex541, aes(x = x, y = y), color = "red", size = 1.5, show.legend = TRUE) +
    #   theme_minimal() +
    #   xlab("Fractions") +
    #   ylab("Intensity")
    
    # Print the plot
    print(PL)
    plot(baseline_corrected)
    print(name)
    # Calculate area
    area_c <- trapz(Complex541$x, Complex541$y)
    Complex540_nc <- Complex539 %>% group_by(x) %>% summarise(y = mean(y))
    area_nc <- trapz(Complex540_nc$x, Complex540_nc$y)
    result_list[[name]] <- list(area_corrected = area_c, area_no_corrected= area_nc, plot1= PL, plot = baseline_corrected)
  }
  
  # Close the PDF device after all plots have been printed
  dev.off()
  
  # Extract names and areas
  names_list <- names(result_list)
  areas_list_corrected <- sapply(result_list, function(z) z$area_corrected)
  areas_list_non_corrected <- sapply(result_list, function(z) z$area_no_corrected)
  # Combine into a data frame
  results_df <- data.frame(names = names_list, area_corrected = areas_list_corrected, area_non_corrected= areas_list_non_corrected)
  
  # Write the data frame to a CSV file
  write.csv(results_df, file = area_name_file, row.names = FALSE)
  
  return(result_list)
}
