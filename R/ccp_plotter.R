

# library(readr)
# corum_database <- read_csv("~/Downloads/New_Corum_Mouse_09032022.csv")
# experiment_data <- read_csv("~/Downloads/pepQuant_Dig_C12E8_Triton_NP-40_Corrected.csv")
# 
# t0 <- Sys.time()
# complexes_list <- mcp_list(corum_database, experiment_data)
# print(Sys.time() - t0)



#' Title
#'
#' @param complex_list 
#' @param N_fractions 
#' @param filter 
#' @param output_name 
#' @param format 
#' @param display_weights 
#' @param standard_weights 
#'
#' @return
#' @export
#'
#' @examples
cpp_plotter <- function(
  complex_list,
  N_fractions = 35,
  filter = 0.93,
  output_name = paste0("complexes_detected_", Sys.Date()),
  format = "pdf", 
  display_weights = FALSE,
  standard_weights = list(
    list(x = 9, label = "1236 KDa"), 
    list(x = 13, label = "720 KDa")
  )
) {
  
  # Input validation
  assertthat::assert_that(is.list(standard_weights))
  assertthat::assert_that(all(sapply(standard_weights, function(standard)
    length(standard) == 2)))
  assertthat::assert_that(all(sapply(standard_weights, function(standard)
    names(standard) == c("x", "label"))))
  assertthat::assert_that(all(sapply(standard_weights, function(standard)
    is.numeric(standard$x))))
  
  
  c_counter <- 0
  pdf(file = paste0(output_name, ".", format), width = 8, height = 6)
  for (i in seq_along(complex_list)) {
    
    data <- complex_list[[i]]$data
    corMat <- complex_list[[i]]$corMat
    tri <- corMat[upper.tri(corMat)]
    
    if (nrow(data) > N_fractions & !all(data["Intensity"] == 0)) {
      if (any(tri > filter)){
        c_counter <- c_counter + 1        
        
        # plot
        p <- ggplot2::ggplot(data,
                             ggplot2::aes(
                               x = SEC_FR,
                               y = Intensity,
                               ggplot2::ggtitle(complex_name),
                               col = prot_name
                             )) +
          ggplot2::geom_line() +
          ggplot2::scale_x_continuous(name = "Fractions", breaks = seq(1, N_fractions, 2)) +
          ggplot2::ggtitle(data$complex_name) +
          ggplot2::theme_minimal()
        
        if (display_weights) {
          standard_x <- sapply(standard_weights, function(standard) standard$x)
          standard_labels <- sapply(standard_weights, function(standard) standard$label)
            p <- p +
              ggplot2::geom_vline(xintercept = standard_x,
                                  colour = "grey",
                                  linetype = "dashed") +
              ggplot2::annotate("text", x = standard_x - 0.5, y = mean(data$Intensity),
                                label = standard_labels, angle = 90, 
                                color = "grey",)
              
          }
        }
        
        complex_list[[i]]$p <- p
        print(p)
      }
  }
  dev.off()
  
  print(paste0(c_counter, " Complexes were detected"))
  return(p)
}
