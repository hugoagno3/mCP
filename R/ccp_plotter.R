

# library(readr)
# corum_database <- read_csv("~/Downloads/New_Corum_Mouse_09032022.csv")
# experiment_data <- read_csv("~/Downloads/pepQuant_Dig_C12E8_Triton_NP-40_Corrected.csv")
# 
# t0 <- Sys.time()
# complexes_list <- mcp_list(corum_database, experiment_data)
# print(Sys.time() - t0)



cpp_plotter <- function(complex_list, filter = 0.93, format = "pdf") {
  
  c_counter <- 0
  pdf(file = "complexes_detected_P16-25_miniBNE_35UG_WT_Digitonin_0.93.pdf", width = 8, height = 6)
  for (i in seq_along(complex_list)) {
    
    data <- complex_list[[i]]$data
    corMat <- complex_list[[i]]$corMat
    tri <- corMat[upper.tri(corMat)]
    
    if (nrow(data) > 37 & !all(data["Intensity"] == 0)) {
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
          ggplot2::scale_x_continuous(name = "Fractions", breaks = seq(1, 37, 2)) +
          ggplot2::ggtitle(data$complex_name) +
          ggplot2::geom_vline(xintercept = 9,
                              colour = "grey",
                              linetype = "dashed") +
          ggplot2::geom_vline(xintercept = 13,
                              colour = "grey",
                              linetype = "dashed") +
          ggplot2::geom_text(
            ggplot2::aes(x = 9 - 0.5, label = "1236 KDa" , y = 200000),
            angle = 90,
            text = ggplot2::element_text(size = 11)
          ) +
          ggplot2::geom_text(
            ggplot2::aes(x = 13 - 0.5, label = "720 KDa", y = 200000),
            angle = 90,
            text = ggplot2::element_text(size = 11)
          )
        
        complex_list[[i]]$p <- p
        print(p)
      }
    }
  }
  dev.off()
  
  print(paste0(c_counter, " Complexes were detected"))
  return(complex_list)
}