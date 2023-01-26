

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
#' @param relative could be = FALSE or TRUE. When TRUE the plots will be relative to the maximum intensity of each protiens. When FALSE, the plot will be the intensities detected in the experiment without relativization. This option is useful when some of the proteins in a protein protein complex are not easy to detect in the Mass spectrometer
#' @param heat_map could be = FALSE or TRUE. When TRUE it will plot heat maps of all protein complexes detected after filter them by Pearson correlation.
#' @param display_weights 
#' @param standard_weights 
#'
#' @return
#' @export
#'
#' @examples
cpp_plotter <- function (complex_list, N_fractions = 35, filter = 0.93, 
                         output_name = paste0("complexes_detected_", Sys.Date()), format = "pdf", 
                         relative= FALSE, heat_map= FALSE, display_weights = FALSE, 
                         standard_weights = list(list(x = 9, label = "1236 KDa"), list(x = 13, label = "720 KDa"))) 
{
  assertthat::assert_that(is.list(standard_weights))
  assertthat::assert_that(all(sapply(standard_weights, function(standard) length(standard) == 
                                       2)))
  assertthat::assert_that(all(sapply(standard_weights, function(standard) names(standard) == 
                                       c("x", "label"))))
  assertthat::assert_that(all(sapply(standard_weights, function(standard) is.numeric(standard$x))))
  c_counter <- 0
  plots_list <- list()
  cp_names <- list()
  while (!is.null(dev.list())) {
    dev.off()
  }
  if (tolower(format) == "pdf") {
    pdf(file = paste0(output_name, ".", format), width = 8, 
        height = 6)
  }
  for (i in seq_along(complex_list)) {
    data <- complex_list[[i]]$data
    corMat <- complex_list[[i]]$corMat
    tri <- corMat[upper.tri(corMat)]
    tri[is.na(tri)] <- 0
    if (nrow(data) > N_fractions & !all(data["Intensity"] == 
                                        0)) {
      if (any(tri > filter)) {
        c_counter <- c_counter + 1
        if (relative) {data <- data %>%group_by(protein_id) %>%
          mutate(Intensity = Intensity/max(Intensity, na.rm=TRUE))
        }
        p <- ggplot2::ggplot(data, ggplot2::aes(x = SEC_FR, 
                                                y = Intensity, ggplot2::ggtitle(complex_name), 
                                                col = prot_name)) + ggplot2::geom_line() + 
          ggplot2::scale_x_continuous(name = "Fractions", 
                                      breaks = seq(1, N_fractions, 5)) + ggplot2::ggtitle(data$complex_name) + 
          ggplot2::theme_minimal()
        if (display_weights) {
          for (weight in standard_weights) {
            p <- p + ggplot2::geom_vline(xintercept = weight$x, 
                                         colour = "grey", linetype = "dashed") + 
              ggplot2::annotate("text", x = weight$x - 
                                  0.5, y = mean(data$Intensity), label = weight$label, 
                                angle = 90, color = "grey")
          }
        }
        plots_list <- c(plots_list, list(p))
        cp_names <- c(cp_names, as.character(data$complex_name[1]))
        if (tolower(format) == "pdf") {
          print(p)
        }
      }
    }
  }
  if (tolower(format) == "pdf") {
    dev.off()
  }
  names(plots_list) <- cp_names
  ######## heatmaps part
  while (!is.null(dev.list())) {
    dev.off()
  }
  if (tolower(format) == "pdf" & heat_map) {
    pdf(file = paste0("heatmap_",output_name, ".", format), width = 8, 
        height = 6)
    plots_list_heatmaps <- list()
  }
  for (i in seq_along(complex_list)) {
    data <- complex_list[[i]]$data
    corMat <- complex_list[[i]]$corMat
    tri <- corMat[upper.tri(corMat)]
    tri[is.na(tri)] <- 0
    if (nrow(data) > N_fractions & !all(data["Intensity"] == 
                                        0)) {
      if (any(tri > filter)) {
        p <-heatmap(corMat, scale = "none", main = unique(data$complex_name))
        
        plots_list_heatmaps <- c(plots_list_heatmaps, list(p))
        cp_names <- c(cp_names, as.character(data$complex_name[1]))
        if (tolower(format) == "pdf" & heat_map) {
          print(p)
        }
      }
    }
  }
  if (tolower(format) == "pdf" & heat_map) {
    dev.off()
  }
  ##### End of heatmap 
  print(paste0(c_counter, " Complexes were detected"))
  return(plots_list)
}