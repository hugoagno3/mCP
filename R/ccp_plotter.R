#'
#' cpp_plotter
#'
#' @description The cpp_plotter() function takes a list of potential protein complexes as input. The function first filters out any complexes composed of less than 2 components with non-zero intensities. Then, it applies a filter based on the minimum definition of a protein complex, keeping only complexes with at least one significant co-elution hit. The function also filters the list of complexes based on the presence of co-eluting proteins, requiring at least one binary interaction higher than the filter for proteins within the candidate complex. Finally, cpp_plotter() generates complexome profiling plots, heatmaps, and network plots of the proteins within the selected complexes.
#' @param complex_list a list output from mcp_list function. It contains all potential protein complexes presents in the experiment and CORUM database. Each element of the list is has its corresponding complex name.
#' @param N_fractions The number of fractions obtained in the co-fractionation experiment.
#' @param filter targeted value of minimum accepted binary interaction hits between 2 proteins within a protein complex (we recomend a 0.93 value).
#' @param output_name main core name for the files to be printed, it is a string: see example. 
#' @param format should be "pdf" to plots pdf. There is an additional possible value is ".": this will omit pdf. In both cases the user will have a list of plot for the detected protein complexes in R. 
#' @param relative could be = FALSE or TRUE. When TRUE the plots will be relative to the maximum intensity of each protiens. When FALSE, the plot will be the intensities detected in the experiment without relativization. This option is useful when some of the proteins in a protein protein complex are not easy to detect in the Mass spectrometer
#' @param heat_map could be = FALSE or TRUE. When TRUE it will plot heat maps of all protein complexes detected after filter them by Pearson correlation. When the algorithm of correlation is not Pearson is should be always FALSE. 
#' @param network could be = FALSE or TRUE. When TRUE it will plot networks heatmaps of all protein complexes detected after filter them by Pearson correlation. 
#' @param display_weights TRUE or FALSE statement to display molecular weight markers in the complexome profiling plots. 
#' @param standard_weights following core= list(list(x =11, label= "1049KDa"), list(x = 13, label ="720 KDa"))). It is possible to add many markers. you have to extend the code, for example to add a thrid marker= list(list(x =11, label= "1049KDa"), list(x=12, label="900 kDa"), list(x = 13, label ="720 KDa"))). Display_weights muss be TRUE. 
#' @param mw true false statement for mw conversion (requires std weights in kDa)
#'
#' @return A list of protein complexes filter by pearson correlation, a pdf file with the detected as protein complexes profiles. A pdf with heatmaps of the detected protein complexes. A txt file with numbers about general false positive when atleast 1 hit is consider as filter. A CVS file containing all protein complexes detected, hits of binary interactions inside the protein complexes, FDR detected by MonteCarloSimulation.
#' @import gprofiler2
#' @import ggplot2
#' @import gdata
#' @import corrr
#' @import dplyr
#' @import assertthat
#' @import  readr
#'  
#' @export
#'
#' @examples For the exaple load the 2 datasets and run mcp_list function and then cpp_ploter 
#'
#'#'out_Hek_P2_1 <- cpp_plotter(complex_list = CL_hek_P2_1,
#'                            format = "pdf", 
#'                            output_name = "m_CP_analysis",
#'                            filter = 0.81,
#'                            N_fractions = 35,
#'                            heat_map = TRUE,
#'                            mw = TRUE,
#'                            relative = FALSE,
#'                            display_weights = TRUE,
#'                            standard_weights = list(list(x =6, label= "2700 KDa"), 
#'                                               list(x = 11, label ="950 KDa"),
#'                                               list(x = 14, label = "750 KDa"), 
#'                                               list(x =27, label ="146 KDa"),
#'                                               list(x =30, label ="60 KDa")))
#'                                               
#'### To generate an example
#'data(Hek293_P2_1)
#'
#'data(Corum_Humans_Database) 
#'
#'CL_hek_P2_1<- mcp_list(corum_database =  Corum_Humans_Database,
#'experiment_data = Hek293_P2_1, 
#'N_fractions = 35, 
#'specie = "hsapiens",
#'method_cor = "pearson",
#'network = TRUE)
#'
#'##### Run the output of mcp_list into the cpp_ploter function. 
#'out_Hek_P2_1 <- cpp_plotter(complex_list = CL_hek_P2_1,
#'                            format = "pdf", 
#'                            output_name = "m_CP_analysis",
#'                            filter = 0.81,
#'                            N_fractions = 35,
#'                            heat_map = TRUE,
#'                            relative = FALSE,
#'                            mw = TRUE,
#'                            display_weights = TRUE,
#'                            standard_weights = list(list(x =6, label= "2700 KDa"), 
#'                                               list(x = 11, label ="950 KDa"),
#'                                               list(x = 14, label = "750 KDa"), 
#'                                               list(x =27, label ="146 KDa"),
#'                                               list(x =30, label ="60 KDa")))



cpp_plotter <- function (relative= FALSE, heat_map= FALSE, network= FALSE,
                         complex_list, N_fractions = 35, filter = 0.81,
                         output_name = paste0("complexes_detected_", 
                         Sys.Date()), format = "pdf", mw= FALSE, display_weights = FALSE,
                         standard_weights = list(list(x = 9, label = "1236 KDa"),
                                                 list(x = 13, label = "720 KDa"))) 
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
  if ((mw & !is.list(standard_weights)) | (mw & length(standard_weights)<2)){
    stop("standard_weights parameter missed or you have less than 2 points")
  }
  
  if (tolower(format) == "pdf") {
    pdf(file = paste0(output_name, ".", format), width = 8, 
        height = 6)
  }
  for (i in seq_along(complex_list)) {
    data <- complex_list[[i]]$data
    corMat <- complex_list[[i]]$corMat
    ########### interactomics
    if(sum(dim(corMat))==0 | nrow(data) < N_fractions | all(data$Intensity[which(!is.na(data$Intensity))] == 0)){
      next
    }
    pectores<- corMat
    pectores<- gdata::unmatrix(pectores)
    pere<- as.data.frame(pectores[pectores>filter & pectores<1])
    colnames(pere)[1]<- "Cor"
    pere_1<- pere %>% dplyr::distinct(Cor, .keep_all = TRUE) 
    pere_1$pairs<- rownames(pere_1)
    perek_2<- list(interactions=pere_1$pairs, Cor= pere_1$Cor)
    ###########
    tri <- corMat[upper.tri(corMat)]
    tri[is.na(tri)] <- 0
    # if (nrow(data) > N_fractions & !all(data["Intensity"] == 
    #                                     0)) {
    if (any(tri > filter)) {
      vere<- data.frame(Hits= length(which(tri > filter)))
      c_counter <- c_counter + 1
    if (relative) {data <- data %>%group_by(protein_id) %>%
        mutate(Intensity = Intensity/max(Intensity, na.rm=TRUE))
      }
    if (mw){ # if plotting with molecular weight is prefered
        x_values = sapply(standard_weights,"[[","x")
        y_values = sapply(standard_weights,"[[","label")
        y_values = log10(parse_number(y_values))
        mw_model = lm(y_values~x_values)
 
        p <- ggplot2::ggplot(data, ggplot2::aes(x = SEC_FR , 
                                                y = Intensity, ggplot2::ggtitle(complex_name), 
                                                col = prot_name)) + ggplot2::geom_line() + 
          ggplot2::scale_x_continuous(name = "Molecular Weight [kDa]", labels = round(10**(mw_model$coefficients["x_values"] * seq(1, N_fractions, 5) + mw_model$coefficients["(Intercept)"])),
                                      breaks = seq(1, N_fractions, 5))+ 
          ggplot2::ggtitle(data$complex_name) + 
          ggplot2::theme_light() + ggplot2::theme(axis.text.x = element_text(angle = 90))
      }
      
      else{ #if ploting with fraction labeling is prefered
        p <- ggplot2::ggplot(data, ggplot2::aes(x = SEC_FR, 
                                                y = Intensity, ggplot2::ggtitle(complex_name), 
                                                col = prot_name)) + ggplot2::geom_line() + 
          ggplot2::scale_x_continuous(name = "Fractions", 
                                      breaks = seq(1, N_fractions, 5)) + ggplot2::ggtitle(data$complex_name) + 
          ggplot2::theme_light()  
        if (display_weights) {
          for (weight in standard_weights) {
            p <- p + ggplot2::geom_vline(xintercept = weight$x, 
                                         colour = "grey", linetype = "dashed") + 
              ggplot2::annotate("text", x = weight$x - 
                                  0.5, y = mean(data$Intensity), label = weight$label, 
                                angle = 90, color = "grey")
          }
        }
      }
      
      plots_list <- c(plots_list, list(c(list(p),list(vere),list(perek_2))))
      cp_names <- c(cp_names, as.character(data$complex_name[1]))
      if (tolower(format) == "pdf") {
        print(p)
      }
    }
    #}
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
    
    for (i in seq_along(complex_list)) {
      data <- complex_list[[i]]$data
      corMat <- complex_list[[i]]$corMat
      if(sum(dim(corMat))==0 | nrow(data) < N_fractions | all(data$Intensity[which(!is.na(data$Intensity))] == 0)){
        next
      } 
      
      tri <- corMat[upper.tri(corMat)]
      tri[is.na(tri)] <- 0
      # if (nrow(data) > N_fractions & !all(data["Intensity"] == 
      #                                     0)) {
      if (any(tri > filter)) {
        corMat[is.na(corMat)] <- 0 ##
        g <-heatmap(corMat, scale = "none", main = unique(data$complex_name))
        if (network){
          complex_list[[i]][["CorMat_rrr"]][is.na(complex_list[[i]][["CorMat_rrr"]])]<-0 ##
          
          if(sum(apply(complex_list[[i]][["CorMat_rrr"]][, -1], 1, duplicated))==0){
            d<- try( corrr::network_plot (complex_list[[i]][["CorMat_rrr"]], min_cor = 0.3))  
          }
          else{
            colname<- colnames(complex_list[[i]][["CorMat_rrr"]])[1]
            d<-  as.data.frame(complex_list[[i]][["CorMat_rrr"]]) %>% ggplot(aes(y={{colname}}))+ geom_bar()+ ggtitle("Error on the matrix, please try manually")
          }
          
        }
        cp_names <- c(cp_names, as.character(data$complex_name[1]))
        if (tolower(format) == "pdf" & heat_map) {
          print(g)
          if(network){
            
            print(d)
          }
        }
        if (network){
          plots_list_heatmaps <- c(plots_list_heatmaps, list(c(list(d))))
        }
      }
      
    }
    dev.off()
    if(network){
      names(plots_list_heatmaps) <- names(plots_list)
      plots_list <- Map(c, plots_list, plots_list_heatmaps)}
  }
  
  ##### End of heatmap 
  print(paste0(c_counter, " Complexes were detected"))
  return(plots_list)
}
