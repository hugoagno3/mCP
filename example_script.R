library(readr)
New_Corum_Mouse_09032022 <- read_csv("~/Documents/m_CP_Analysis/Documents to do DATA Analysis with m_CP/A single Profile fraccionation experiment/New_Corum_Mouse_09032022.csv")
pepQuant_Dig_C12E8_Triton_NP_40_Corrected <- read_csv("~/Documents/m_CP_Analysis/pP16-25_TEMF_Detergents/P16-25_DETERGENT_EMF_TEST/pepQuant_Dig_C12E8_Triton_NP-40_Corrected.csv")
View(pepQuant_Dig_C12E8_Triton_NP_40_Corrected)
str(pepQuant_Dig_C12E8_Triton_NP_40_Corrected)

devtools::load_all()
# make complexes list
library(dplyr)
Triton<- pepQuant_Dig_C12E8_Triton_NP_40_Corrected %>% select(1,36:137)  
CL1 <- mcp_list(corum_database =  New_Corum_Mouse_09032022,
               experiment_data = Triton)
##debug(cpp_plotter)
out_triton <- cpp_plotter(CL1, display_weights = TRUE, format = "", output_name = "Test1triron", filter = 0.93)
getwd()
out$`Sarcoglycan-sarcospan-syntrophin-dystrobrevin complex`
library(ggplot2)
d <- list()
for (i in 1:10) {
  g <- ggplot2::ggplot(data.frame(x = 1:10, y = i:i + 10)) + 
    ggplot2::geom_line(aes(x = x, y= y))
  d <- c(d, list(g))
  names(d)[i] <- i
}
