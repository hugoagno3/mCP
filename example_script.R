library(readr)
New_Corum_Mouse_09032022 <- read_csv("~/Downloads/New_Corum_Mouse_09032022.csv")
pepQuant_Dig_C12E8_Triton_NP_40_Corrected <- read_csv("~/Downloads/pepQuant_Dig_C12E8_Triton_NP-40_Corrected.csv")

str(pepQuant_Dig_C12E8_Triton_NP_40_Corrected)


# make complexes list
CL <- mcp_list(corum_database =  New_Corum_Mouse_09032022,
               experiment_data = pepQuant_Dig_C12E8_Triton_NP_40_Corrected)
out <- cpp_plotter(CL)
