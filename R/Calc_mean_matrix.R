#' Calc_mean_matrix makes average of matrixes as input for mCP
#'
#' @description The function Calc_mean_matrix() recognizes the replicates by a tag like _A_ and _B_ as part of the names for the MS file and the *Frac_index*. It is the position in which you can find the number of fraction in the column name separated by underscores. for example the frac_index=5 is date_Surname_Measurement_A_01 because the number for the fraction is 5 spaces between underscores, like a_b_c_d_FractionNumber. 
#' @param NAmatrix is a Matrix input. It is an experimental matrix wide format first column protein_id and the other columns fractions with replicates 1 and 2. The matrix should be imported into R te function read.table("my_experimental_matrix.csv",sep =",",dec = ".", header= T)
#' @param pattern_group_A String that is use to differenciate replicates in this case _A_ is Replicate 1.  The function calc mean matrix will recognize the replicates by a tag like _A_ and _B_ as part of the names for the MS file 
#' @param pattern_group_B String that is use to differenciate replicates in this case _B_ is Replicate 2.
#' @param Protein_ID_column a numeric value that indicates where is the columns that contains the protein Ids.  
#' @param frac_index a numeric value that indicates the position in which you can find the number of fraction in the column name separated by underscores. for example the frac_index=5 is date_Surname_Measurement_A_01 because the number for the fraction is 5 spaces between underscores, like a_b_c_d_FractionNumber. 
#' @param save_file A logical TRUE or FALSE variant to save the file.
#' @param save_name an string that indicate the name of the file to be saved.
#'
#' @return Matrix with first column protein_id and the other colums the average of the fractions. 
#' @export
#'
#' @examples
#' data(NAmatrix_P2_1)
#' ### columns looks like these 
#' ###
#'[1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_A_01_GA1_1_6588.d"
#'[1] "D:\\Current Projects\\Projects_2022\\2022_77_H_Amede_mCP_Hek293_Cells\\03_Raw Data timsTOF Pro\\P2\\H_Amedei_04112022_P31_HEK293_P2_1_B_01_GA1_1_6626.d"
#'####
#'Hek293_input_for_mCP <- Calc_mean_matrix(NAmatrix = NAmatrix_P2_1,
#'                        pattern_group_A = "_A_",
#'                        pattern_group_B = "_B_",
#'                        frac_index=17,
#'                        Protein_ID_column = 1,
#'                        save_file = TRUE,
#'                        save_name = "Hek293_P2_1.csv")
#' 
#' 


Calc_mean_matrix <- function(NAmatrix,
                             pattern_group_A = "_A_",
                             pattern_group_B = "_B_",
                             Protein_ID_column=1,
                             frac_index=7,
                             save_file = FALSE,
                             save_name = "protpepQuant.csv"){
  #### data import ###############################################################
  #NAmatrix <- read.table(Input_path,sep =",",dec = ".", header= T)
  #### work flow #################################################################
  # select first Protein Id
  NAmatrix[,Protein_ID_column] <- sapply(strsplit(NAmatrix[,Protein_ID_column],split = ";"),"[[",1)
  
  # extract Identifiers
  Peptid_ID_table <- NAmatrix[,1:2]
  
  # split Matrix in a und B
  Mat_A <- as.matrix(split_mat(NAmatrix[,3:ncol(NAmatrix)],pattern = pattern_group_A,Column = T))
  Mat_A[is.na(Mat_A)] <- 0
  Mat_B <- as.matrix(split_mat(NAmatrix[,3:ncol(NAmatrix)],pattern = pattern_group_B,Column = T))
  Mat_B[is.na(Mat_B)] <- 0

  Fracs_A <- as.numeric(sapply(strsplit(sapply(strsplit(colnames(Mat_A),split="_"),"[[",frac_index),split=".",fixed=T),"[[",1))
  Fracs_B <- as.numeric(sapply(strsplit(sapply(strsplit(colnames(Mat_B),split="_"),"[[",frac_index),split=".",fixed=T),"[[",1))
  
  Fracs <- unique(c(Fracs_A,Fracs_B))
  
  # create new Dataframe with all Fracions:
  Fractions <- unique(c(Fracs_A,Fracs_B))
  Mean_Mat <- matrix(numeric(length(Fractions)*nrow(NAmatrix)),nrow = nrow(NAmatrix),ncol=length(Fractions),
                     dimnames = list(Peptid_ID_table$PEP.GroupingKey,Fractions))
  
  # Add up values in Mean_mat
  Indeces <- match(Fracs_A,colnames(Mean_Mat))
  Mean_Mat[,Indeces] <- Mean_Mat[,Indeces] + Mat_A
  
  Indeces <- match(Fracs_B,colnames(Mean_Mat))
  Mean_Mat[,Indeces] <- Mean_Mat[,Indeces] + Mat_B
  
  Indeces <- match(na.omit(Fracs_A[match(Fracs_B,Fracs_A)]),colnames(Mean_Mat))
  Mean_Mat[,Indeces] <- Mean_Mat[,Indeces]/2
  
  # sort Fractions in Mean Matrix 
  Mean_Mat <- Mean_Mat[,match(sort(as.numeric(colnames(Mean_Mat))),colnames(Mean_Mat))]
  
  # combine peptide rows to protein rows
  Protein_mat <- combine_1_n_mapping(Peptid_ID_table[,Protein_ID_column],Mean_Mat)
  
  # plot if neccessary:
  # heatmap((G[1:300,1:300]>0.70)*1,scale = "none",keep.dendro = FALSE)
  Protein_mat <- data.frame(protein_id=rownames(Protein_mat),as.data.frame(Protein_mat))
  # save resultmatrix
  if(save_file){write.csv(Protein_mat,file = save_name, row.names = FALSE)}
  
  # return Matrix
  return(Protein_mat)
}

# Example
# NAmatrix <- read.table("/home/hugo/Documents/m_CP_Analysis/P29/283vs282/P29_282_WT_283.csv",sep =",",dec = ".", header= T)
# getwd()
# Mouse_1 <- NAmatrix[,1:179]
# WT_282 <- Calc_mean_matrix(NAmatrix = Mouse_1,
#                      pattern_group_A = "_A_",
#                      pattern_group_B = "_B_",
#                      frac_index=7,
#                      Protein_ID_column = 1,
#                      save_file = FALSE,
#                      save_name = "protpepQuant.csv")
# View(WT_282)
#colnames(NAmatrix)<-gsub("__","_",colnames(NAmatrix))



# ToDo save with protein_id as first column
#write.csv(A,"P29_282_WT_283_Niels.csv",row.names = TRUE)
#write.table(NAmatrix,"P29_282_WT_283.csv",sep=",",row.names = FALSE)

