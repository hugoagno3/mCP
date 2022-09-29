#' Calc_mean_matrix
#'
#' @param NAmatrix 
#' @param pattern_group_A 
#' @param pattern_group_B 
#' @param Protein_ID_column 
#' @param frac_index 
#' @param save_file 
#' @param save_name 
#'
#' @return
#' @export
#'
#' @examples

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

