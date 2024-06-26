% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Calc_mean_matrix.R
\name{Calc_mean_matrix}
\alias{Calc_mean_matrix}
\title{Calc_mean_matrix built the mCP matrix input by averaco of replicates fractions}
\usage{
Calc_mean_matrix(
  NAmatrix,
  pattern_group_A = "_A_",
  pattern_group_B = "_B_",
  Protein_ID_column = 1,
  frac_index = 7,
  save_file = FALSE,
  save_name = "protpepQuant.csv"
)
}
\arguments{
\item{NAmatrix}{is a Matrix input. It is an experimental matrix wide format first column protein_id and the other columns fractions with replicates 1 and 2. The matrix should be imported into R te function read.table("my_experimental_matrix.csv",sep =",",dec = ".", header= T)}

\item{pattern_group_A}{String that is use to differenciate replicates in this case \emph{A} is Replicate 1.  The function calc mean matrix will recognize the replicates by a tag like \emph{A} and \emph{B} as part of the names for the MS file}

\item{pattern_group_B}{String that is use to differenciate replicates in this case \emph{B} is Replicate 2.}

\item{Protein_ID_column}{a numeric value that indicates where is the columns that contains the protein Ids.}

\item{frac_index}{a numeric value that indicates the position in which you can find the number of fraction in the column name separated by underscores. for example the frac_index=5 is date_Surname_Measurement_A_01 because the number for the fraction is 5 spaces between underscores, like a_b_c_d_FractionNumber.}

\item{save_file}{A logical TRUE or FALSE variant to save the file.}

\item{save_name}{an string that indicate the name of the file to be saved.}
}
\value{
Matrix with first column protein_id and the other colums the average of the fractions.
}
\description{
The function Calc_mean_matrix() recognizes the replicates by a tag like \emph{A} and \emph{B} as part of the names for the MS file and the \emph{Frac_index}. It is the position in which you can find the number of fraction in the column name separated by underscores. for example the frac_index=5 is date_Surname_Measurement_A_01 because the number for the fraction is 5 spaces between underscores, like a_b_c_d_FractionNumber.
}
\examples{
data(NAmatrix_P2_1)
### columns looks like these 
###
[1] "D:\\\\Current Projects\\\\Projects_2022\\\\2022_77_H_Amede_mCP_Hek293_Cells\\\\03_Raw Data timsTOF Pro\\\\P2\\\\H_Amedei_04112022_P31_HEK293_P2_1_A_01_GA1_1_6588.d"
[1] "D:\\\\Current Projects\\\\Projects_2022\\\\2022_77_H_Amede_mCP_Hek293_Cells\\\\03_Raw Data timsTOF Pro\\\\P2\\\\H_Amedei_04112022_P31_HEK293_P2_1_B_01_GA1_1_6626.d"
#### _A_ and _B_ are patterns for recognition groups and the numeric position that separate the number of fraction is 17 between underscore. 
Hek293_input_for_mCP <- Calc_mean_matrix(NAmatrix = NAmatrix_P2_1,
                       pattern_group_A = "_A_",
                       pattern_group_B = "_B_",
                       frac_index=17,
                       Protein_ID_column = 1,
                       save_file = TRUE,
                       save_name = "Hek293_P2_1.csv")


}
