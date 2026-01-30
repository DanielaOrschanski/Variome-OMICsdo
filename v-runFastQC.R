#' @title Run FastQC
#' @description Executes the FastQC for R1 and R2.
#' @param patent_dir Path of the directory that contains R1 and R2 fastq files. It can be either the "trimmed" folder or the original folder.
#' @import stringr
#' @import viridis
#' @import reshape2
#' @import readr
#' @export
runFastQC <- function(patient_dir) {
  
  file_list <- list.files(patient_dir)
  
  #Para que se pueda poner como entrada la carpeta de los trimmeados o la carpeta original:
  if  (startsWith(basename(patient_dir), "trimmed")) {
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "val_1.fq.gz")])) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("val_1.fq%s", gzip))], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("val_2.fq%s", gzip))], sep="")
  } else {
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))], sep="")
  }
  
  if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
    stop("There are no fastq files in this directory")
  }
  
  #Evita repetir el analisis si ya fue hecho
  if (startsWith(basename(patient_dir), "trimmed")) {
    if (!(length(nchar(file_list[endsWith(file_list, "val_1_fastqc")])) == 0) & !(length(nchar(file_list[endsWith(file_list, "val_2_fastqc")])) == 0)) {
      message("The FastQC for this sample has already been done.")
      return(paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1_fastqc")], sep=""))
    }
  } else {
    if ((length(nchar(file_list[endsWith(file_list, "R1_fastqc")])) != 0) & (length(nchar(file_list[endsWith(file_list, "R2_fastqc")])) != 0)) {
      message("The FastQC for this sample has already been done.")
      return(paste0(patient_dir, "/", file_list[endsWith(file_list, "R1_fastqc")], sep=""))
    }
  }
  
  #Ejecuta el FastQC R1
  system2(FastQC, fileR1)
  file_list <- list.files(patient_dir)
  
  #Extrae el zip R1
  file_fastqc_zip <- paste0(patient_dir, "/", file_list[endsWith(file_list, "1_fastqc.zip")], sep="")
  unzip(file_fastqc_zip, exdir = patient_dir)
  
  #Ejecuta el FastQC R2
  system2(FastQC, fileR2)
  file_list <- list.files(patient_dir)
  
  #Extrae el zip R2
  file_fastqc_zip <- paste0(patient_dir, "/", file_list[endsWith(file_list, "2_fastqc.zip")], sep="")
  unzip(file_fastqc_zip, exdir = patient_dir)
  
  message("FastQC's analysis has finished!")
  
}