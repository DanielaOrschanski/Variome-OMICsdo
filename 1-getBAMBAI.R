#path_dir = path indicando el directorio donde están R1 y R2 (archivos .fastq)

getBAMBAI <- function(path_dir, ref = "HG19") {
  if (dir.exists(paste(path_dir, "/OMICsdo", sep = ""))) {
    file_list <- list.files(paste(path_dir, "/OMICsdo", sep=""))
    bam_file <- file_list[endsWith(file_list, "sortedR1R2.bam")]
    
    if (!(length(nchar(bam_file)) == 0)) {
      bam.control <- sprintf("%s/OMICsdo/%s", path_dir, bam_file)
    } else {
      bam.control <- fasta_to_bam(path_dir, ref = ref)
    }
    bai_file <- file_list[endsWith(file_list, "sortedR1R2.bam.bai")]
    if (length(nchar(bai_file)) == 0) { indexBam(bam.control) }
    
  } else {
    bam.control <- fasta_to_bam(path_dir, ref = ref)
    indexBam(bam.control)
  }
  
  return(bam.control)
}
