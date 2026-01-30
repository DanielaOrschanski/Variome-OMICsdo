#' @title Creates a VCF file from a BAM file
#' @description Selecting the path where the BAM file is saved, bam_to_vcf creates a VCF (Variant Call Format) file.
#' It uses a basic pipeline to create it, using the softwares: TrimGalore, BWA, PICARD, GATK, SamTools.
#' The default filters used for keeping or removing variants are:
#' For SNPs: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
#' For indels: 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
#' @return Returns the path of the created VCF file named "<patient>_filtered_MitoR.vcf".
#' Several files are created and saved in a new directory called "MitoR" at the same directory as the R1 and R2 reads were located:
#' 1 The VCF file
#' 2 The VCF index file. It is saved
#' @import magrittr
#' @import dplyr
#' @export
#' @examples
#' bam_to_vcf("~/Github/MitochondriaAnalysis/MitoR/Patients/37019.bam")

bam_to_vcf <- function(sorted_BAM, ref, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {
  
  path_to_MitoR <- substr(sorted_BAM, start=0, stop=(nchar(sorted_BAM)-nchar(basename(sorted_BAM))))
  patient_ID <- strsplit(path_to_MitoR, "/") %>% unlist()
  patient_ID <- patient_ID[length(patient_ID) - 1]
  
  if(file.exists(sprintf("%s%s_OMICsdo.vcf", path_to_MitoR, patient_ID))) {
    message("The VCF has already been generated for this patient")
    return(sprintf("%s%s_OMICsdo.vcf", path_to_MitoR, patient_ID))
  }
  
  #######################################3
  #Agregado por mi
  omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
 
  #Checks all the softwares needed are downloaded ----
  linea_software <- grep("(?i)Samtools", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    Samtools <- downloadSamtools()
  } else {
    Samtools <-  sub("^Samtools ", "", linea_software)
  }
  
  linea_software <- grep("(?i)BWA", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    BWA <- downloadBWA()
  } else {
    BWA <-  sub("^BWA ", "", linea_software)
  }
  
  linea_software <- grep("(?i)PICARD", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    PICARD <- downloadPICARD()
  } else {
    PICARD <-  sub("^PICARD ", "", linea_software)
  }
  
  linea_software <- grep("(?i)GATK", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    GATK <- downloadGATK()
  } else {
    GATK <-  sub("^GATK ", "", linea_software)
  }
  
  if (ref == "HG38") {
    linea_software <- grep("(?i)HG38FASTA", softwares, ignore.case = TRUE, value = TRUE)
    if(length(nchar(linea_software)) == 0) {
      reference <- downloadHG38()
    } else {
      reference <- sub("^HG38FASTA ", "", linea_software)
    }
  } else {
    linea_software <- grep("(?i)HG19FASTA", softwares, ignore.case = TRUE, value = TRUE)
    if(length(nchar(linea_software)) == 0) {
      reference <- downloadHG19()
    } else {
      reference <- sub("^HG19FASTA ", "", linea_software)
    }
  }
  
  #########################################33
  
  # STEP 5 MarkDuplicates
  if(!(file.exists(sprintf("%s%s_dedup_reads.bam", path_to_MitoR, patient_ID )))) {
    system2("java", sprintf("-jar %s MarkDuplicates -I %s --VALIDATION_STRINGENCY SILENT --CREATE_INDEX True --ASSUME_SORTED True -M %s%s_marked_dup_metrics.txt -O %s%s_dedup_reads.bam", PICARD, sorted_BAM, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)
  }
  
  #-------------------------------------------------------------------
  #Intentando arreglar el problema de que no se genera nada en el step 6:
  #samtools view -H /media/4tb1/Daniela/CNVs/MuestrasCNVs/46432/OMICsdo/46432_dedup_reads.bam
  
  args <- c("AddOrReplaceReadGroups",
            #sprintf("I=%s%s_dedup_reads.bam", path_to_MitoR, patient_ID),
            sprintf("I=%s", sorted_BAM),
            sprintf("O=%s%s_dedup_reads_with_RG.bam", path_to_MitoR, patient_ID),
            "RGID=1",
            "RGLB=lib1",
            "RGPL=illumina",
            "RGPU=unit1",
            sprintf("RGSM=%s", patient_ID))
  
  if(!(file.exists(sprintf("%s%s_dedup_reads_with_RG.bam", path_to_MitoR, patient_ID )))) {
    system2("java", args = c("-jar", PICARD, args))
    system(sprintf("%s index %s%s_dedup_reads_with_RG.bam", Samtools, path_to_MitoR, patient_ID))
    
  }
  
  # STEP 6 HaplotypeCaller
  #system2("java", sprintf("-jar %s HaplotypeCaller -I %s%s_dedup_reads.bam -O %s%s_raw_variants.vcf -ip 100 -R %s", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference), stdout = TRUE, wait = TRUE)
  
  comando_gatk <- paste(sprintf("java -jar %s HaplotypeCaller", GATK),
                        "-R", reference,
                        "-I", sprintf("%s%s_dedup_reads_with_RG.bam", path_to_MitoR, patient_ID),
                        "-O", sprintf("%s%s_raw_variants.vcf", path_to_MitoR, patient_ID),
                        sep = " ")
  if(!(file.exists(sprintf("%s%s_raw_variants.vcf", path_to_MitoR, patient_ID )))) {
    system(comando_gatk)
  }
  
  
  # STEP 7 SelectVariants - SNPs
  if(!(file.exists(sprintf("%s%s_raw_snps.vcf", path_to_MitoR, patient_ID )))) {
    system2("java", sprintf("-jar %s SelectVariants -V %s%s_raw_variants.vcf -O %s%s_raw_snps.vcf -R %s -select-type SNP", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference), stdout = TRUE, wait = TRUE)
  }
  
  # STEP 8 SelectVariants - Indels
  if(!(file.exists(sprintf("%s%s_raw_indels.vcf", path_to_MitoR, patient_ID )))) {
    system2("java", sprintf("-jar %s SelectVariants -V %s%s_raw_variants.vcf -O %s%s_raw_indels.vcf -R %s -select-type INDEL", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference), stdout = TRUE, wait = TRUE)
  }
  
  # STEP 9 VariantFiltration - SNPs
  system2("java", sprintf("-jar %s VariantFiltration -V %s%s_raw_snps.vcf -O %s%s_filtered_snps.vcf  -R %s --filter-expression 'QD %s || FS %s || MQ %s || MQRankSum %s || ReadPosRankSum %s' --filter-name 'mitor_indel_filter'", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS), stdout = TRUE, wait = TRUE)
  
  # STEP 10
  system2("java", sprintf("-jar %s VariantFiltration -V %s%s_raw_indels.vcf -O %s%s_filtered_indels.vcf  -R %s --filter-expression 'QD %s || FS %s || ReadPosRankSum %s' --filter-name 'mitor_snp_filter'", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS), stdout = TRUE, wait = TRUE)
  
  # STEP 11
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s%s_filtered_indels.vcf -O %s%s_filtered_indels_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s%s_filtered_snps.vcf -O %s%s_filtered_snps_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)
  
  # STEP 12
  if(!(file.exists(sprintf("%s%s_MitoR.vcf", path_to_MitoR, patient_ID )))) {
    system2("java", sprintf("-jar %s MergeVcfs -I %s%s_filtered_snps_tomerge.vcf -I %s%s_filtered_indels_tomerge.vcf -O %s%s_OMICsdo.vcf", PICARD, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)
  }
  
  # Deletes the unnecessary files
  #file.remove(sprintf("%s%s_raw_variants.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_dedup_reads_with_RG.bam", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_dedup_reads_with_RG.bam.bai", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_dedup_reads.bam", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_dedup_reads.bai", path_to_MitoR, patient_ID))
  
  file.remove(sprintf("%s%s_raw_snps.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_indels.vcf", path_to_MitoR, patient_ID))
  #file.remove(sprintf("%s%s_marked_dup_metrics.txt", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_indels.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_snps.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_indels_tomerge.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_snps_tomerge.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_variants.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_snps.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_indels.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_indels.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_snps.vcf.idx", path_to_MitoR, patient_ID))
  #file.remove(sprintf("%s%s_MitoR.vcf.idx", path_to_MitoR, patient_ID))
  
  
  return(sprintf("%s%s_OMICsdo.vcf", path_to_MitoR, patient_ID))
}
