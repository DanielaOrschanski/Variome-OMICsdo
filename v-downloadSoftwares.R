downloadFastQC <- function(soft_directory) {
  #soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  
  tryCatch(
    expr = {
      softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
      linea_software <- grep("(?i)FastQC", softwares, ignore.case = TRUE, value = TRUE)
      
      if(length(nchar(linea_software))==0) {
        stop()
      } else {
        FastQC <<- strsplit(linea_software, " ")[[1]][[2]]
        system2(FastQC, "--help")
        message("FastQC was already downloaded")
        return(FastQC)
      }
      
    },
    error = function(e) {
      message("FastQC download and installation is about to begin...
              Please be patient, it may take a while.")
      print(e)
      
      # En caso de haberlo descargado y hay algun problema con el ejecutable,
      # eliminamos y descargamos de nuevo
      if (file.exists(sprintf('%s/FastQC', soft_directory))) {
        system2("rm", sprintf('-r %s/FastQC', soft_directory))
      }
      
      tryCatch(
        expr = {
          # Proceso de Descarga de FastQC
          dir.create(sprintf("%s/FastQC", soft_directory))
          setwd(sprintf("%s/FastQC", soft_directory))
          URL <- ("https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip")
          
          system2("wget", URL, wait = TRUE, stdout = NULL, stderr = NULL)
          
          file <- basename(URL)
          file_dir <- file.path(sprintf("%s/FastQC", soft_directory), file)
          filedc <- substr(file, start = 0, stop = (nchar(file) - 4))
          
          # Descomprimo y elimino el ZIP
          system2("unzip", sprintf("%s -d %s/FastQC", file_dir, soft_directory), wait = TRUE)
          system2("rm", file_dir)
          
          # Definicion de la variable FastQC con el ejecutable
          FastQC <<- sprintf('%s/FastQC/FastQC/fastqc', soft_directory)
          system2(FastQC, "--help")
          
          # En caso de que la descarga haya sido exitosa, agregamos el path al archivo TXT
          softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
          
          if(length(nchar(softwares[-grep("FastQC", softwares, ignore.case = TRUE)])) != 0 ) {
            softwares <- softwares[-grep("FastQC", softwares, ignore.case = TRUE)]
          }
          softwares_actualizado <- c(softwares, sprintf("FastQC %s", FastQC))
          # Reescribimos el archivo
          write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
          
          
          message("FastQC download and installation completed successfully")
          
          return(FastQC)
        },
        error = function(e) {
          message("An error occured while performing the FastQC download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install wget unzip
------------------------------------------------------")
          print(e)
        },
        
        finally = {
          message("-.Message from FastQC")
        }
      )
    }
    
  )
}


#################################################################################################

downloadTrimGalore <- function(soft_directory) {
  
  #soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  tryCatch(
    expr = {
      softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
      linea_software <- grep("(?i)TrimGalore", softwares, ignore.case = TRUE, value = TRUE)
      
      if(length(nchar(linea_software))==0) {
        stop()
      } else {
        TrimGalore <<- strsplit(linea_software, " ")[[1]][[2]]
        system2(TrimGalore, "--help")
        message("TrimGalore was already downloaded")
        return(TrimGalore)
      }
      
    },
    error = function(e) {
      # En caso de haberlo descargado y hay algun problema con el ejecutable, eliminamos y descargamos de nuevo
      if (file.exists(sprintf('%s/TrimGalore', soft_directory))) {
        system2("rm", sprintf('-r %s/TrimGalore', soft_directory))
        print("There is a problem with the TrimGalore exe file. It will be removed and download again")
      }
      
      message("TrimGalore download and installation is about to begin...
              Please be patient, it may take a while.")
      print(e)
      
      tryCatch(
        expr = {
          # Primero hay que verificar que este descargado FastQC y Cutadapt.
          downloadFastQC()
          
          # Proceso de Descarga de TrimGalore
          dir.create(sprintf("%s/TrimGalore", soft_directory))
          URL <- "https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz"
          setwd(sprintf("%s/TrimGalore", soft_directory))
          system2("wget" ,sprintf("https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -O trim_galore.tar.gz"), wait = TRUE, stdout = NULL, stderr = NULL)
          system2("tar", "xvzf trim_galore.tar.gz")
          
          TrimGalore <<- sprintf("%s/TrimGalore/TrimGalore-0.6.10/trim_galore", soft_directory)
          system2(TrimGalore, "--help")
          
          # En caso de que la descarga haya sido exitosa, agregamos el path al archivo TXT
          softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
          
          if(length(nchar(softwares[-grep("TrimGalore", softwares, ignore.case = TRUE)])) != 0 ) {
            softwares <- softwares[-grep("TrimGalore", softwares, ignore.case = TRUE)]
          }
          softwares_actualizado <- c(softwares, sprintf("TrimGalore %s", TrimGalore))
          write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
          
          return(TrimGalore)
        },
        
        error = function(e) {
          message("An error occured while performing the TrimGalore download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install cutadapt wget tar
------------------------------------------------------")
          print(e)
        },
        
        finally = {
          message("-.Message from TrimGalore")
        }
      )
    },
    finally = {
      message("TrimGalore download and installation completed successfully")
    }
  )
}



########################################################################################
downloadSamtools <- function(soft_directory) {
  
  tryCatch(
    {
      system2(sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory))
    },
    error = function(e) {
      message("The installation of Samtools will begin now. Check if all required packages are downloaded.")
      print(e)
      
      dir.create(sprintf("%s/Samtools", soft_directory))
      dir <- sprintf("%s/Samtools", soft_directory)
      URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
      system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
      
      system2("bzip2", sprintf("-d %s/Samtools/samtools-1.16.1.tar.bz2", soft_directory))
      system2("tar", c("-xvf", sprintf("%s/Samtools/%s -C %s", soft_directory, list.files(sprintf("%s/Samtools", soft_directory)), dir)))
      
      file.remove(sprintf("%s/samtools-1.16.1.tar", dir))
      
      samtools_dir <- sprintf("%s/samtools-1.16.1", dir)
      system(paste("cd", shQuote(samtools_dir), "&& ./configure"))
      system(paste("cd", shQuote(samtools_dir), "&& make"))
      
      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", soft_directory))
      
      Samtools <- sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory)
      softwares_actualizado <- c(softwares, sprintf("Samtools %s", Samtools))
      write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
      
    },
    warning = function(w) {
      message("The installation of Samtools will begin now. Check if all required packages are downloaded.")
      
      dir.create(sprintf("%s/Samtools", soft_directory))
      dir <- sprintf("%s/Samtools", soft_directory)
      URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
      system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
      
      system2("bzip2", sprintf("-d %s/Samtools/samtools-1.16.1.tar.bz2", soft_directory))
      system2("tar", c("-xvf", sprintf("%s/Samtools/%s -C %s", soft_directory, list.files(sprintf("%s/Samtools", soft_directory)), dir)))
      
      file.remove(sprintf("%s/samtools-1.16.1.tar", dir))
      
      samtools_dir <- sprintf("%s/samtools-1.16.1", dir)
      system(paste("cd", shQuote(samtools_dir), "&& ./configure"))
      system(paste("cd", shQuote(samtools_dir), "&& make"))
      
      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
      
      Samtools <- sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory)
      softwares_actualizado <- c(softwares, sprintf("Samtools %s", Samtools))
      write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
      
    },
    
    finally = {
      message("-.Message from Samtools")
      Samtools <<- sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory)
    }
  )
  
  return(sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory))
}

downloadBWA <- function(soft_directory) {
  
  omicsdo_sof <- soft_directory
  print(omicsdo_sof)
  
  
  tryCatch(
    expr = {
      system(sprintf('%s/BWA/usr/bin/bwa', soft_directory))
    },
    error = function(e) {
      message("Installation of BWA will now begin. Check if all the required packages are downloaded.")
      print(e)
      
      dir.create(sprintf("%s/BWA", omicsdo_sof))
      
      bwa_url1 <- "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm"
      bwa_dir1 <- file.path(omicsdo_sof, "BWA")
      system2("wget", args = c(bwa_url1, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)
      
      dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", omicsdo_sof))
      bwa_url2 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm"
      bwa_dir2 <- file.path(omicsdo_sof, "BWA/bwa-0.7.17-lp154.6.1.src")
      system2("wget", args = c(bwa_url2, "-P", bwa_dir2), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.src.rpm | cpio -D %s -idmv", bwa_dir2, bwa_dir2), wait = TRUE)
      
      bwa_url3 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm"
      system2("wget", args = c(bwa_url3, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.i586.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)
      
      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      BWA <- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("BWA %s", BWA))
      #write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
      
    },
    warning = function(w) {
      message("Installation of BWA will now begin. Check if all the required packages are downloaded.")
      
      dir.create(sprintf("%s/BWA", omicsdo_sof))
      
      bwa_url1 <- "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm"
      bwa_dir1 <- file.path(omicsdo_sof, "BWA")
      system2("wget", args = c(bwa_url1, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)
      
      dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", omicsdo_sof))
      bwa_url2 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm"
      bwa_dir2 <- file.path(omicsdo_sof, "BWA/bwa-0.7.17-lp154.6.1.src")
      system2("wget", args = c(bwa_url2, "-P", bwa_dir2), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.src.rpm | cpio -D %s -idmv", bwa_dir2, bwa_dir2), wait = TRUE)
      
      bwa_url3 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm"
      system2("wget", args = c(bwa_url3, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.i586.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)
      
      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
      BWA <- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("BWA %s", BWA))
      #write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
      
    },
    
    finally = {
      BWA <<- sprintf('%s/BWA/usr/bin/bwa', soft_directory)
      message("-.Message from BWA")
    }
  )
  
  return(sprintf('%s/BWA/usr/bin/bwa', soft_directory))
}

# GATK ##############################################################################

downloadGATK <- function(soft_directory) {
  
  omicsdo_sof <- soft_directory
  
  if(file.exists(sprintf('%s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof))) {
    system(sprintf('java -jar %s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof))
  } else {
    message("Installation of GATK will now begin. Check if all the required packages are downloaded.")
    dir.create(sprintf("%s/GATK", omicsdo_sof))
    URL <- 'https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip'
    system2("wget", args = c(URL, "-P", paste(omicsdo_sof, "/GATK", sep="")), wait = TRUE, stdout = NULL, stderr = NULL)
    file <- basename(URL)
    file_dir <- paste(omicsdo_sof, "/GATK/", file , sep = "")
    filedc <- substr(file, start = 0, stop = (nchar(file) - 4))
    unzip(zipfile = file_dir, exdir = paste(omicsdo_sof, "/GATK", sep=""))
    
    #Escribo el path en el txt
    softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
    GATK <- sprintf('%s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof)
    softwares_actualizado <- c(softwares, sprintf("GATK %s", GATK))
    #write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
    write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
    
  }
  
  GATK <<- sprintf('%s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof)
  message("-.Message from GATK")
  
  return(GATK)
}



# PICARD ##################################################################################


downloadPICARD <- function(soft_directory) {
  
  #omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  omicsdo_sof <- soft_directory
  
  if (file.exists(sprintf('%s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof))) {
    sprintf('java -jar %s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof)
    
  } else {
    message("Installation of PICARD will now begin. Check if all the required packages are downloaded.")
    dir.create(sprintf("%s/PICARD", omicsdo_sof))
    URL <- "https://github.com/broadinstitute/picard/archive/refs/tags/2.27.5.tar.gz"
    system2("wget", args = c(URL, "-P", paste(omicsdo_sof, "/PICARD", sep="")), wait = TRUE, stdout = NULL, stderr = NULL)
    system2("gzip" , sprintf("-d %s/PICARD/2.27.5.tar.gz", omicsdo_sof))
    system2("tar" , sprintf("-xvf %s/PICARD/2.27.5.tar -C %s/PICARD", omicsdo_sof, omicsdo_sof))
    file.remove(sprintf("%s/PICARD/2.27.5.tar", omicsdo_sof))
    
    URL2 <- "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar"
    dir2 <- sprintf("%s/PICARD/picard-2.27.5", omicsdo_sof)
    system2("wget", args = c(URL2, "-P", dir2), wait = TRUE, stdout = NULL, stderr = NULL)
    
    #Escribo el path en el txt
    softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
    PICARD <- sprintf('%s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof)
    softwares_actualizado <- c(softwares, sprintf("PICARD %s", PICARD))
    #write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
    write(softwares_actualizado, file = sprintf("%s/path_to_soft.txt", soft_directory))
    
  }
  
  PICARD <<- sprintf('%s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof)
  message("-.Message from PICARD")
  
  return(PICARD)
}
