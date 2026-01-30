fix_filter_values <- function(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS,
                              QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS) {
  params <- c(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  i=1
  for (i in seq_along(params)) {
    params[i] <- stringr::str_replace(params[i], ",", ".")
    if (!grepl("\\.", params[i])) {
      params[i] <- paste0(params[i], ".0")
    }
  }
  return(params)
}
