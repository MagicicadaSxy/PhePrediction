

high_lm_cv <- function(m){
  print(paste("trait",m,"loading"))

  phe <- SNP_tep_phe(m,SNPclass = "high_SNP",setwd=setwd) %>%
    SNP_dos_base(.)

  lm <- model_SNP_lm(phe)
  print(paste("MLR_R2:",lm))
  cv <- model_SNP_cv(phe)
  print(paste("ENR_R2:",cv))

  R2 <- c(lm,cv)
  return(R2)
}