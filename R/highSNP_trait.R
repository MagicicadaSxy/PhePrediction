highSNP_trait <- function(m){
  print(paste("trait",m,"loading"))
  n <- m
  phe <- SNP_tep_phe(m,SNPclass = "high_SNP",setwd=setwd) %>%
    SNP_dos_base(.) %>%
    set_SNP_trait(m=n,.)

  lm <- model_SNP_Trait_lm(phe)
  print(paste("MLR_R2:",lm))
  cv <- model_SNP_Trait_cv(phe)
  print(paste("ENR_R2:",cv))

  R2 <- c(lm,cv)
  return(R2)
}