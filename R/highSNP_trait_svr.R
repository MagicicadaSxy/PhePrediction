highSNP_trait_svr <- function(m){
  print(paste("trait",m,"loading"))
  n <- m
  phe <- SNP_tep_phe(m,SNPclass = "high_SNP",setwd=setwd) %>%
    SNP_dos_base(.) %>%
    set_SNP_trait(m=n,.)

  svr <- model_SNP_Trait_svr(phe)
  print("========================================================================")
  return(svr)
}