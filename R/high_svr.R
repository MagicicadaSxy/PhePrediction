high_svr <- function(m){
  print(paste("trait",m,"loading"))
  n <- m

  phe <- SNP_tep_phe(m,SNPclass = "high_SNP",setwd=setwd)
  phe <- SNP_dos_base(phe)

  svr <- model_SNP_svr(phe)
  print("========================================================================")
  return(svr)
}