exp<-openxlsx2::read_xlsx("formatted/merged_eqtl.xlsx")
out<-extract_outcome_data(
  snps = exp$SNP,
  outcomes = "ieu-a-1120",
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.3,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  splitsize = 10000,
  proxy_splitsize = 500
)
out <- TwoSampleMR::format_data(
  dat=Ovarian_cancer_full_outcome,
  type = "outcome",
  snps = exp$SNP,
  header = TRUE,
  phenotype_col = "outcome",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "pval.outcome",
  chr_col = "chr.outcome",
  samplesize_col = "samplesize.outcome",
  eaf_col = "eaf.outcome",
  pos_col = "pos.outcome",
  ncase_col = "ncase.outcome",
  ncontrol_col = "ncontrol.outcome"
)
  
dat1<-harmonise_data(exp, out, action = 2)
res_mi1=generate_odds_ratios(mr(dat,method_list = c("mr_wald_ratio","mr_ivw")))
res_mi1$B<-p.adjust(res_mi1$pval,"bonferroni")
d<-subset(res_mi1,B<0.05)
writexl::write_xlsx(res_mi, "res_mi_oc,xlsx")




