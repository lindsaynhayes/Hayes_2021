exp_genes <- read.csv('data/unique_consensus_expressed.csv')
DEG <- read.csv('data/DEGgenes.csv')
notDEG <- setdiff(exp_genes$gene, DEG$gene)
cbind(nrow(DEG), length(notDEG), nrow(exp_genes))
message("3441 total mouse DEG mapped to ", length(unique(DEG$gene)), " unique genes")

GWAS_cat <- read.csv('data/GWAS_catalog_gene_centered.csv')
gwas.all <- GWAS_cat[grep("autism", GWAS_cat$disease),]
gwas <- gwas.all[which(gwas.all$gene %in% exp_genes$gene), ]
message(nrow(gwas.all), " total gwas risk genes mapped to ", nrow(gwas), " expressed mouse genes")

DEG.gwas <- length(intersect(DEG$gene, gwas$gene))
DEG.notGwas <- length(setdiff(DEG$gene, gwas$gene))
notDEG.gwas <- length(intersect(notDEG, gwas$gene))
notDEG.notGwas <- length(setdiff(notDEG, gwas$gene))

pct.DEG <- length(intersect(DEG$gene, gwas$gene))/length(DEG$gene)*100
pct.notDEG <- length(intersect(notDEG, gwas$gene))/length(notDEG)*100
mat <- rbind(c(DEG.gwas, DEG.notGwas), c(notDEG.gwas, notDEG.notGwas))
rownames(mat) <- c("DEG & GWAS    |  DEG & notGWAS", "notDEG & GWAS | notDEG & notGWAS")
mat
fit <- chisq.test(mat)
fit
print(paste0(round(pct.DEG, 2), '% of DEG are associated with Autism'))
print(paste0(round(pct.notDEG, 2), '% of the expressed but not DEG are associated with Autism')) 
if(pct.DEG > pct.notDEG & fit$p.value < 0.05){
  print(paste0('Risk variations ARE significantly enriched in DEG P-Value = ', round(fit$p.value,3)))
} else {
  print(paste0('Risk variations are NOT significantly enriched in DEG P-Value = ', round(fit$p.value,3)))
}