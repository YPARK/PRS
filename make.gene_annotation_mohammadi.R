#!/usr/bin/env Rscript

source("Util.R")
library(tidyverse)
library(data.table)

data.file = "data/final_scores.RData"

if(!file.exists(data.file)) {
    log.msg("couldn't find: %s\n", data.file)
    q()
}

.dat = load.data(data.file)
pval.tab = .dat[["DE.logPvals.limmma"]] # three m's?
logFC.tab = .dat[["DE.logFC.limma"]]

info.tab = read.gene.info(rownames(pval.tab))
pval.tab = annotate.gene.tab(pval.tab, info.tab)
logFC.tab = annotate.gene.tab(logFC.tab, info.tab)

write_tsv(logFC.tab, "data/gene_annot_limma_logFC.tab.gz")
write_tsv(pval.tab, "data/gene_annot_limma_pval.tab.gz")
