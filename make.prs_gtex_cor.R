#!/usr/bin/env Rscript

argv = commandArgs(trailingOnly = TRUE)

if(length(argv) < 7) {
    q()
}

## PRS.FILE = 'PRS_gtex_hg38/ctg_ad.prs.gz'
## EXPR.FILE = 'data/gtex_expr_hg38/Brain_Cortex.v8.normalized_expression.bed.gz'
## COV.FILE = 'data/gtex_cov_hg38/Brain_Cortex.v8.covariates.txt'
## P.VAL.CUTOFF = 1e-2
## N.CONTROL = 1000
## K = 5
## OUT.FILE = 'temp.txt.gz'

PRS.FILE = argv[1]
EXPR.FILE = argv[2]
COV.FILE = argv[3]
P.VAL.CUTOFF = as.numeric(argv[4])
N.CONTROL = as.integer(argv[5])
K = as.integer(argv[6])
OUT.FILE = argv[7]

library(dplyr)
library(tidyr)
source('Util.R')

prs.tab = readr::read_tsv(PRS.FILE)
expr.tab = readr::read_tsv(EXPR.FILE)
cov.tab = readr::read_tsv(COV.FILE)

take.gtex.iid <- function(s) {
    ret = strsplit(s, split = '-')[[1]]
    ret[2]
}

gene.info = expr.tab %>%
    select(`#chr`, end, gene_id) %>%
    rename(tss = end) %>%
    mutate(idx = 1:n())

X = expr.tab %c% (-(1:4)) %>%
    as.matrix() %>%
    t()

C = cov.tab %c% (-1) %>%
    as.matrix() %>%
    t()

x.info = tibble(iid = rownames(X)) %>%
    mutate(x.pos = 1:n()) %>%
    mutate(iid = sapply(iid, take.gtex.iid)) %>%
    mutate(iid = as.character(iid))

c.info = tibble(iid = rownames(C)) %>%
    mutate(c.pos = 1:n()) %>%
    mutate(iid = sapply(iid, take.gtex.iid)) %>%
    mutate(iid = as.character(iid))

y.info = prs.tab %>%
    select(iid) %>%
    mutate(iid = sapply(iid, take.gtex.iid)) %>%
    mutate(iid = as.character(iid)) %>% 
    mutate(y.pos = 1:n())

matched = x.info %>%
    left_join(c.info, by = 'iid') %>% 
    left_join(y.info, by = 'iid') %>%
    na.omit()

X = X %r% matched$x.pos %>%
    scale()

X[is.na(X)] = 0

C = C %r% matched$c.pos %>%
    scale()

C[is.na(C)] = 0

Y = (prs.tab %r% matched$y.pos) %>%
    mutate(y = mu/sig) %>%
    select(y) %>%
    scale() %>% 
    as.matrix()

## 0. Initial round of correlation to establish control genes
raw.tab = zqtl::calc.qtl.stat(X, Y) %>%
    select(-y.col, -resid.se) %>% 
    rename(p.val.raw = p.val,
           beta.raw = beta,
           se.raw = se)

## 1. Remove confounders
lm.out = lm(X ~ C + 1)
X.res = residuals(lm.out) %>%
    scale()

## alternative RUV approach
control.genes = raw.tab %>%
     filter(p.val.raw >= P.VAL.CUTOFF) %>% 
     top_n(N.CONTROL, 1 - p.val.raw)

C.ruv = X %c% control.genes$x.col
c.svd = zqtl::take.ld.svd(C.ruv)
lm.out = lm(X ~ (c.svd$U %c% 1:K) + 1)
X.ruv = residuals(lm.out) %>%
    scale()

## 2. Next round of correlation
out.1 = zqtl::fast.z.cov(X.res, Y) %>%
    (function(z) { colnames(z) = 'z.res'; z }) %>% 
    as.data.frame() %>% 
    mutate(p.val.res = sapply(z.res, zqtl::zscore.pvalue)) %>% 
    mutate(x.col = 1:n())

out.2 = zqtl::fast.z.cov(X.ruv, Y) %>%
    (function(z) { colnames(z) = 'z.ruv'; z }) %>% 
    as.data.frame() %>% 
    mutate(p.val.ruv = sapply(z.ruv, zqtl::zscore.pvalue)) %>% 
    mutate(x.col = 1:n())

## output
out.tab = out.2 %>%
    left_join(out.1) %>% 
    left_join(raw.tab) %>%
    (function(x) gene.info %>% rename(x.col=idx) %>% left_join(x)) %>%
    select(-x.col)

readr::write_tsv(out.tab, path = OUT.FILE)
log.msg('Done')
