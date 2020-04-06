#!/usr/bin/env Rscript

argv = commandArgs(trailingOnly = TRUE)

if(length(argv) < 8) {
    q()
}

## PRS.FILE = 'PRS_gtex_hg38/ctg_ad.prs.gz'
## EXPR.FILE = 'data/gtex_expr_hg38/Brain_Cortex.v8.normalized_expression.bed.gz'
## ANNOT.FILE = 'data/annotation.txt.gz'
## P.VAL.CUTOFF = 1e-2
## N.CONTROL = 1000
## K = 10
## SUBCAT = 'BP'
## OUT.FILE = 'temp.txt.gz'

PRS.FILE = argv[1]
EXPR.FILE = argv[2]
ANNOT.FILE = argv[3]
P.VAL.CUTOFF = as.numeric(argv[4])
N.CONTROL = as.integer(argv[5])
K = as.integer(argv[6])
SUBCAT = argv[7]
OUT.FILE = argv[8]

library(dplyr)
library(tidyr)
source('Util.R')

dir.create(dirname(OUT.FILE), recursive=TRUE, showWarnings = FALSE)

prs.tab = readr::read_tsv(PRS.FILE)
expr.tab = readr::read_tsv(EXPR.FILE)

tis.name = basename(EXPR.FILE) %>%
    gsub(pattern='.v8.normalized_expression.bed.gz',
         replacement='')

## .cn = c('ensembl_gene_id',
##         'entrezgene',
##         'hgnc_symbol',
##         'transcript_start',
##         'transcript_end',
##         'percentage_gc',
##         'gs_name',
##         'gs_id',
##         'gs_cat',
##         'gs_subcat')

.ct = 'cccdddcccc'
annot.tab = readr::read_tsv(ANNOT.FILE, col_types = .ct) %>%
    mutate(transcript_start = as.integer(transcript_start)) %>% 
    mutate(transcript_end = as.integer(transcript_end))    
gc()

take.gtex.iid <- function(s) {
    ret = strsplit(s, split = '-')[[1]]
    ret[2]
}

gene.info = expr.tab %>%
    select(`#chr`, end, gene_id) %>%
    tidyr::separate('gene_id', c('ensembl_gene_id', 'remove'), sep='[.]') %>%
    select(-remove) %>%
    rename(tss = end) %>%
    mutate(idx = 1:n())

X = expr.tab %c% (-(1:4)) %>%
    as.matrix() %>%
    t()

x.info = tibble(iid = rownames(X)) %>%
    mutate(x.pos = 1:n()) %>%
    mutate(iid = sapply(iid, take.gtex.iid)) %>%
    mutate(iid = as.character(iid))

y.info = prs.tab %>%
    select(iid) %>%
    mutate(iid = sapply(iid, take.gtex.iid)) %>%
    mutate(iid = as.character(iid)) %>%
    mutate(y.pos = 1:n())

matched = x.info %>%
    left_join(y.info, by = 'iid') %>%
    na.omit()

X = X %r% matched$x.pos %>%
    scale()

X[is.na(X)] = 0

Y = (prs.tab %r% matched$y.pos) %>%
    mutate(y = mu/sig) %>%
    select(y) %>%
    scale() %>%
    as.matrix()

## Pathway annotation
A.tab = annot.tab %>% filter(gs_subcat == SUBCAT)

C.annot = gene.info %>%
    left_join(A.tab) %>%
    filter(!is.na(hgnc_symbol))

gs.tab = C.annot %>%
    mutate(len = transcript_end - transcript_start) %>% 
    group_by(gs_subcat, gs_name) %>%
    summarize(n = length(unique(ensembl_gene_id)),
              gc = round(mean(percentage_gc), 2),
              len = round(mean(len))) %>%
    ungroup()

C.mat = C.annot %>%
    mutate(len = transcript_end - transcript_start) %>% 
    select(idx, len, percentage_gc, gs_name) %>%
    mutate(val = 1) %>% 
    tidyr::spread(key = gs_name, value = val, fill = 0) %>%
    ungroup() %>%
    mutate(len = exp(scale(log(len)))) %>% 
    mutate(percentage_gc = exp(scale(log(percentage_gc))))

gc()

## 0. Initial round of correlation to establish control genes
raw.tab = zqtl::calc.qtl.stat(X, Y) %>%
    select(-y.col, -resid.se) %>%
    rename(p.val.raw = p.val,
           beta.raw = beta,
           se.raw = se)

## 1. Remove confounders
control.genes = raw.tab %>%
     filter(p.val.raw >= P.VAL.CUTOFF) %>%
     top_n(N.CONTROL, 1 - p.val.raw)

C.ruv = X %c% control.genes$x.col
c.svd = zqtl::take.ld.svd(C.ruv)
lm.out = lm(X ~ (c.svd$U %c% 1:K) + 1)
X.ruv = residuals(lm.out) %>%
    scale()

x.idx = C.mat$idx %>% unlist(use.names = FALSE)
xx = X.ruv %c% x.idx

C = C.mat[, -(1:3)] %>% as.matrix()
C.names = colnames(C)
colnames(C) = NULL

## gene by gene z-score
## zz = zqtl::fast.z.cov(X.ruv, Y) %>%
##     (function(z) { colnames(z) = 'z'; z })

## 2. Test statistics for the meanshift
test.mean.shift <- function(.xx, .y,
                            .C, .C.names,
                            perm.size = 200,
                            n.perm = 1000) {

    fast.z.cov = zqtl::fast.z.cov
    N.j = apply(.C, 2, sum)

    zz.obs = fast.z.cov(.xx, .y)
    stat.obs = sweep(t(.C) %*% zz.obs, 1, sqrt(N.j), `/`)

    nfalse = rep(1, length(stat.obs))
    ntot = 1

    .yy = rep(.y, perm.size) %>%
        matrix(ncol=perm.size, byrow=FALSE)

    stat.perm.mean = stat.obs * 0
    stat.perm.sd = abs(stat.obs) * 0

    for(b in 1:n.perm) {

        .yy = .yy %>% apply(MARGIN=2, FUN=sample)

        .zz.perm = fast.z.cov(.xx, .yy)
        .stat.perm = sweep(t(.C) %*% .zz.perm, 1, sqrt(N.j), `/`)

        .nfalse = sweep(abs(.stat.perm), 1, abs(stat.obs), `>`) %>%
            apply(MARGIN=1, FUN=sum)

        nfalse = nfalse + .nfalse
        ntot = ntot + perm.size

        stat.perm.mean = stat.perm.mean + apply(.stat.perm, 1, sum)
        stat.perm.sd = stat.perm.sd + apply(.stat.perm^2, 1, sum)

        min.pval = min(nfalse / ntot)

        log.msg('Perm = %d, min p-value = %.2e', ntot, min.pval)
    }

    stat.perm.mean = stat.perm.mean / ntot
    stat.perm.sd = sqrt(stat.perm.sd / ntot - stat.perm.mean^2)

    ret = tibble(gs_name = .C.names,
                 z = as.numeric(stat.obs),
                 null.mean = as.numeric(stat.perm.mean),
                 null.sd = as.numeric(stat.perm.sd),
                 pval = as.numeric(nfalse) / ntot)

}

out.tab = test.mean.shift(xx, Y, C, C.names, 500, 1000)

out.tab = out.tab %>%
    left_join(gs.tab) %>%
    mutate(tis = tis.name)

readr::write_tsv(out.tab, OUT.FILE)
log.msg('Done')
