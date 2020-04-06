#!/usr/bin/env Rscript

argv = commandArgs(trailingOnly = TRUE)

if(length(argv) < 7) {
    q()
}

source('Util.R')
library(tidyverse)
library(data.table)

PLINK.DIR = argv[1]            # e.g., "data/Ruzicka_Plink/"
LD.FILE = argv[2]              # e.g., "LD.info.txt"
GWAS.FILE = argv[3]            # e.g., "data/GWAS/clozuk_scz.bed.gz"
ANNOT.FILE = argv[4]           # e.g., "data/gene_annot_limma_pval.tab.gz"
CUTOFF = as.numeric(argv[5])   # e.g., 3
CIS.DIST = as.numeric(argv[6]) # e.g., 1e5
OUT.FILE = argv[7]             # e.g.,  "out.txt.gz"

################################################################

annot.pval.tab = read_tsv(ANNOT.FILE) %>%
    gather(key="celltype", value="pval",
           -`#chr`, -tss, -tes, -hgnc_symbol) %>%
    filter(pval > CUTOFF) %>%
    mutate(lb = pmax(tss - CIS.DIST, 0)) %>%
    mutate(ub = tes + CIS.DIST)

effect.tab = fread(GWAS.FILE)

REF.PLINK.DIR = "data/1KG_EUR/"

################################################################
TEMP.DIR = ('temp_' %&&% OUT.FILE %&&% '_' %&&% GWAS.FILE) %>%
    gsub(pattern = '/', replacement = '_')

dir.create(TEMP.DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(dirname(OUT.FILE), recursive = TRUE, showWarnings = FALSE)

ld.tab = readr::read_tsv(LD.FILE, col_types='cii') %>%
    mutate(query = chr %&&% ':' %&&% (start + 1) %&&% '-' %&&% (stop - 1)) %>%
    as.data.frame()

effect.tab = fread(GWAS.FILE) # this is fast enough

.names = colnames(effect.tab)

if("lodds" %in% .names) {
    effect.tab[, beta := lodds]
}

if(!("chr" %in% .names) && ("#CHR" %in% .names)) {
    effect.tab[, chr := `#CHR`]
}

if(!("chr" %in% .names) && ("#chr" %in% .names)) {
    effect.tab[, chr := `#chr`]
}

if(!("snp.loc" %in% .names)){
    effect.tab[, snp.loc := `stop`]
}

effect.tab = effect.tab[]


#' @param .svd SVD result
#' @param tau regularization parameter
take.DinvVt <- function(.svd, tau) {
    D = .svd$D
    U = .svd$U
    V.t = .svd$V.t
    .tau = tau / sqrt(nrow(U))
    sweep(V.t, 1, D + .tau, `/`)
}

#' @param .svd SVD result
#' @param zz z-score vector/matrix
#' @param tau regularization parameter
pred.coeff <- function(.svd, zz, tau) {
    W.t = take.DinvVt(.svd, tau)
    .beta = t(W.t) %*% (W.t %*% zz)
}

#' @param .svd SVD result
#' @param zz z-score vector/matrix
#' @param tau regularization parameter
pred.prs <- function(.svd, zz, tau) {
    W.t = take.DinvVt(.svd, tau)
    .svd$U %*% (W.t %*% zz)
}

#' Local polygenic score prediction with regularization parameter
#' tau tuned by cross-comparison with additional genotype matrix
#' @param zz z-score vector
#' @param X primary genotype matrix
#' @param X.test additional genotype matrix
#' @return risk score
pred.prs.cv <- function(zz, X, X.test){

    nn = nrow(X)
    pp = ncol(X)
    nn.test = nrow(X.test)

    stopifnot(ncol(X) == ncol(X.test))
    stopifnot(nrow(zz) == pp)

    zz = zz + 1e-8 * RNORM(pp, ncol(zz))
    xx = X + 1e-8 * RNORM(nn, pp) # add small noise
    xx.test = X.test + 1e-8 * RNORM(nn.test, pp) # add small noise

    svd.train = zqtl::take.ld.svd(xx, eigen.reg = 0)
    svd.test = zqtl::take.ld.svd(xx.test, eigen.reg = 0)
    xx.test = sweep(svd.test$U, 2, svd.test$D, `*`) %*% svd.test$V.t
    y.test = pred.prs(svd.test, zz, tau=0)

    cor.score <- function(tau) {
        beta = pred.coeff(svd.train, zz, tau)
        y.pred = xx.test %*% beta
        .cor = cor(y.test, y.pred, method="spearman") %>% as.numeric()
        sum(1 - .cor)
    }

    opt.out = optim(0, cor.score, method="Brent", lower=0, upper=nn)
    tau.opt = opt.out$par

    log.msg("Found optimal tau %.2f with score %.2f",
            tau.opt, 1 - opt.out$value)

    y.opt = pred.prs(svd.train, zz, tau=tau.opt) %>%
        as.data.table()

    colnames(y.opt) = "score"

    list(y = y.opt, tau = tau.opt, cor = 1 - opt.out$value)
}

#' @param .snps tab with beta, se
take.z.vector <- function(.snps) {
    .snps %>%
        mutate(z = beta/se) %>%
        select(z) %>%
        as.matrix()
}


#' Build PRS within r-th local LD block
#' @param r index
#' @return risk score with iid
take.prs.ld.ct <- function(r) {

    .chr = ld.tab[r, "chr"]
    .chr.num = str_remove(.chr, "chr") %>% as.integer()
    .lb = ld.tab[r, "start"]
    .ub = ld.tab[r, "stop"]


    .effect.1 = effect.tab[chr == .chr & snp.loc >= .lb & snp.loc <= .ub]
    .effect.2 = effect.tab[chr == .chr.num & snp.loc >= .lb & snp.loc <= .ub]

    .effect = rbind(.effect.1, .effect.2)
    .effect = .effect[, chr := NULL][]

    log.msg("Take an LD block: %d SNPs", nrow(.effect))


    .plink = subset.plink(PLINK.DIR %&% "/" %&% .chr, .chr, .lb, .ub, TEMP.DIR)

    .plink.ref = subset.plink(REF.PLINK.DIR %&% "/" %&% .chr, .chr, .lb, .ub, TEMP.DIR)

    log.msg("Read genotypes: n=%d, n=%d", nrow(.plink$BED), nrow(.plink.ref$BED))

    .matched = match.plink(.plink, .plink.ref)
    .plink = .matched$lhs
    .plink.ref = .matched$rhs

    log.msg("After matching --> n=%d, n=%d", nrow(.plink$BED), nrow(.plink.ref$BED))

    valid.snps = .plink$BIM %>%
        select(-missing, -rs) %>%
        mutate(x.pos = 1:n()) %>%
        (function(x) left_join(.effect, x, by = 'snp.loc')) %>%
        filter(!is.na(x.pos)) %>%
        filter((plink.a1 == a1 & plink.a2 == a2) | (plink.a2 == a1 & plink.a1 == a2)) %>%
        mutate(gwas.flip = if_else(plink.a1 == a1, 1.0, -1.0)) %>%
        mutate(beta = gwas.flip * beta) %>%
        select(-gwas.flip)

    if(nrow(valid.snps) == 0) {
        log.msg("No SNPs that satisfy our criteria")
        return(list(cor = NA))
    }

    valid.snps = valid.snps %>%
        group_by(rs) %>%
        slice(which.max(abs(beta))) %>%
        ungroup()

    .iid = .plink$FAM

    ## Take subset of genotype matrix
    X = .plink$BED %c% valid.snps$x.pos
    X.ref = .plink.ref$BED %c% valid.snps$x.pos

    stopifnot(nrow(X) == nrow(.plink$FAM))

    zz = take.z.vector(valid.snps)
    prs.tot = pred.prs.cv(zz, X, X.ref)

    prs.tot$y = prs.tot$y[, iid := .iid$iid]
    prs.tot$y = prs.tot$y[, celltype := ".all"]
    prs.tot = prs.tot$y

    ## construct PRS vectors annotated by cell types
    .pval.ld.tab = annot.pval.tab %>%
        filter(`#chr` == .chr.num) %>%
        filter(tss >= .lb, tes <= .ub)

    take.ct.snps <- function(.ct) {

        .annot.ct = .pval.ld.tab %>%
            filter(celltype == .ct)

        ret = tibble()
        for(j in 1:nrow(.annot.ct)){
            .lb.ct = .annot.ct$lb[j]
            .ub.ct = .annot.ct$ub[j]

            .ret.j = valid.snps %>%
                filter(snp.loc > .lb.ct, snp.loc < .ub.ct)

            ret = bind_rows(ret, .ret.j)
        }

        ret %>%
            group_by(snp.loc) %>%
            slice(which.max(abs(beta))) %>%
            ungroup()
    }

    prs.ct = NULL

    for(.ct in unique(.pval.ld.tab$celltype)) {
        .snps.ct = take.ct.snps(.ct)

        if(nrow(.snps.ct) < 10) next # at least 10 SNPs

        .zz.ct = take.z.vector(.snps.ct)

        .xx.ct = .plink$BED[ ,.snps.ct$x.pos, drop = FALSE]
        .xx.test.ct = .plink.ref$BED[ ,.snps.ct$x.pos, drop = FALSE]

        .prs.ct = pred.prs.cv(.zz.ct, .xx.ct, .xx.test.ct)
        .prs.ct = .prs.ct$y
        .prs.ct = .prs.ct[, iid := .iid$iid]
        .prs.ct = .prs.ct[, celltype := .ct]

        prs.ct = rbind(prs.ct, .prs.ct)
    }

    return(rbind(prs.tot, prs.ct))
}

prs.list = lapply(1:nrow(ld.tab), take.prs.ld.ct)

out.tab = bind_rows(prs.list) %>% as.data.table
out.tab = out.tab[, .(score = sum(score)), by = .(iid, celltype)]

write_tsv(out.tab, path = OUT.FILE)
