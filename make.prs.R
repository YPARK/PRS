#!/usr/bin/env Rscript

argv = commandArgs(trailingOnly = TRUE)

if(length(argv) != 4) {
    q()
}

source('Util.R')
library(tidyverse)
library(data.table)

PLINK.DIR = argv[1]              # e.g., PLINK.DIR = "data/Ruzicka_Plink/"
LD.FILE   = argv[2]              # e.g., LD.FILE = "LD.info.txt"
GWAS.FILE = argv[3]              # e.g., GWAS.FILE = "data/GWAS/ssgac_speeding.bed.gz"
OUT.FILE  = argv[4]              # e.g., OUT.FILE = "out.txt.gz"

REF.PLINK.DIR = "data/1KG_EUR/"

if(file.exists(OUT.FILE)) {
    log.msg("File exists: %s", OUT.FILE)
    q()
}

dir.create(dirname(OUT.FILE), recursive=TRUE, showWarning=FALSE)

################################################################
TEMP.DIR = ('temp_' %&&% OUT.FILE %&&% '_' %&&% GWAS.FILE) %>%
    gsub(pattern = '/', replacement = '_')

dir.create(TEMP.DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(dirname(OUT.FILE), recursive = TRUE, showWarnings = FALSE)

ld.tab = readr::read_tsv(LD.FILE, col_types='cii') %>%
    mutate(query = chr %&&% ':' %&&% (start + 1) %&&% '-' %&&% (stop - 1)) %>%
    as.data.frame()

#' @param .chr
#' @param .lb
#' @param .ub
read.gwas.tab <- function(.chr, .lb, .ub) {

    ret = fread("tabix -h " %&% GWAS.FILE %&% " " %&%
                .chr %&% ":" %&% .lb %&% "-" %&% .ub)

    .names = colnames(ret)

    if("lodds" %in% .names) {
        ret[, beta := lodds]
    }

    if(!("chr" %in% .names) && ("#CHR" %in% .names)) {
        ret[, chr := `#CHR`]
    }

    if(!("chr" %in% .names) && ("#chr" %in% .names)) {
        ret[, chr := `#chr`]
    }

    if(!("snp.loc" %in% .names)){
        ret[, snp.loc := `stop`]
    }

    ret = ret[]
    return(ret)
}

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

#' Build PRS within r-th local LD block
#' @param r index
#' @return risk score with iid
take.prs.ld <- function(r) {

    .chr = ld.tab[r, "chr"]
    .chr.num = str_remove(.chr, "chr") %>% as.integer()
    .lb = ld.tab[r, "start"]
    .ub = ld.tab[r, "stop"]

    ## .effect.1 = effect.tab[chr == .chr & snp.loc >= .lb & snp.loc <= .ub]
    ## .effect.2 = effect.tab[chr == .chr.num & snp.loc >= .lb & snp.loc <= .ub]

    .effect.1 = read.gwas.tab(.chr, .lb, .ub)
    .effect.2 = read.gwas.tab(.chr.num, .lb, .ub)

    .effect = rbind(.effect.1, .effect.2)
    .effect = .effect[, chr := NULL][]

    log.msg("Take an LD block: %d SNPs", nrow(.effect))

    if(nrow(.effect) == 0) return(data.table())

    .plink = subset.plink(PLINK.DIR %&% "/" %&% .chr, .chr, .lb, .ub, TEMP.DIR)
    .plink.ref = subset.plink(REF.PLINK.DIR %&% "/" %&% .chr, .chr, .lb, .ub, TEMP.DIR)

    log.msg("Read genotypes: n=%d, n=%d", nrow(.plink$BED), nrow(.plink.ref$BED))

    .matched = match.plink(.plink, .plink.ref)

    .plink = .matched$lhs
    .plink.ref = .matched$rhs

    log.msg("After matching --> n=%d, n=%d", nrow(.plink$BED), nrow(.plink.ref$BED))

    if(is.null(.plink)) return(data.table())
    if(is.null(.plink.ref)) return(data.table())

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
        return(data.table())
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

    zz = valid.snps %>%
        mutate(z = beta/se) %>%
        select(z) %>%
        as.matrix()

    out = pred.prs.cv(zz, X, X.ref)
    out$y[, iid := .iid$iid]
    out$iid = .iid
    return(out$y)
}

prs.list = lapply(1:nrow(ld.tab), take.prs.ld)
.dt = bind_rows(prs.list) %>% as.data.table
out.tab = .dt[, .(score = sum(score)), by = iid]
write_tsv(out.tab, path = OUT.FILE)
