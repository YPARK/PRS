options(stringsAsFactors = FALSE)

`%c%` <- function(mat, cols) mat[, cols, drop = FALSE]
`%r%` <- function(mat, rows) mat[rows, , drop = FALSE]
`%&%` <- function(a,b) paste(a, b, sep = '')
`%&&%` <- function(a,b) paste(a, b, sep = '')

.unlist <- function(...) unlist(..., use.names = FALSE)
.zeros <- function(n1, n2) matrix(0, n1, n2)

sigmoid <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x) - log(1- x)

.eval <- function(str) eval(parse(text = str))

log.msg <- function(...) {
    ss = as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
    flush(stderr())
}

load.data <- function(fileName){
    load(fileName)
    mget(ls()[ls() != "fileName"])
}

RNORM <- function(d1, d2) {
    matrix(rnorm(d1 * d2), nrow = d1, ncol = d2)
}


#' @param plink.hdr
#' @param chr
#' @param plink.lb
#' @param plink.ub
#' @param temp.dir
subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

    require(dplyr)

    .error <- function(e) {
        print(e)
        log.msg('No QTL here!\n')
        return(NULL)
    }

    chr = unlist(chr)
    plink.lb = unlist(plink.lb) %>% as.integer()
    plink.ub = unlist(plink.ub) %>% as.integer()

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num = gsub(pattern = 'chr', replacement = '', chr) %>%
            as.integer()

        out.hdr = temp.dir %&% '/plink'

        plink.cmd = sprintf('rm -f %s.*; ./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', out.hdr, plink.hdr, chr.num, plink.lb, plink.ub, out.hdr)

        ## print(plink.cmd)

        system(plink.cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

        plink = zqtl::read.plink(out.hdr)
        colnames(plink$BIM) = c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
        colnames(plink$FAM) = c('fam', 'iid', 'father', 'mother', 'sex.code', '.pheno')

        if(any(is.logical(plink$BIM$plink.a1))) {
            plink$BIM$plink.a1 = 'T'
        }

        if(any(is.logical(plink$BIM$plink.a2))) {
            plink$BIM$plink.a2 = 'T'
        }
        return(plink)
    }

    plink = tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                      error = .error)
    return(plink)
}

# Match the right-hand-side PLINK to the left-hand-side
#' @param .plink.lhs
#' @param .plink.rhs
match.plink <- function(.plink.lhs, .plink.rhs) {

    if(is.null(.plink.lhs) && !is.null(.plink.rhs)) { return(.plink.rhs) }
    if(!is.null(.plink.rhs) && is.null(.plink.rhs)) { return(.plink.lhs) }
    if(is.null(.plink.lhs) && is.null(.plink.rhs)) { return(NULL) }

    lhs = .plink.lhs$BIM %>%
        as_tibble() %>%
        mutate(lhs.pos = 1:n()) %>%
        rename(lhs.plink.a1 = plink.a1, lhs.plink.a2 = plink.a2) %>%
        select(-missing, -rs)

    rhs = .plink.rhs$BIM %>%
        as_tibble() %>%
        mutate(rhs.pos = 1:n()) %>%
        rename(rhs.plink.a1 = plink.a1, rhs.plink.a2 = plink.a2) %>%
        select(-missing, -rs)

    .matched = left_join(lhs, rhs, by = c("chr", "snp.loc")) %>% 
        na.omit()

    if(nrow(.matched) < 1) return(NULL)

    .matched = .matched %>%
        dplyr::filter(((lhs.plink.a1 == rhs.plink.a1) & (lhs.plink.a2 == rhs.plink.a2)) |
                      ((lhs.plink.a2 == rhs.plink.a1) & (lhs.plink.a1 == rhs.plink.a2))) %>%
        arrange(chr, snp.loc)

    if(nrow(.matched) < 1) return(NULL)

    ret.lhs = .plink.lhs
    ret.rhs = .plink.rhs

    ret.lhs$BIM = ret.lhs$BIM[.matched$lhs.pos, , drop = FALSE]
    ret.lhs$BED = ret.lhs$BED[ , .matched$lhs.pos, drop = FALSE]

    ret.rhs$BIM = ret.rhs$BIM[.matched$rhs.pos, , drop = FALSE]
    ret.rhs$BED = ret.rhs$BED[ , .matched$rhs.pos, drop = FALSE]

    flip.tab = ret.lhs$BIM %>%
        mutate(lhs.pos = 1:n()) %>%
        left_join(ret.rhs$BIM %>% mutate(rhs.pos = 1:n()),
                  by = c("chr", "snp.loc"),
                  suffix = c(".lhs", ".rhs")) %>%
        filter(plink.a1.lhs != plink.a1.rhs)

    ret.rhs$BIM[flip.tab$rhs.pos, ] <- ret.lhs$BIM[flip.tab$lhs.pos, ]

    flip.bed = ret.rhs$BED[, flip.tab$rhs.pos]
    zero.idx = flip.bed <= 0.5
    two.idx = flip.bed >= 1.5
    flip.bed[two.idx] = 0
    flip.bed[zero.idx] = 2
    ret.rhs$BED[, flip.tab$rhs.pos] = flip.bed

    list(lhs = ret.lhs, rhs = ret.rhs)
}

################################################################
#' calculate univariate effect sizes and p-values
#'
#' @name calc.qtl.stat
#' @usage calc.qtl.stat(xx, yy)
#' @param xx n x p genotype matrix
#' @param yy n x t phenotype matrix
#' @param se.min mininum standard error (default: 1e-8)
#' @param verbose (default: FALSE)
#'
#' @return summary statistics matrix
#'
#' @export
#'
calc.qtl.stat <- function(xx, yy, se.min = 1e-8, verbose = FALSE) {

    loadNamespace("dplyr")
    loadNamespace("tidyr")

    zscore.pvalue <- function(x) {
        2 * pnorm(abs(x), lower.tail=FALSE)
    }

    .xx = scale(xx)
    .yy = scale(yy)

    rm.na.zero <- function(xx) {
        return(replace(xx, is.na(xx), 0))
    }

    ## cross-product is much faster than covariance function
    n.obs = crossprod(!is.na(.xx), !is.na(.yy))
    beta.mat = crossprod(.xx %>% rm.na.zero(), .yy %>% rm.na.zero()) / n.obs

    if(verbose){
        log.msg('Computed cross-products')
    }

    ## residual standard deviation
    resid.se.mat = matrix(NA, ncol(.xx), ncol(.yy))

    for(k in 1:ncol(.yy)) {

        beta.k = beta.mat[, k]
        yy.k = .yy[, k]
        err.k = sweep(sweep(.xx, 2, beta.k, `*`), 1, yy.k, `-`)
        se.k = apply(err.k, 2, sd, na.rm = TRUE)

        if(verbose) {
            log.msg('Residual on the column %d', k)
        }
        resid.se.mat[, k] = se.k + se.min
    }

    ## organize as consolidated table
    y.cols = 1:ncol(yy)
    colnames(beta.mat) = y.cols
    colnames(n.obs) = y.cols
    colnames(resid.se.mat) = y.cols

    beta.tab = beta.mat %>%
        as.data.frame() %>%
        dplyr::mutate(x.col = 1:n()) %>%
        tidyr::gather(key = 'y.col', value = 'beta', y.cols)

    resid.se.tab = resid.se.mat %>%
        as.data.frame() %>%
        dplyr::mutate(x.col = 1:n()) %>%
        tidyr::gather(key = 'y.col', value = 'resid.se', y.cols)

    nobs.tab = n.obs %>%
        as.data.frame() %>%
        dplyr::mutate(x.col = 1:n()) %>%
        tidyr::gather(key = 'y.col', value = 'n', y.cols)

    out.tab = beta.tab %>%
        dplyr::left_join(nobs.tab, by = c('x.col', 'y.col')) %>%
        dplyr::left_join(resid.se.tab, by = c('x.col', 'y.col')) %>%
        dplyr::mutate(se = resid.se/sqrt(n)) %>%
        dplyr::mutate(p.val = zscore.pvalue(beta/se))

    out.tab = out.tab %>%
        dplyr::mutate(x.col = as.integer(x.col)) %>%
        dplyr::mutate(y.col = as.integer(y.col))

    return(out.tab)
}

################################################################
#' @param genes a vector of gene names
read.gene.info <- function(genes) {

    ensembl = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                               host="grch37.ensembl.org",
                               path="/biomart/martservice",
                               dataset="hsapiens_gene_ensembl")

    ensembl.hs = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)

    .attr = c("hgnc_symbol",
              "chromosome_name",
              "transcription_start_site",
              "transcript_start",
              "transcript_end")

    .dt = biomaRt::getBM(attributes=.attr,
                         filters="hgnc_symbol",
                         values=genes,
                         mart=ensembl.hs) %>%
        as.data.table

    info.tab = .dt[, .(tss=min(transcript_start),
                      tes=max(transcript_end)),
                  by=.(hgnc_symbol, chromosome_name)]

    return(info.tab)
}

#' @param .tab
#' @param info.tab
annotate.gene.tab <- function(.tab, info.tab) {

    .name.dt = tibble(hgnc_symbol = rownames(.tab)) %>%
        mutate(r = 1:n()) %>%
        left_join(info.tab, by = "hgnc_symbol") %>%
        na.omit()

    rr = .name.dt$r

    .name.dt = .name.dt %>%
        dplyr::select(chromosome_name, tss, tes, hgnc_symbol) %>%
        as.data.table

    ret = cbind(.name.dt, .tab[rr, ]) %>%
        dplyr::filter(chromosome_name %in% 1:22) %>%
        arrange(chromosome_name, tss) %>%
        rename(`#chr` = chromosome_name)
}
