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

RNORM <- function(d1, d2) {
    matrix(rnorm(d1 * d2), nrow = d1, ncol = d2)
}

################################################################
## clean potential genetic signals
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

## read.geno.data <- function(ld.idx,
##                            cis.dist,
##                            geno.dir,
##                            temp.dir,
##                            ld.file = 'LD.info.txt',
##                            maf.cutoff = .01) {

##     ld.info = read.table(ld.file, header = TRUE)

##     if(ld.idx < 0) return(NULL)
##     if(ld.idx > nrow(ld.info)) return(NULL)

##     ld.info = ld.info %r% ld.idx
##     ld.lb = ld.info[1, 'start']
##     ld.ub = ld.info[1, 'stop']
##     chr.int = ld.info[1, 'chr'] %>%
##         gsub(pattern = 'chr', replacement = '')

##     plink = subset.plink(geno.dir %&&% '/chr' %&&% chr.int,
##                          chr.int, ld.lb, ld.ub, temp.dir)

##     maf = apply(plink$BED, 2, mean, rm.na = TRUE) %>%
##         (function(x) data.frame(maf = x)) %>%
##         mutate(pos = 1:n())

##     bim.cols = c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')

##     BIM = plink$BIM %>%
##         (function(x) { colnames(x) = bim.cols; return(x); }) %>%
##         mutate(chr=gsub(chr, pattern='chr', replacement='')) %>%
##         mutate(chr='chr' %&&% chr) %>% 
##         mutate(pos = 1:n()) %>%
##         left_join(maf, by = 'pos') %>%
##         filter(snp.loc >= ld.lb, snp.loc <= ld.ub) %>%
##         na.omit() %>%
##         filter(maf >= maf.cutoff, maf <= (1 - maf.cutoff))

##     X = (plink$BED %c% BIM$pos) %>% scale() %>%
##         (function(x) { x[is.na(x)] = 0; return(x) })

##     ret = list(X = X, BIM = BIM)
##     return(ret)
## }

## ################################################################
## ## standardize z-scores: First, fitting
## ## z ~ N(R*(mu*I), tau^2*R)
## ## then standardize
## ## z = (z - mu*I) / tau
## scale.zscore <- function(Z, V.t, D, stdize = TRUE) {

##     .p = ncol(V.t)
##     .K = nrow(V.t)

##     if(.K < 10) {
##         return(Z)
##     }

##     DV.t = sweep(V.t, 1, D, `*`)
##     .y = sweep(V.t, 1, D, `/`) %*% Z
##     .x = DV.t %*% matrix(1, ncol(V.t), 1)

##     .xx = .x / sqrt(.p) ## to scale down for numerical stability
##     .num = t(.xx) %*% .y
##     .denom = sum(.xx * .x)

##     .mu = .num / .denom

##     z.mean = t(DV.t) %*% (.x %*% .mu)
##     if(stdize) {
##         z.sd = (.y - .x %*% .mu) %>%
##             (function(.c) apply(.c^2, 2, sum) / (.K - 1)) %>%
##             sqrt()
##         ret = sweep(Z - z.mean, 2, pmax(z.sd, 1e-8), `/`)
##     } else {
##         ret = Z - z.mean
##     }
##     return(ret)
## }

## ################################################################
## ## project z score matrix onto different LD block
## rotate.Z <- function(.svd.1, .svd.0, z) {

##     U.1 = .svd.1$U
##     D.1 = .svd.1$D
##     Vt.1 = .svd.1$V.t

##     U.0 = .svd.0$U
##     D.0 = .svd.0$D
##     Vt.0 = .svd.0$V.t

##     ret = (Vt.0 %*% z)
##     ret = sweep(U.0, 2, D.0, `/`) %*% ret
##     ret = t(U.1) %*% ret
##     ret = sweep(t(Vt.1), 2, D.1, `*`) %*% ret

##     ret = scale.zscore(ret, Vt.1, D.1)
## }

