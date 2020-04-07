#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly = TRUE)

if(length(argv) != 4) { q() }

source('Util.R')
library(tidyverse)
library(data.table)

DATA.DIR = argv[1] # e.g., "docs/share/PARTITION_Ruzicka/"
CUTOFF = argv[2]   # e.g., "3"
CIS.DIST = argv[3] # e.g., "1e4"
OUT.FILE = argv[4] # e.g., "out.txt.gz"

read.aov <- function(ff) {

    .trait = basename(ff) %>%
        str_split(pattern="[.]") %>%
        (function(x) x[[1]][1])

    .dt = fread("gzip -cd " %&% ff) %>%
        spread(key="celltype", value="score")

    .dt = .dt[, iid := NULL]
    .dt[is.na(.dt)] = 0
    .lm = lm(.all ~ ., data=.dt)

    .trim <- function(x) as.character(x) %>%
                             str_remove_all("`") %>% 
                             str_trim("both")

    .beta = .lm %>% summary() %>%
        coefficients()

    .r1 = rownames(.beta) %>% .trim()

    .aov = .lm %>% aov() %>% summary() %>%
        (function(x) x[[1]])

    .pve = .aov %>%
        select(`Sum Sq`) %>%
        unlist(use.names=FALSE)

    vtot = sum(.pve)
    .pve = .pve / vtot

    .r2 = rownames(.aov) %>% .trim()

    .beta = .beta %>% as_tibble() %>% 
        mutate(celltype = .r1)

    .aov = .aov %>% as_tibble() %>%
        mutate(celltype = .r2) %>%
        mutate(pve = .pve)

    ret = left_join(.beta, .aov, by = "celltype") %>%
        mutate(trait = .trait) %>%
        na.omit() %>% 
        rename(beta = Estimate, se = `Std. Error`,
               pval.t = `Pr(>|t|)`, pval.pve = `Pr(>F)`) %>%
        select(trait, celltype, beta, se, pve, pval.t, pval.pve)

}

.files = list.files(DATA.DIR,
                    paste(CUTOFF, CIS.DIST, "prs.gz", sep="."),
                    full.names=TRUE)

aov.tab = lapply(.files, read.aov) %>%
    bind_rows()

write_tsv(aov.tab, OUT.FILE)
