#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly = TRUE)

if(length(argv) < 1) {
    q()
}

out.file = argv[1]
library(dplyr)

## Take gene names
ensembl = biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL',
                           host='useast.ensembl.org',
                           path='/biomart/martservice', dataset='hsapiens_gene_ensembl')

ensembl.hs = biomaRt::useDataset('hsapiens_gene_ensembl', mart=ensembl)

.attr = c('ensembl_gene_id',
          'entrezgene',
          'hgnc_symbol',
          'transcript_start',
          'transcript_end',
          'percentage_gene_gc_content')

coding.genes = biomaRt::getBM(attributes=.attr,
                              filters='biotype',
                              values=c('protein_coding'),
                              mart=ensembl.hs)

coding.genes = coding.genes %>%
    group_by(ensembl_gene_id, entrezgene) %>% 
    summarize(hgnc_symbol = paste(unique(hgnc_symbol), collapse='|'),
              transcript_start = min(transcript_start),
              transcript_end = max(transcript_end),
              percentage_gc = mean(percentage_gene_gc_content)) %>%
    ungroup() %>%
    filter(nchar(hgnc_symbol) > 0) %>%
    na.omit()

gc()

## Pathway and GO annotation
gs.sub = c('TFT', 'MF', 'BP', 'CC',
           'CP:KEGG','CP:REACTOME', 'CP:BIOCARTA')

msigdb = msigdbr::msigdbr(species = 'Homo sapiens') %>%
    rename(entrezgene = entrez_gene) %>%
    rename(hgnc_symbol = human_gene_symbol) %>%
    filter(gs_subcat %in% gs.sub) %>%
    select(-sources) %>% 
    select(-gene_symbol) %>% 
    select(-species_name)

## Match ensembl -> pathway
out.tab = coding.genes %>%
    left_join(msigdb) %>%
    filter(!is.na(gs_id))

readr::write_tsv(out.tab, path = out.file)

