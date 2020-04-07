Traits := $(shell cat traits.txt)
ODIR := docs/share

all:
	@echo "Try ... make step1"
	@echo "Then ... make step2"
	@echo "Then ... make step3"

# Produce PRS matrix on the Ruzicka data
step1: $(foreach trt, $(Traits), $(ODIR)/PRS_Ruzicka/$(trt).prs.gz)

# % = $(GENO)/$(TRT)
$(ODIR)/PRS_%.prs.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript make.prs.R data/$(shell echo $* | awk -F'/' '{ print $$1 }')_Plink LD.info.txt data/GWAS/$(shell echo $* | awk -F'/' '{ print $$2 }').bed.gz $@

# PRS matrix with partitioning by DEG with different p-value cutoffs
step2: $(foreach trt, $(Traits), $(foreach cutoff, 2 3 4, $(foreach cis, 1e3 1e4 1e5, $(ODIR)/PARTITION_Ruzicka/$(trt).$(cutoff).$(cis).prs.gz)))

# % = $(GENO)/$(TRT)_$(CUTOFF)_$(CIS)
$(ODIR)/PARTITION_%.prs.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript make.prs.partition.R data/$(shell echo $* | awk -F'/' '{ print $$1 }')_Plink LD.info.txt data/GWAS/$(shell echo $* | awk -F'/' '{ print $$2 }' | awk -F'.' '{ print $$1 }').bed.gz data/gene_annot_limma_pval.tab.gz $(shell echo $* | awk -F'/' '{ print $$2 }' | awk -F'.' '{ print $$2 " " $$3 }')  $@
