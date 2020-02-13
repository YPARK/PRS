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

