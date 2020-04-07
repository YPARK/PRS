Traits := $(shell cat traits.txt)
ODIR := docs/share

all:
	@echo "Try ... make step1"
	@echo "Then ... make step2"
	@echo "Then ... make step3"

################################################################
# Produce PRS matrix on the Ruzicka data
step1: $(foreach trt, $(Traits), $(ODIR)/PRS_Ruzicka/$(trt).prs.gz)

# % = $(GENO)/$(TRT)
$(ODIR)/PRS_%.prs.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript make.prs.R data/$(shell echo $* | awk -F'/' '{ print $$1 }')_Plink LD.info.txt data/GWAS/$(shell echo $* | awk -F'/' '{ print $$2 }').bed.gz $@

################################################################
# PRS matrix with partitioning by DEG with different p-value cutoffs
step2: jobs/step2/partition.short.gz

JOBS_STEP2 := $(foreach geno, Ruzicka UK10K, $(foreach trt, $(Traits), $(foreach cutoff, 2 3 4, $(foreach cis, 1e3 1e4 1e5, jobs/step2/partition_$(geno)/$(trt).$(cutoff).$(cis).job))))

jobs/step2/partition.short.gz: $(JOBS_STEP2)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=8:00:00 -b y -j y -N STEP2_PARTITION -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

# % = $(GENO)/$(TRT).$(CUTOFF).$(CIS) -> $(ODIR)/PARTITION_%.prs.gz:
jobs/step2/partition_%.job:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ -d $(dir $(ODIR)/PARTITION_$*) ] || mkdir -p $(dir $(ODIR)/PARTITION_$*)
	@echo Rscript make.prs.partition.R data/$(shell echo $* | awk -F'/' '{ print $$1 }')_Plink LD.info.txt data/GWAS/$(shell echo $* | awk -F'/' '{ print $$2 }' | awk -F'.' '{ print $$1 }').bed.gz data/gene_annot_limma_pval.tab.gz $(shell echo $* | awk -F'/' '{ print $$2 }' | awk -F'.' '{ print $$2 " " $$3 }') $(ODIR)/PARTITION_$*.prs.gz > $@
