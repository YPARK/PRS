Traits := $(shell cat traits.txt)
ODIR := docs/share

all:
	@echo "Try ... make step1"
	@echo "Then ... make step2"
	@echo "Then ... make step3"

CUTOFF := 2 3 4

CIS_WINDOW := 1e3 1e4 1e5

################################################################
# Produce PRS matrix on the Ruzicka data
step1: jobs/step1/prs.short.gz

step1_long: jobs/step1/prs.long.gz

JOBS_STEP1 := $(foreach geno, Ruzicka UK10K, $(foreach trt, $(Traits), jobs/step1/prs_$(geno)/$(trt).job))

jobs/step1/prs.short.gz: $(JOBS_STEP1)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N STEP1_PRS -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

# % = $(GENO)/$(TRT) --> $(ODIR)/PRS_%.prs.gz
jobs/step1/prs_%.job:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo Rscript make.prs.R data/$(shell echo $* | awk -F'/' '{ print $$1 }')_Plink LD.info.txt data/GWAS/$(shell echo $* | awk -F'/' '{ print $$2 }').bed.gz $(ODIR)/PRS_$*.prs.gz > $@

################################################################
# PRS matrix with partitioning by DEG with different p-value cutoffs
step2: jobs/step2/partition.short.gz

step2_long: jobs/step2/partition.long.gz

JOBS_STEP2 := $(foreach geno, Ruzicka UK10K, $(foreach trt, $(Traits), $(foreach cutoff, $(CUTOFF), $(foreach cis, $(CIS_WINDOW), jobs/step2/partition_$(geno)/$(trt).$(cutoff).$(cis).job))))

jobs/step2/partition.short.gz: $(JOBS_STEP2)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=16:00:00 -b y -j y -N STEP2_PARTITION -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

# % = $(GENO)/$(TRT).$(CUTOFF).$(CIS) -> $(ODIR)/PARTITION_%.prs.gz:
jobs/step2/partition_%.job:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ -d $(dir $(ODIR)/PARTITION_$*) ] || mkdir -p $(dir $(ODIR)/PARTITION_$*)
	@echo Rscript make.prs.partition.R data/$(shell echo $* | awk -F'/' '{ print $$1 }')_Plink LD.info.txt data/GWAS/$(shell echo $* | awk -F'/' '{ print $$2 }' | awk -F'.' '{ print $$1 }').bed.gz data/gene_annot_limma_pval.tab.gz $(shell echo $* | awk -F'/' '{ print $$2 }' | awk -F'.' '{ print $$2 " " $$3 }') $(ODIR)/PARTITION_$*.prs.gz > $@

################################################################
# Calculate anova
step3: jobs/step3/aov.short.gz

JOBS_STEP3 := $(foreach geno, Ruzicka UK10K, $(foreach cutoff, $(CUTOFF), $(foreach cis, $(CIS_WINDOW), jobs/step3/aov_$(geno)/$(cutoff)_$(cis).job)))

jobs/step3/aov.short.gz: $(JOBS_STEP3)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N STEP3_AOV -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

# % = $(geno)/$(cutoff)_$(dist)
jobs/step3/aov_%.job:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ -d $(dir $(ODIR)/AOV_$*) ] || mkdir -p $(dir $(ODIR)/AOV_$*)
	@echo Rscript make.prs.partition.aov.R $(ODIR)/PARTITION_$(shell echo $* | awk -F'/' '{ print $$1 }')/ $(shell echo $* | awk -F'/' '{ print $$2 }' | sed 's/_/ /g') $(ODIR)/AOV_$(shell echo $* | awk -F'/' '{ print $$1 FS "aov_" $$2 ".txt.gz" }') > $@

################################################################
jobs/%.long.gz: jobs/%.short.gz
	gzip -cd $< | awk 'system(" ! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=36:00:00 -b y -j y -N long_$(shell basename $* | sed 's/\//_/g') -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@
