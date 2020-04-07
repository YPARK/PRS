VCF := $(foreach chr, $(shell seq 1 22), $(foreach ext, vcf.gz vcf.gz.tbi, $(foreach data, 50 109 merged, data/$(data)/chr$(chr).$(ext))))

all: human_g1k_v37.fasta human_g1k_v37.fasta.fai $(VCF)

human_g1k_v37.fasta:
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -O $@.gz
	gzip -d $@.gz

human_g1k_v37.fasta.fai: human_g1k_v37.fasta
	samtools faidx $<

%-sort.vcf.gz: %.vcf.gz
	vcf-sort $< | bgzip -c > $@

%-sort.vcf.gz.tbi: %-sort.vcf.gz
	tabix $<

%.vcf.gz.tbi: %.vcf.gz
	tabix $<

data_50 := Ruzicka_50_samples_OmniExpress_10-2018

data/50/chr%.vcf.gz: $(data_50)-sort.vcf.gz $(data_50)-sort.vcf.gz.tbi 
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	bcftools view --regions $* $< -Ou | bcftools norm --fasta-ref human_g1k_v37.fasta -Ou | bcftools annotate -Oz -o $@ -I '%CHROM:%POS:%ID'

data_109 := Ruzicka_109_samples_OmniExpress_9-2018

data/109/chr%.vcf.gz: $(data_109)-sort.vcf.gz $(data_109)-sort.vcf.gz.tbi 
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	bcftools view --regions $* $< -Ou | bcftools norm --fasta-ref human_g1k_v37.fasta -Ou | bcftools annotate -Oz -o $@ -I '%CHROM:%POS:%ID'

data/merged/chr%.vcf.gz: data/109/chr%.vcf.gz data/50/chr%.vcf.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	bcftools merge --merge all $^ -O b | bcftools sort -O z > $@

## Download data from the Michigan Imputation Server
DOWNLOAD := $(foreach chr, $(shell seq 1 22), imputed/chr$(chr).zip)
PASS := "7rl*P$fOAB5vEv"

download: imputed imputed/qcreport.html $(DOWNLOAD)

imputed/qcreport.html: imputed
	wget https://imputationserver.sph.umich.edu/share/results/d32118e34938472becae64690594ec1c/qcreport.html -O $@

imputed:
	mkdir -p $@

imputed/chr1.zip:
	wget https://imputationserver.sph.umich.edu/share/results/3ccb1b8291051968346fc7e84dfdf83/chr_1.zip -O $@

imputed/chr10.zip:
	wget https://imputationserver.sph.umich.edu/share/results/71b9f31c585f3432f86bc4a33b22a1e9/chr_10.zip -O $@

imputed/chr11.zip:
	wget https://imputationserver.sph.umich.edu/share/results/d818d491eec808d999e16cfc5dc2b4c3/chr_11.zip -O $@

imputed/chr12.zip:
	wget https://imputationserver.sph.umich.edu/share/results/d9023629e4c5c66d9fe0291102bab18e/chr_12.zip -O $@

imputed/chr13.zip:
	wget https://imputationserver.sph.umich.edu/share/results/156357ab8758ab57ac0ee38b84f33a06/chr_13.zip -O $@

imputed/chr14.zip:
	wget https://imputationserver.sph.umich.edu/share/results/c2874f670ac502a4509e799d1b0d7870/chr_14.zip -O $@

imputed/chr15.zip:
	wget https://imputationserver.sph.umich.edu/share/results/ad6dc4d2568fac3f5cd1924a6b924a4b/chr_15.zip -O $@

imputed/chr16.zip:
	wget https://imputationserver.sph.umich.edu/share/results/5c3e8392ebd552355540c179bc49be44/chr_16.zip -O $@

imputed/chr17.zip:
	wget https://imputationserver.sph.umich.edu/share/results/d369f59944003d38376deb29acbbd44d/chr_17.zip -O $@

imputed/chr18.zip:
	wget https://imputationserver.sph.umich.edu/share/results/2d3a172538a23bd1a33734d82d8b43db/chr_18.zip -O $@

imputed/chr19.zip:
	wget https://imputationserver.sph.umich.edu/share/results/fe8b931c0356ee07b7c63454ce2901ce/chr_19.zip -O $@

imputed/chr2.zip:
	wget https://imputationserver.sph.umich.edu/share/results/83dc42a82bfe1a553e0cb8943a96334b/chr_2.zip -O $@

imputed/chr20.zip:
	wget https://imputationserver.sph.umich.edu/share/results/f304182849f215926d0d4b73356d0ded/chr_20.zip -O $@

imputed/chr21.zip:
	wget https://imputationserver.sph.umich.edu/share/results/78e58f51566d2b45c4a7f5d2be4d8e12/chr_21.zip -O $@

imputed/chr22.zip:
	wget https://imputationserver.sph.umich.edu/share/results/78ae7a0a93940019b0b18dfef7210557/chr_22.zip -O $@

imputed/chr3.zip:
	wget https://imputationserver.sph.umich.edu/share/results/e8c3a718b01a08b62745a5935e902df9/chr_3.zip -O $@

imputed/chr4.zip:
	wget https://imputationserver.sph.umich.edu/share/results/df1fe2f7e9c38833bb8101cc220ec117/chr_4.zip -O $@

imputed/chr5.zip:
	wget https://imputationserver.sph.umich.edu/share/results/f5028d386837046937d7a602819d87b1/chr_5.zip -O $@

imputed/chr6.zip:
	wget https://imputationserver.sph.umich.edu/share/results/1809c61463502314588a6e5a716d00ed/chr_6.zip -O $@

imputed/chr7.zip:
	wget https://imputationserver.sph.umich.edu/share/results/569ea25a88a582662dc8df6220640a8a/chr_7.zip -O $@

imputed/chr8.zip:
	wget https://imputationserver.sph.umich.edu/share/results/c72df023c8928d73ec439a887fbcd08f/chr_8.zip -O $@

imputed/chr9.zip:
	wget https://imputationserver.sph.umich.edu/share/results/2ae8e4a4d1100ab4cec8f2aedd715250/chr_9.zip -O $@

################################################################
UNZIP := $(foreach chr, $(shell seq 1 22), imputed/chr$(chr).vcf.gz)

unzip: $(UNZIP)

imputed/chr%.vcf.gz: imputed/chr%.zip
	unzip $< -d $(dir $@)

################################################################
## Q/C and generate PLINK
getPlink : $(foreach chr, $(shell seq 1 22), plink/chr${chr}.bed)


## Rename snp ID to avoid duplicate var problem
uniq/chr%.vcf.gz: imputed/chr%.dose.vcf.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	bcftools norm -m-any imputed/chr$*.dose.vcf.gz | bcftools norm --check-ref w -f human_g1k_v37.fasta | bcftools norm --rm-dup any | bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' -o $@

PLINK = ./bin/plink

## build plink with Q/C
plink/chr%.bed: uniq/chr%.vcf.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)

	$(PLINK) --vcf $< --double-id --snps-only --make-bed --out plink/raw_chr$*

	$(PLINK) --bfile plink/raw_chr$* --double-id --list-duplicate-vars ids-only --out plink/temp_chr$*

	$(PLINK) --bfile plink/raw_chr$* --exclude plink/temp_chr$*.dupvar --make-bed --out plink/nodup_chr$*

	cat imputed/chr$*.info.gz | gzip -d | tail -n+2 | awk -F '\t' '{ rs = $$1 ":" $$2 ":" $$3; scores[rs] = (rs in scores) ? (scores[rs] > $$7 ? $$7 : scores[rs]) : $$7; } END { for(rs in scores) print rs FS scores[rs] }' > plink/chr$*.info

	$(PLINK) --bfile plink/nodup_chr$* --qual-scores plink/chr$*.info 2 1 0 --qual-threshold 0.6 --make-bed --maf 0.05 --geno 0.1 --hwe 1e-4 --out plink/chr$*
	rm plink/temp_chr$** plink/nodup_chr$** plink/raw_chr$** plink/chr$*.info
