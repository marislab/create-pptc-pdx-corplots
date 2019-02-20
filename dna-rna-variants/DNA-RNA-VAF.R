#!/usr/bin/env Rscript 

# Usage: ./DNA-RNA-VAF.R <DNA_MAF.rda> <RNA_MAF.maf> <clinical_file.txt>

if(!require(data.table)){
  install.packages("data.table", repos='http://cran.us.r-project.org')
}

args <- commandArgs(trailingOnly = TRUE)


# Loading dna.maf
print("loading dna maf...")
load(args[1]) # loads dna.maf as pptc.merge

# reading the RNA.MAF file
print("reading rna maf file...")
rna <- fread(args[2], header = T,sep="\t")

# Calculating RNA VAF values
rna.maf <- setDT(rna)[,RNA.VAF:= (t_alt_count / (t_ref_count + t_alt_count) )]

rm(rna)
gc()

# subsetting columns in the rna.maf dataframe to perform merge
rna.maf.subset <- subset(rna.maf,select = c("Chromosome", "Start_Position","HGVSc","Tumor_Sample_Barcode","RNA.VAF"))
rna.maf.subset <- setDT(rna.maf.subset)[,c("Tumor_Sample_Barcode_RNA") := tstrsplit(Tumor_Sample_Barcode, ".variants.maf", fixed=TRUE) ]
rna.maf.subset <- setDT(rna.maf.subset)[,c("RNA.human.bam.filename") := paste0(Tumor_Sample_Barcode_RNA, ".bam") ]
rm(rna.maf)
gc()

# reading clinical file to map models to RNA MAF
clinical <- fread(args[3], header = T, sep = "\t")
clinical.subset <- subset(clinical, select = c("RNA.human.bam.filename", "Model"))

# merging RNA MAF file with clinical file
print("performing merge with RNA MAF file...")
rna.model <- merge(rna.maf.subset,clinical.subset, all.x = T, by = "RNA.human.bam.filename")


# creating concat column in RNA and RNA MAF
pptc.merge <- setDT(pptc.merge)[,concat_col:= paste0(Chromosome,":",Start_position,":",cDNA_Change,":",Tumor_Sample_Barcode)]
rna.model <- setDT(rna.model)[,concat_col:= paste0(Chromosome,":",Start_Position,":",HGVSc,":",Model)]

# merging on common columns between RNA and DNA MAF 
print("performing merge with pptc-RNA-VEPpass.maf file...")
resulting_matches <- merge(pptc.merge,rna.model, all.x = T, by = "concat_col")


# Save new DNA MAF with RNA VAF information as RDA and MAF file
save(rna.model, file = paste0(Sys.Date(),"-allpdx-clean-maf-241.MAF.rda") )
write.table(rna.model, file = paste0(Sys.Date(),"-pptc-RNA-VEPpass-nonsilent.maf"), sep = "\t", row.names = F, quote = F, col.names = T)

