# Pulling out Dx-Relapse matches from MAF file

load("2019-02-14-allpdx-clean-maf-240.rda")
inData <- pptc.merge
#colnames(inData)


info_1 <- readLines("all_models.txt")
for(x in info_1)
{
  line = strsplit(x,"\t")[[1]]
  #print(line)
  
  file_name = paste0(line[1],"_",line[4],".txt")
  
  new <- subset(inData,TSB == line[3] | TSB == line[6] ,c(TSB,Hugo_Symbol,Variant_Classification,cDNA_Change,Protein_Change,VAF))
  write.table(new,file_name,sep="\t",row.names=FALSE)
  
  tmp <- read.table(file_name,sep="\t", header = TRUE)
  gsub("","0.0",tmp$Protein_Change)
  
  tmp$Phase <- ifelse(tmp$TSB == line[3],"Diagnosis","Relapse")
  write.table(tmp,file_name,sep="\t",row.names=FALSE)
  
}

# For Same Phase Samples
info <- readLines("same_phase_all_models.txt")
for(x in info)
{
  line = strsplit(x,"\t")[[1]]
  #print(line)
  
  file_name = paste0(line[1],"_",line[4],".txt")
  
  new <- subset(inData,TSB == line[3] | TSB == line[6] ,c(TSB,Hugo_Symbol,Variant_Classification,cDNA_Change,Protein_Change,VAF))
  write.table(new,file_name,sep="\t",row.names=FALSE)
  
  tmp <- read.table(file_name,sep="\t", header = TRUE)
  
  tmp$Phase <- ifelse(tmp$TSB == line[3],"Tumor Sample 1","Tumor Sample 2")
  write.table(tmp,file_name,sep="\t",row.names=FALSE)
  
  
}

# for ALL-105/115 
new <- subset(inData,TSB == "PPTC-AF03-XTP1-A-1-0-D-human" | TSB == "PPTC-AF08-XTP1-A-1-0-D-human" ,c(TSB,Hugo_Symbol,Variant_Classification,cDNA_Change,Protein_Change,VAF))
write.table(new,"ALL-105_ALL-115.txt",sep="\t",row.names=FALSE)

tmp <- read.table("ALL-105_ALL-115.txt",sep="\t", header = TRUE)

tmp$Phase <- ifelse(tmp$TSB == "PPTC-AF03-XTP1-A-1-0-D-human","Relapse 1","Relapse 3")
write.table(tmp,"ALL-105_ALL-115.txt",sep="\t",row.names=FALSE)
