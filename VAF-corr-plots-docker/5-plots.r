# .libPaths("/Users/patelk26/anaconda3/lib/R/library") # Change lib paths according to your Terminal lib path, that's where all the packages associated with the script will be installed.
set.seed(12345)

# Checking for required packages, if not found, will be installed

if(!require(ggplot2)){
  install.packages("ggplot2", repos='http://cran.us.r-project.org')
}
if(!require(ggbeeswarm)){
  install.packages("ggbeeswarm", repos='http://cran.us.r-project.org')
}
if(!require(gghighlight)){
  install.packages("gghighlight", repos='http://cran.us.r-project.org')
}
if(!require(ggrepel)){
  install.packages("ggrepel", repos='http://cran.us.r-project.org')
}
if(!require(ggthemes)){
  install.packages("ggthemes", repos='http://cran.us.r-project.org')
}
if(!require(ggpubr)){
  install.packages("ggpubr", repos='http://cran.us.r-project.org')
}



pptc.folder <- "~/Box Sync/PPTC-genomics-collaboration/"
source(paste0(pptc.folder, "Manuscript/scripts/theme.R"))

# To Fetch Histologies from Histology.Detailed column 
df_hist <- read.table("~/Box Sync/PPTC-genomics-collaboration/Data/clinical/2019-06-21-pdx-clinical-final-for-paper.txt", sep = "\t", header = T)


# Storing all Dx-Relapse Sample pairs in a list so each file can be processed iteratively
file_list <- list.files(path="./Dx-Relapse/", pattern="*-vaf.txt") 


for (x in file_list){
  
  base = strsplit(x,"-vaf.txt")[[1]]  # splitting to save Dx and relapse sample names
  words = strsplit(base,"_")[[1]] 
  
  axis_x = paste0(words[1]) # Model name 
  axis_y = paste0(words[2]) # Model name
  
  # Reading files iteratively and storing in a dataframe
  dat1 <- read.table(file = paste0("./Dx-Relapse/",x), header=TRUE, sep = "\t")
  
  dat1[is.na(dat1)] = 0  # Replace all na's with 0
  dat1$Protein_Change[dat1$Protein_Change== "0"] <- " "
  
  # Matching the model name and fetching corresponding Histology.Detail for plot_title
  inx = which(df_hist$Model %in% axis_x)
  plot_title = df_hist$Histology.Detailed[inx]
  
  # Plotting scatterplot for all Dx-Relapse Samples 
  
  p <- ggplot(dat1, aes(dat1$Diagnosis, dat1$Relapse, color = group)) + 
    geom_point( size = 10,fill = 4, alpha = 1/6) + 
    scale_colour_manual(values = c("gray34", "dodgerblue3", "firebrick3")) + 
    labs(title = paste(plot_title), x = axis_x, y = axis_y) + 
    geom_vline(xintercept = 0.1, linetype = "dashed") + # Adding vertical intercept
    geom_hline(yintercept = 0.1, linetype = "dashed") + # Adding horizontal intercept line
    geom_text_repel(aes(label =ifelse(dat1$label > 0, paste("",dat1$Hugo_Symbol,"",dat1$label_PC,""),'')),size = 3.5,hjust = 0,vjust = 0, nudge_x = 0.005,point.padding = NA,segment.color = NA, show.legend = FALSE, xlim = c(0.02,NA),ylim = c(0.025,0.96)) + # if label column is one, then label the point with gene name and it's corresponding Protein Change value
    theme_Publication() + 
    theme(plot.title = element_text(hjust = 0.5))+ 
    xlim(0,1) + 
    ylim(0,1)
  
  # Save as PDF
  ggsave(paste0("./plots/",Sys.Date(),"-",x,"-vaf.pdf"), p, width=8, height=6, device = "pdf")
  
}


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# For all Same Phase Samples
file_list <- list.files(path="./Same-Phase/", pattern="*-vaf.txt") 

for (x in file_list){
 
  
  if (x == "ALL-105_ALL-115-vaf.txt"){
    
    dat1 <- read.table(file = paste0("./Same-Phase/",x), header=TRUE,sep = "\t")
    
    dat1[is.na(dat1)] = 0
    dat1$Protein_Change[dat1$Protein_Change== "0"] <- " "
    
    # Plotting scatterplot for ALL-105 and ALL-115
    
    p <- ggplot(dat1, aes(dat1$Relapse.1, dat1$Relapse.3, color = group)) + 
      geom_point( size = 10,fill = 4, alpha = 1/8) + 
      scale_colour_manual(values = c("gray34", "firebrick3", "firebrick3")) +
      labs(title = "BCP-ALL", x = "ALL-105", y = "ALL-115") + 
      geom_vline(xintercept = 0.1, linetype = "dashed")+ 
      geom_hline(yintercept = 0.1, linetype = "dashed") + 
      geom_text_repel(aes(label =ifelse(dat1$label > 0, paste("",dat1$Hugo_Symbol,"",dat1$label_PC,""),'')),size = 3.5,hjust = 0,vjust = 0, nudge_x = 0.005,point.padding = NA,segment.color = NA, show.legend = FALSE, xlim = c(0.02,NA),ylim = c(0.025,0.96)) + 
      theme_Publication() + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      xlim(0,1) + 
      ylim(0,1)

    # Save as PDF
    ggsave(paste0("./plots/",Sys.Date(),"-same-phase-ALL-105-115-vaf.pdf"), p, width=8, height=6, device = "pdf")
 
  }
  
  else if (x == "OS-34_OS-34-SJ-vaf.txt"){
    
    dat1 <- read.table(file = paste0("./Same-Phase/",x), header=TRUE,sep = "\t")
    
    dat1[is.na(dat1)] = 0
    dat1$Protein_Change[dat1$Protein_Change== "0"] <- " "
    
    # Plotting Scatterplot for OS-34 and OS-34-SJ
    
    p <- ggplot(dat1, aes(dat1$Tumor.Sample.1, dat1$Tumor.Sample.2, color = group)) + 
      geom_point( size = 10,fill = 4, alpha = 1/8) + 
      scale_colour_manual(values = c("gray34", "dodgerblue3", "dodgerblue3"))
    p <- p + labs(title = "Osteosarcoma", x = "OS-34", y = "OS-34-SJ") + 
      geom_vline(xintercept = 0.1, linetype = "dashed")+ 
      geom_hline(yintercept = 0.1, linetype = "dashed") + 
      geom_text_repel(aes(label =ifelse(dat1$label > 0, paste("",dat1$Hugo_Symbol,"",dat1$label_PC,""),'')),size = 3.5,hjust = 0,vjust = 0, nudge_x = 0.005,point.padding = NA,segment.color = NA, show.legend = FALSE, xlim = c(0.02,NA),ylim = c(0.025,0.96)) + 
      theme_Publication() + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      xlim(0,1) + 
      ylim(0,1)
    
    # Save as PDF
    ggsave(paste0("./plots/",Sys.Date(),"-same-phase-Relapse-OS-34-OS-34-SJ-vaf.pdf"), p, width=8, height=6, device = "pdf")
    
    
    }# end of else if block
  
  else {
  
  
  base = strsplit(x,"-vaf.txt")[[1]]
  words = strsplit(base, "_")[[1]] 
  axis_x = paste0(words[1])
  axis_y = paste0(words[2])
  

  dat1 <- read.table(file = paste0("./Same-Phase/",x), header=TRUE,sep = "\t")
  
  dat1[is.na(dat1)] = 0
  dat1$Protein_Change[dat1$Protein_Change== "0"] <- " "
  
  # Saving plot-title from Histology.Detailed column to the corresponding Model
  inx = which(df_hist$Model %in% axis_x)
  plot_title = df_hist$Histology.Detailed[inx]
  
  # Plotting ScatterPlot for the rest of the Same Phase Samples
  
  p <- ggplot(dat1, aes(dat1$Tumor.Sample.1, dat1$Tumor.Sample.2, color = group)) + 
    geom_point( size = 10,fill = 4, alpha = 1/6) + 
    scale_colour_manual(values = c("gray34", "firebrick3", "firebrick3")) + 
    labs(title = paste(plot_title), x = axis_x, y = axis_y) + 
    geom_vline(xintercept = 0.1, linetype = "dashed")+ 
    geom_hline(yintercept = 0.1, linetype = "dashed") + 
    geom_text_repel(aes(label =ifelse(dat1$label > 0, paste("",dat1$Hugo_Symbol,"",dat1$label_PC,""),'')),size = 3.5,hjust = 0,vjust = 0, nudge_x = 0.005,point.padding = NA,segment.color = NA, show.legend = FALSE, xlim = c(0.02,NA),ylim = c(0.025,0.96)) + 
    theme_Publication() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlim(0,1) + 
    ylim(0,1)
  
  # Save as PDF
  ggsave(paste0("./plots/",Sys.Date(),"-same-phase-",x,"-vaf.pdf"), p, width=8, height=6, device = "pdf")

  } # end of else block

  
}


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Total mutations for OS-36-SJ/OS-36 (barplot)
dat2 <- read.delim("barplot_OS-36-SJ-36.csv", header=TRUE)

dat2[is.na(dat2)] = 0

q <-ggplot(dat2, aes(x= factor(dat2$Model, levels = c("OS-36-SJ","OS-36")), y=dat2$Total.Mutation, fill = group)) +
  geom_bar(stat = "identity",width = 0.3, position = "dodge") + 
  scale_fill_manual(values = c("firebrick3","firebrick3")) + 
  labs(title = "Osteosarcoma", x = "Model", y = "Total Mutations") + 
  theme_Publication() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.title = element_blank())
# Save as PDF
ggsave(paste0("./plots/",Sys.Date(),"-OS-36-36SJ-tot_mut-barplot.pdf"), q, width=15, height=15, device = "pdf", units = "cm")


# Total mutations for ALL-105/115 (barplot)
data <- read.csv("barplot_ALL-102-105-115.csv", header=TRUE)

data[is.na(data)] = 0

# Plotting Barplot to compare total # of mutations between ALL-102, ALL-105 and ALL-115

r <- ggplot(data, aes(x=data$Model, y=data$Total.Mutation, fill = group)) +
  geom_bar( stat="identity",width = 0.3) + 
  scale_fill_manual(values = c("dodgerblue3", "firebrick3", "gray")) + 
  labs(title = "", x = "Models", y = "Total Mutations") + 
  theme_Publication() + 
  theme(plot.title = element_text(hjust = 0.5))

# Save as PDF
ggsave(paste0("./plots/",Sys.Date(),"-ALL-105-115-tot_mut-barplot.pdf"), r, width=5, height=5, device = "pdf")


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Mutations across all Dx-Relapse Samples (Stacked barplot)
dat3 <- read.csv("./total_Mutations-All-Dx-Relapse-Models.csv", header=TRUE)

dat3[is.na(dat3)] = 0

dat3$Phase <- factor(dat3$Phase, levels = c("Relapse 3","Relapse 1","Diagnosis"))

# Plotting Stacked-barplot for all Dx-Relapse samples

s <- ggplot(dat3, aes(x = factor(dat3$Flag), y = Total.Mutations,fill = Phase))+  
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("gray", "firebrick3", "dodgerblue3")) +
  geom_bar( stat="identity", width = 0.9) + 
  scale_x_discrete(labels= c('ALL-46 \n ALL-121','ALL-58 \n ALL-123','ALL-80 \n ALL-81','ALL-102 \n ALL-105 \n ALL-115','ALL-25 \n ALL-61','ALL-82 \n ALL-83','COG-N-603x \n COG-N-623x','ICb-9850PNET \n IC-22909PNET-rIII','NCH-CA-1 \n NCH-CA-2','OS-32 \n OS-36','OS-32 \n OS-36-SJ','Rh-30 \n Rh-30R')) +
  theme_Publication() + 
  theme(axis.text.x = element_text(angle = 85, hjust = 1, vjust =1)) + 
  labs(title = "", x = "Models", y = "Total Mutations") +
  ylim(0,1000)
ggsave(paste0("./plots/",Sys.Date(),"-stacked-barplot-Dx-Relapse-AllSamples.pdf"), s, width=10, height=8, device = "pdf")

dat3$logMut <- log(dat3$Total.Mutations,2)
dat3$facet <- "Matched Pairs"
v <- print(ggviolin(dat3, x = "Phase_1", y = "logMut", fill = "Phase_1",
               palette = dxrelcol, alpha = 0.8,add = "boxplot",
               add.params = list(fill = "white"))+
        stat_compare_means(comparisons = list(c("Diagnosis", "Relapse")), label.y = c(15), label = "p.format")+ # Add significance levels
        #stat_compare_means(label.y.npc = "top") + ##global p
        facet_grid(~facet) +##global p
        theme_Publication() + 
        xlab("") + ylab('log2[Total Mutations]')+
        scale_y_continuous(limits=c(0,20)))
ggsave(paste0("./plots/",Sys.Date(),"-boxplot-Dx-Relapse-AllSamples.pdf"), v, width=6, height=6, device = "pdf")



