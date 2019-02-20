#!/usr/bin/Python
# labelling set of genes

import pandas as pd
import numpy as np
import glob
from os.path import basename

files1 = glob.glob("./Dx-Relapse/*-vaf.txt")

for x in files1:
    base = basename(x).split(".")[0]
    #print(base) # ALL-58_ALL-123-vaf
    
    # resulting file name with location
    res = "./Dx-Relapse/"+base+".txt"
    
    df1 = pd.read_table(x,sep = "\t")

    #----------- final gene lists --------------------#

    ALL = []
    with open("./final-gene-lists/leukemia-goi-list.txt") as fh1:
        for line in fh1:
            l = line.split()
            ALL.append(l[0])

    OS = []
    with open("./final-gene-lists/osteosarcoma-goi-list.txt") as fh2:
        for line in fh2:
            l = line.split()
            OS.append(l[0])

    IC = []
    with open("./final-gene-lists/brain-goi-list.txt") as fh3:
        for line in fh3:
            l = line.split()
            IC.append(l[0])

    NBL = []
    with open("./final-gene-lists/neuroblastoma-goi-list.txt") as fh4:
        for line in fh4:
            l = line.split()
            NBL.append(l[0])


    sarco = []
    with open("./final-gene-lists/sarcomacarcinoma-goi-list.txt") as fh5:
        for line in fh5:
            l = line.split()
            sarco.append(l[0])   


#----------- # Find if Hugo_Symbol present in the array of gene list --------------------#

    if(base.startswith("ALL")):

        df1.loc[df1['Hugo_Symbol'].isin(ALL)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(ALL),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")


    elif(base.startswith("OS")):

        df1.loc[df1['Hugo_Symbol'].isin(OS)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(OS),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")

    elif(base.startswith("COG")):

        df1.loc[df1['Hugo_Symbol'].isin(NBL)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(NBL),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")


    elif(base.startswith("IC")):
        
        df1.loc[df1['Hugo_Symbol'].isin(IC)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(IC),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")

    else:
        df1.loc[df1['Hugo_Symbol'].isin(sarco)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(sarco),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")



files2 = glob.glob("./Same-Phase/*-vaf.txt")

for x in files2:
    base = basename(x).split(".")[0]
    #print(base) # ALL-58_ALL-123-vaf
    
    # resulting file name with location
    res = "./Same-Phase/"+base+".txt"
    
    df1 = pd.read_table(x,sep = "\t")

    #----------- final gene lists --------------------#

    ALL = []
    with open("./final-gene-lists/leukemia-goi-list.txt") as fh1:
        for line in fh1:
            l = line.split()
            ALL.append(l[0])

    OS = []
    with open("./final-gene-lists/osteosarcoma-goi-list.txt") as fh2:
        for line in fh2:
            l = line.split()
            OS.append(l[0])

    IC = []
    with open("./final-gene-lists/brain-goi-list.txt") as fh3:
        for line in fh3:
            l = line.split()
            IC.append(l[0])

    NBL = []
    with open("./final-gene-lists/neuroblastoma-goi-list.txt") as fh4:
        for line in fh4:
            l = line.split()
            NBL.append(l[0])


    sarco = []
    with open("./final-gene-lists/sarcomacarcinoma-goi-list.txt") as fh5:
        for line in fh5:
            l = line.split()
            sarco.append(l[0])   



#----------- # Find if Hugo_Symbol present in the array of gene list --------------------#

    if(base.startswith("ALL")):

        df1.loc[df1['Hugo_Symbol'].isin(ALL)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(ALL),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")


    elif(base.startswith("OS")):

        df1.loc[df1['Hugo_Symbol'].isin(OS)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(OS),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")

    elif(base.startswith("COG")):

        df1.loc[df1['Hugo_Symbol'].isin(NBL)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(NBL),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")


    elif(base.startswith("IC")):
        
        df1.loc[df1['Hugo_Symbol'].isin(IC)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(IC),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")

    else:
        df1.loc[df1['Hugo_Symbol'].isin(sarco)]

        df1['label'] = 0
        df1['label_PC'] = ""

        # Set label as 1 at those indexes that matches 
        df1.loc[df1['Hugo_Symbol'].isin(sarco),'label'] = 1
        
        # set label_pc with corresponding protein change value
        df1['label_PC'] = np.where(df1['label'] == 1, df1['Protein_Change'].str[2:], df1['label_PC'])

        # Write to file
        df1.to_csv(res,index = False,sep = "\t")


print("Done labelling!")


