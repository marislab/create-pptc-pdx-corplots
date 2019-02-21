#!/usr/bin/Python

import pandas as pd
import numpy as np
import glob
from os.path import basename

# for Dx-Relapse pairs
dx_rel_files = glob.glob("./Dx-Relapse/*.txt")

for x in dx_rel_files:
    base = basename(x).split(".")[0]
    #print(base)
    res = "./Dx-Relapse/"+base+"-vaf.txt"
    
    df = pd.read_csv(x, sep = "\t")
    df.fillna(0, inplace = True)

    # pivot table using protein change as index
    new_df = pd.pivot_table(df,index = ['Protein_Change','Hugo_Symbol'],columns = ['Phase'],values = 'VAF').reset_index().rename_axis(None,1)

    new_df.fillna(0, inplace = True)

    # Defining thresholds
    def f(row):
        if ( (row['Diagnosis'] > 0.10) & (row['Relapse'] > 0.10) ):
            val = 'common'
        elif ( (row['Diagnosis'] == 0.00) & (row['Relapse'] > 0.10) ):
            val = 'relapse-specific'
        else:
            val = 'early-specific'
        return val

    new_df['group'] = new_df.apply(f, axis=1)

    new_df = pd.merge(new_df, df, how='left')

    # adding columns to final dataframe that will be written to the file
    final = new_df.filter(['Hugo_Symbol','Diagnosis','Relapse','group','Variant Classification','Protein_Change','cDNA_Change'], axis=1)
    final = final.drop_duplicates()
    final.fillna(0, inplace = True)

    final.to_csv(res, sep = "\t")

print("Done for Dx-Relapse Samples!")



# for pairs in Same-Phase
same_phase_files = glob.glob("./Same-Phase/*.txt")

for x in same_phase_files:
    if(x == "./Same-Phase/ALL-105_ALL-115.txt"):

        base = basename(x).split(".")[0]
        
        res = "./Same-Phase/"+base+"-vaf.txt"
        df = pd.read_csv(x, sep = "\t")
        df.fillna(0, inplace = True)


        # pivot table using protein change as index
        new_df = pd.pivot_table(df,index = ['Protein_Change','Hugo_Symbol'],columns = ['Phase'],values = 'VAF').reset_index().rename_axis(None,1)

        new_df.fillna(0, inplace = True)

        # Defining thresholds
        def f(row):
            if ( (row['Relapse 1'] > 0.10) & (row['Relapse 3'] > 0.10) ):
                val = 'common'
            elif ( (row['Relapse 1'] == 0.00) & (row['Relapse 3'] > 0.10) ):
                val = 'relapse-3-specific'
            else:
                val = 'relapse-1-specific'
            return val

        new_df['group'] = new_df.apply(f, axis=1)


        new_df = pd.merge(new_df, df, how='left')
        #new_df

        # adding columns to final dataframe that will be written to the file
        final = new_df.filter(['Hugo_Symbol','Relapse 1','Relapse 3','group','Variant Classification','Protein_Change','cDNA_Change'], axis=1)
        final = final.drop_duplicates()
        final.fillna(0, inplace = True)

        final.to_csv("./Same-Phase/ALL-105_ALL-115-vaf.txt", sep = "\t")







    else:
    
        base = basename(x).split(".")[0]
        #print(base)
        res = "./Same-Phase/"+base+"-vaf.txt"
        
        df = pd.read_csv(x, sep = "\t")
        df.fillna(0, inplace = True)
        
        # pivot table using protein change as index
        new_df = pd.pivot_table(df,index = ['Protein_Change','Hugo_Symbol'],columns = ['Phase'],values = 'VAF').reset_index().rename_axis(None,1)

        new_df.fillna(0, inplace = True)

        # Defining thresholds
        def f(row):
            if ( (row['Tumor Sample 1'] > 0.10) & (row['Tumor Sample 2'] > 0.10) ):
                val = 'common'
            elif ( (row['Tumor Sample 1'] == 0.00) & (row['Tumor Sample 2'] > 0.10) ):
                val = 'Tumor-Sample-2-specific'
            else:
                val = 'Tumor-Sample-1-specific'
            return val

        new_df['group'] = new_df.apply(f, axis=1)


        new_df = pd.merge(new_df, df, how='left')
        #new_df

        # adding columns to final dataframe that will be written to the file
        final = new_df.filter(['Hugo_Symbol','Tumor Sample 1','Tumor Sample 2','group','Variant Classification','Protein_Change','cDNA_Change'], axis=1)
        final = final.drop_duplicates()
        final.fillna(0, inplace = True)
        final.to_csv(res, sep = "\t")

      

print("Done for Same-Phase Samples!")

