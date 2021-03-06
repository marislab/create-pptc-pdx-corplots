#!/bin/bash
# Following files needs to be present in the same folder as this script:
# 1. 2019-02-14-allpdx-clean-maf-240.rda (maf file)
# 2. 2019-02-20-pdx-clinical-for-web.txt (clinical file)
# 3. all_models.txt
# 4. same_phase_all_models.txt
# 5. barplot_OS-36-SJ-36.txt (for total mutations for OS - barplot)
# 6. total_Mutations-All-Dx-Relapse-Models.csv (for total mutations across all samples - stacked barplot)
# 7. barplot_ALL-102-105-115.csv (for total mutations for ALL-105/115 - barplot)
# 8. A folder named final-gene-lists containing all final gene lists
# 9. Scripts 1-5 in the same folder

# NOTE: change .libPaths() in the 5-plots.r script to the respective libPath of the terminal
# Requirement: Python 3.6
# Usage: ./pipeline.sh

echo "Fetching data files..."
mkdir data/ final-gene-lists/
cd data/

# 1. Download DNA MAF file
wget --output-document='2019-02-14-allpdx-clean-maf-240.rda' https://ndownloader.figshare.com/files/14414198

# 2. Clinical file
wget --output-document='pptc-pdx-clinical-web.txt' https://ndownloader.figshare.com/files/14508536

cd ../final-gene-lists/
wget --output-document='leukemia-goi-list.txt' https://ndownloader.figshare.com/files/14541185
wget --output-document='osteosarcoma-goi-list.txt' https://ndownloader.figshare.com/files/14541191
wget --output-document='renal-goi-list.txt' https://ndownloader.figshare.com/files/14541197
wget --output-document='brain-goi-list.txt' https://ndownloader.figshare.com/files/14541182
wget --output-document='neuroblastoma-goi-list.txt' https://ndownloader.figshare.com/files/14541188
wget --output-document='rare-goi-list.txt' https://ndownloader.figshare.com/files/14541194
wget --output-document='sarcoma-goi-list.txt' https://ndownloader.figshare.com/files/14541200

cd ..

echo "Pulling the data from MAF file..."

Rscript 1-pull-data.r

echo "Organizing files into Dx-Relapse and Same-phase folders..."

mkdir Dx-Relapse Same-Phase plots packages

python 2-organize-files.py

echo "Restructuring Data..."

python 3-restructure-data.py

echo "Labelling genes..."

python 4-labelling-genes.py

echo "Plotting samples..."

Rscript 5-plots.r

echo "Finished!"
