#!/bin/bash
# Data files for DNA-RNA-VAF.R

cd ~

# 1. Download DNA MAF file
wget --output-document='2019-02-14-allpdx-clean-maf-240.rda' https://ndownloader.figshare.com/files/14414198

# 2. Download RNA MAF file
#wget --output-document='' https://ndownloader.figshare.com/files/14414198

# 3. Clinical file
wget --output-document='pptc-pdx-clinical-web.txt' https://ndownloader.figshare.com/files/14508536