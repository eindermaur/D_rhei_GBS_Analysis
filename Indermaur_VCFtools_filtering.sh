### Elizabeth J Indermaur ###
### VCFtools initial filtering for Didymella rhei ###

### VCFtools documentation found here: ###
# https://vcftools.sourceforge.net/man_latest.html

# Filtering parameters for D. rhei
# MAF > 0.01
# remove indels
# 5 < mean site depth < 50
# set genotypes with <5x coverage to missing
# missing data <50%

# Call vcftools
/programs/vcftools-v0.1.16/bin/vcftools 

# Apply filters
/programs/vcftools-v0.1.16/bin/vcftools \
--vcf /home/$USER/YOUR_DIRECTORY/YOUR_FILE.vcf \
--maf 0.01 \
--minDP 5 \
--min-meanDP 5 \
--max-meanDP 50 \
--max-missing 0.5 \
--remove-indels \
--out /home/$USER/YOUR_DIRECTORY/YOUR_NEW_FILE.vcf \
--recode

# Between these steps, Tassel was used to inspect new dataset and 
# remove individuals that have >50% missing data

# vcf recoding for -012 format
# This creates necessary .012, .012.pos, and .012.indv files for 
# identity by state matrix ahead of clone correction 
/programs/vcftools-v0.1.16/bin/vcftools \
--vcf /home/$USER/YOUR_DIRECTORY/YOUR_FILE.vcf \
--012 \
--out /home/$USER/YOUR_DIRECTORY/YOUR_NEW_FILE

# Create missing data files to select genotypes with the least
# amount of missing data
# This creates necessary .imiss file for clone correction
/programs/vcftools-v0.1.16/bin/vcftools \
--vcf /home/$USER/YOUR_DIRECTORY/YOUR_FILE.vcf \
--missing-indv \
--out /home/$USER/YOUR_DIRECTORY/YOUR_NEW_FILE

### ALL DONE ###


