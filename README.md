# General Outline of SNP Discovery Workflow 

E.J. Indermaur, 2023

## Summary:

This takes you through each step from aligning your raw reads to analyzing population structure. It requires use of a cloud server and R (on your machine). 
You will need to generate some files on your own (see Step 5). All scripts can be found here, except for custom scripts for Step 1. Please email eji7@cornell.edu if you're interested. 
Each script may contain errors and should serve as a template that you can customize to best handle your data. 

For good examples of how to describe the methods, consider reading: 
- Greg Vogel’s [*P. capsici* GWAS paper](https://apsjournals.apsnet.org/doi/full/10.1094/PHYTO-04-20-0112-FI)
* Martha Sudermann’s [*P. fulva* paper](https://apsjournals.apsnet.org/doi/10.1094/PHYTO-06-21-0244-R)
+ Chase Crowell’s [*M. americana* paper](https://apsjournals.apsnet.org/doi/10.1094/PHYTO-05-21-0201-R)

Before you start, please consider viewing the following resources:
- Elshire et al. 2011, [original GBS publication](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019379)
* Avi Karn’s [SNP calling tutorial](https://avikarn.com/2019-04-20-GBS-SNP-calling-tutorial/)
+ The Buckler Lab’s [TASSEL documentation](https://www.maizegenetics.net/tassel) (also where you download the software)
- Cornell BioHPC [Software user guide](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=445), for information on TASSEL and VCFtools
* The [VCFtools github](https://vcftools.github.io/man_latest.html)

Happy genotyping!

## Workflow:

### 1.	Identify tags and taxa from `fastq` files, align reads to reference, call SNPs, converts to `.vcf`

- **Host:** server
- **Program:** Tassel
- **Script:** Please email eji7@cornell.edu for custom script
- **Input:** `fastq` file(s), key file, reference genome
- **Output:** your variants in `.vcf`

**NOTE:** 
- See Avi Karn’s website for instructions on how to make your key file. 
- If your reference genome is in `fastq.gz`, you’ll need to change the extension to `fastq.bgz`

### 2.	Initial filtering of .vcf

- **Host:** server
- **Program:** vcftools
- **Script:** [Indermaur_VCFtools_filtering.sh](https://github.com/eindermaur/D_rhei_GBS_Analysis/blob/main/Indermaur_VCFtools_filtering.sh)
- **Input:** your `.vcf` from Step 1
- **Output:** your filtered `.vcf`
  
**NOTE:** You can apply filters one by one to see which contribute to most of your site losses.

### 3.	Inspect .vcf and apply additional filters

- **Host:** your machine
- **Program:** Tassel5
- **Input:** your `.vcf` from Step 2
- **Output:** your filtered `.vcf` with any additional changes (e.g. removed individual taxa)

### 4.	Create files necessary for Identify by State (IBS) and clone correction

- **Host:** server
- **Program:** vcftools
- **Script:** [Indermaur_VCFtools_filtering.sh](https://github.com/eindermaur/D_rhei_GBS_Analysis/blob/main/Indermaur_VCFtools_filtering.sh)
- **Input:** your `.vcf` from Step 3
- **Output:** `.012`, `.012.pos`, `.012.indv`, and `.imiss` files

### 5.	Create 2 files with metadata necessary for clone correction and PCA

- **Host:** your machine
- **Program:** Excel

Here are example previews with the proper format:

1. `"NAME_isolates.csv"`
   
| Sample | SZ_for_analysis | In_storage |
| --- | --- | --- |
| 21_10_2 | 123456795 | TRUE |
| 21_10_3 | 123456796 | TRUE |

2. `"Isolate_metadata.csv"`

| SampleSZ | Source | year | Host_cultivar | State | Region | County | Field | UniqueGenotype | InCCDataset |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 21_10_2 | Libby Indermaur | 2021 | GVA_RHU_2021_1008 | NY | CNY | Ontario | Ontario1 | 1 | FALSE |
| 21_10_3 | Libby Indermaur | 2021 | GVA_RHU_2021_1008 | NY | CNY | Ontario | Ontario1 | 2 | TRUE |

### 6.	Calculate IBS  

- **Host:** your machine
- **Program:** R
- **Script:** [D_rhei_IBS_bothplates.R](https://github.com/eindermaur/D_rhei_GBS_Analysis/blob/main/D_rhei_IBS_bothplates.R)
- **Input:** `.012` files from Step 4
- **Output:** IBS matrix

### 7.	Clone correct

- **Host:** your machine
- **Program:** R
- **Script:** [D_rhei_IBS_bothplates.R](https://github.com/eindermaur/D_rhei_GBS_Analysis/blob/main/D_rhei_IBS_bothplates.R)
- **Input:** IBS matrix, `.012.indv`, `.imiss`, and `“NAME_isolates.csv”` from Step 5. 
- **Output:** `“Clone_assignments.txt”`, `“CC_samples_list.csv”` (a list of individuals with the least missing data per clonal group)

NOTE: 
1. At line 79, “x < .97” is the clonality cutoff I determined based on sequencing error rates between technical reps.
   This number may be different for you. To determine your cutoff, you will need to compare the error rates between your
   technical reps in the IBS matrix from Step 7. For example, if the identities between your reps are 96.7, 98.3, 95.8,
   and 98.6%, then your sequencing error rates are the inverse, or 3.3, 1.7, 4.2, and 1.4%. The mean sequencing error
   rate would be 2.5%. Because the “modify_matrix” function cannot take a thousandths place digit, your corresponding
   clonality cutoff would be “x < .97” instead of “x < .975”. This is an error rate cutoff of 3%.

### 8.	Create `.vcf` with ONLY the individuals in your `CC_samples_list.csv`

- **Host:** your machine
- **Program:** Tassel5
- **Input:** your `.vcf` from Step 3
- **Output:** `.vcf` with only the individuals representing unique genotypes

### 9.	Create `.012.indv` clone corrected files

- **Host:** BioHPC
- **Program:** vcftools
- **Script:** [Indermaur_VCFtools_filtering.sh](https://github.com/eindermaur/D_rhei_GBS_Analysis/blob/main/Indermaur_VCFtools_filtering.sh) (only step 2) 
- **Input:** your `.vcf` from Step 8
- **Output:** `.012`, `.012.pos`, `.012.indv` files for your CC dataset

### 10.	Summarize trends from clone correction

- **Host:** your machine
- **Program:** R
- **Script:** [PhenotypeAndCloneTreands.R](https://github.com/gmv23/pcap-gwas/blob/master/PhenotypeAndCloneTrends.R)
- **Input:** `Isolate_metadata.csv` from Step 5, `Clone_assignments.txt` from Step 7, and your clone corrected `.012.indv` file from Step 9
- **Output:** `“isolate_plotting_metadata.csv”`, `“field_plotting_metadata.csv”`, a summary file (not necessary for analysis but provides useful information), and a sample map

### 11.	Analyze population structure

- **Host:** your machine
- **Program:** R
- **Script:** [PopStructure.R](https://github.com/gmv23/pcap-gwas/blob/master/PopStructure.R)
- **Input:** `.012`, `.012.indv`, `.012.pos`, `“isolate_metadata.csv”`, `“isolate_plotting_metadata.csv”`, and `“field_plotting_metadata.csv”` from Step 10.
- **Output:** PC plots, table with Weir and Cockerham’s F<sub>st</sub> estimates

NOTE: 
1. This Fst estimate can only be calculated on populations with more than 5 unique genotypes. 
