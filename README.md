## CONTENTS
THE PSAP PACKAGE CONTAINS FOUR BASH SCRIPTS AND ALL DEPENDANT R SCRIPTS AND LOOKUP TABLES. EACH SCRIPT IS DESCRIBED BELOW.

1) ```family_analysis_pipeline.sh```: Calls ANNOVAR to annotate data, calls an Rscript that performs some basic cleaning steps (mendelian inheritance filter - allows de novos, PSAP calibration filter, missing data filter, allele frequency discrepancy filter) and annotates all individuals with PSAP, calls an R script that performs a family based analysis to identify candidate variants (shared among all affected individuals and pattern of inheritance is consistent with disease model) 

2) ```individual_analysis_pipeline.sh```: Calls ANNOVAR to annotate data, calls an Rscript that performs some basic cleaning steps (mendelian inheritance filter - allows de novos, PSAP calibration filter, missing data filter, allele frequency discrepancy filter) and annotates all individuals with PSAP, and calls an R script that will report out candiate variants (inheritance pattern consistent with disease model)

3) ```annotate_psap.sh```: Assumes data has already been annotated by ANNOVAR, calls an Rscript that performs some basic cleaning steps (mendelian inheritance filter - allows de novos, PSAP calibration filter, missing data filter, allele frequency discrepancy filter) and annotates all individuals with PSAP, and calls an R script that performs a family based analysis to identify candidate variants.  This script is not a stand alone script and will only run when called by one of the above scripts.

4) ```get_annovar_annos.sh```:  Assumes ANNOVAR has been downloaded.  It downloads all required annotation files from ANNOVAR.  This script only needs to be run once and should be run before running any of the above scripts.  This script requires an internet connection to run successfully.

## REQUIRED SOFTWARE
This pipeline uses the R statistical software and ANNOVAR.  Please ensure R (http://r-project.org) and ANNOVAR (http://annovar.openbioinformatics.org) are installed.  Paths to all other accessory softwares/scripts are hard coded to the directories within the PSAP directory.

##### NOTE: Using other annotation software is not recommended because annotations may differ from the PSAP lookup tables and will introduce biases in downstream analyses. 

## PREPARING THE SCRIPT
The local paths to the PSAP directory and ANNOVAR software need to be hard coded into the ```get_annovar_annos.sh```, ```family_analysis_pipeline.sh``` and ```individual_analysis_piepline.sh`` scripts at the lines indicated within each script.

Prior to running any of the analysis or annotation scripts for the first time you will need to run the get_annovar_annos.sh script to download all necessary ANNOVAR annotation files.  
##### NOTE: It may take this script several days to download all of the necessary annotation files as some of them are very large (80 GB)

## RUNING THE SCRIPT
To run the ```family_analysis_pipeline.sh``` or ```individual_analysis_pipeline.sh``` a VCF file, output file name, and pedigree file must be provided as arguments to the script in that order.

EG. ```family_analysis_pipeline.sh example.VCF example example.ped OR individual_analysis_pipeline.sh example.VCF example example.ped```

##### NOTE: This script requires the ```FILENAME.avinput.hg19_multianno.txt``` to be present and it requires the ```FILENAME.avinput.hg19_multianno.txt``` file be annotated with GencodeV19 gene names, CADD scored, and allele frequencies from Sep 2014 release of 1000 Genomes, ESP 6500, and Mac63k_Freq (ExAC)

The ```-h``` argument will provide a list of all the necessary arguments and print an example of the syntax for running the analysis scripts
If the incorrect number of arguments is provided (too many or too few) the script will exit with an "INCORRECT NUMBER OF ARGUMENTS" message and a list of the necessary arguments.

We recommend running this script with at least 4GB of memory and the user should have at least 500GB of disk space available.

## DEPENDANCY GENERATED FILES
```FILENAME.avinput```: This file is created by ANNOVAR from the provided VCF within the pipeline scripts.  The information from the VCF is used to create an annovar formatted file that ANNOVAR will use to annotate with ANNOVAR annotations.

```FILENAME.avinput.hg19_multianno.txt```: This file is created by ANNOVAR from the FILENAME.avinput file within the pipeline scripts.  This is the final ANNOVAR output file and contains all the ANNOVAR annotations required for PSAP analysis and candidate variant identification. More information for this file can be found in the PSAP_OUTPUT_GUIDE.

## SCRIPT GENERATED FILES
```FILENAME.header```: This file contains the VCF header and is created within the pipeline scripts.  It is used to identify data for each of the individuals included in the analysis

```FILENAME_popScore.txt```: This file contains all the data that is annotated with PSAP for all individuals if using the family analysis pipeline or one per individual if using the individual based pipeline.  More information for this file can be found in the PSAP_OUTPUT_GUIDE.

```FILENAME_missing.txt```: This file contains all the data that cannot be annotated with PSAP for all individuals if using the family analysis pipeline or one per individual if using the individual based pipeline

```FILE.report.txt```: This file contains all candidate variants, ordered by PSAP value with the best candidate at the top.  More information for this file can be found in the PSAP_OUTPUT_GUIDE.

### PEDIGREE FILE FORMAT (SPACE SEPARATED, NO HEADER):

```
FAMILY ID
INDIVIDUAL ID
PATERNAL ID (0 IF NO FATHER)
MATERNAL ID (0 IF NO MOTHER)
GENDER (1 FOR MALE, 2 FOR FEMALE)
CASE-CONTROL STATUS (1 FOR UNAFFECTED, 2 FOR AFFECTED)
```
##### NOTE: The individual ID must match the ID used in the VCF header for that individual

EXAMPLE PED FILE FOR A FAMILY:  

```
FAM1 INDV1 PAT1 MAT1 1 2 
FAM1 PAT1 0 0 1 1 
FAM1 MAT1 0 0 2 1 
```

EXAMPLE PED FILE FOR AN INDIVIDUAL: 

```
INDV1 INDV1 0 0 2 2 
```

VCF files with multiple unrelated individuals can be analysed using the ```individual_analysis_pipeline.sh```.  The individuals do not need to be split into separate VCF files and all individuals shoudl be included in a single pedigree file.
  
For best results, all individuals or family members should be included in a single multi-person VCF file and should be jointly genotyped.

This script will write all results to the directory in which it is called.  
##### Please ensure the script is called from a directory in which you have permission to write files.  When the script finishes running it will print the directory containing all results.



```
$ANNOVAR_PATH/annovar_latest/annotate_variation.pl
$ANNOVAR_PATH/annovar_latest/convert2annovar.pl
$ANNOVAR_PATH/annovar_latest/table_annovar.pl
$ANNOVAR_PATH/annovar_latest/humandb/
$PATH_PATH/psap/RScripts/generic_apply_popStat.R
$PATH_PATH/psap/RScripts/generic_candidate_analysis.R
$PATH_PATH/psap/RScripts/individual_apply_popstat.R
$PATH_PATH/psap/RScripts/unrelated_candidate_analysis.R
```

## WHEN USING THIS SCRIPT PLEASE CITE
Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data. Nucleic Acids Research, 38:e164, 2010

Wilfert AB, Chao K, Kaushal M, Jain S, Zöllner S, Adams DR and Conrad DF.  Genome-wide significance testing of variation from single case exomes. Nature Genetics. doi:10.1038/ng.3697. 2016
