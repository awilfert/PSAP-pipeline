#!/bin/bash
# Requires vcf file as input
# $1 = vcf file, $2 = output file, $3 = ped file

PSAP_PATH= #INSERT PATH TO PSAP DIRECTORY HERE - INCLUDE THE PSAP DIRECTORY IN THE PATH
ANNOVAR_PATH= #INSERT PATH TO ANNOVAR DIRECTORY HERE - INCLUDE THE ANNOVAR DIRECTORY IN THE PATH
echo $PWD
module load R
if [ $# -gt 0 ] && [ $1 == "-h" ]
then
	echo "arg1 =  VCF file"
	echo "arg2 = output file name" 
	echo "arg3 = family pedigree file"
	echo "Example: popScore_analysis_pipeline.sh INDV001.vcf INDV.OUTFILE INDV001.ped"
	exit
fi

if [ $# == 3 ]
then
# Check that all required ANNOVAR annotation files are in the humandb directory
        MISSING=0
        for FILE in "hg19_ALL.sites.2014_09.txt" "hg19_cadd.txt" "hg19_esp6500si_all.txt" "hg19_snp137.txt" "hg19_wgEncodeGencodeBasicV19Mrna.fa" "hg19_wgEncodeGencodeBasicV19.txt" "hg19_mac63kFreq_ALL.txt"
        do
                if [ ! -f ${ANNOVAR_PATH}/humandb/$FILE ]
                then
                        MISSING=$(( $MISSING+1 ))
                fi
        done
# If any of the required annotation files are missing, exit and print error  message
        if [ $MISSING -gt 0 ]
        then
                echo "ERROR: Missing required ANNOVAR annotations.  Please run get_annovar_annos.sh prior to running this script."
                exit
        fi

# Extract and move to VCF file directory
        FILE_LOC=${1%/*.vcf} # Extract location of VCF file
        cd $FILE_LOC # Use location of  VCF file as working directory, this is where all output will be written
        echo $PWD
        VCF=${1##/*/} # Extract VCF file name
        OUTFILE=$2 # Name of output file (no directory should be included here)
        PED_FILE=$3 # Name of pedigree file (directory should be included here)

# Convert vcf file to annovar file
        echo "PROGRESS: Converting VCF file to annovar input"
        perl ${ANNOVAR_PATH}/convert2annovar.pl -format vcf4old $VCF -outfile ${OUTFILE}.avinput -includeinfo

# Write column names from VCF file to header file (will be used later)
        grep '#' $VCF | tail -n 1 > ${OUTFILE}.avinput.header # Extract all lines of the VCF header.  The last line of the VCF header contains coumn names - write columna names to .avinput.header file

# If there is no annotated directory create annotated directory
        if [ $(ls -d $PWD/*/ | grep -c -w "annotated") == 0 ]
        then
                echo "Creating directory annotated/ to store annotation and analysis output"
                mkdir annotated
        fi

# Annotate with ANNOVAR
	echo "PROGRESS: Annotating data with ANNOVAR"
	perl ${ANNOVAR_PATH}/table_annovar.pl ${OUTFILE}.avinput -remove -outfile annotated/${OUTFILE}.avinput ${ANNOVAR_PATH}/humandb/ -buildver hg19 -protocol wgEncodeGencodeBasicV19,mac63kFreq_ALL,esp6500si_all,1000g2014sep_all,snp137,cadd -operation g,f,f,f,f,f -nastring NA -otherinfo -argument -separate,,,,,-otherinfo

# EXTRACT INDIVIDUAL IDS
	echo "PROGRESS: Extracting individual IDs"
	IDS=$(awk '{print $2}' $PED_FILE)
	IDX=1

# RUN apply_popStat_individual.R for each individual
	echo "PROGRESS: Starting PSAP annotation" 
	for i in $IDS
	do
		Rscript ${PSAP_PATH}/RScripts/apply_popStat_individual.R ${OUTFILE}.avinput $i $PSAP_PATH &
		if [ `expr $IDX % 10` -eq 0 ]
		then
			echo "PROGRESS: Annotating individuals" $(( $IDX - 1 * 20)) "-" $(($IDX * 20)) "out of" ${#IDS}
			wait # Limit number of individuals annotated to no more than 20 at a time
		fi
		IDX=$(($IDX+1))
	done
	wait

# Generate report file - will look for variants present in two or more affected with PSAP < 1e-3 and not observed in unaffected
	echo "PROGRESS: Generating report file for all individuals"
	Rscript ${PSAP_PATH}/RScripts/unrelated_candidate_analysis.R ${OUTFILE}.avinput $PED_FILE $PSAP_PATH

else
	echo "ERROR: Incorrect number of arguments." $# "arguments provided"
	echo "Please provide a VCF file, output file name, and family pedigree file.  Please use the -h argument for help."
fi
