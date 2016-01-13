#!/bin/bash
# Requires vcf file as input
# $1 = vcf file, $2 = output file, $3 = ped file

PSAP_PATH= #ADD PATH TO PSAP DIRECTORY HERE
ANNOVAR_PATH= # ADD PATH TO ANNOVAR DIRECTORY HERE
echo $PWD

if [ $# -gt 0 ] && [ $1 == "-h" ]
then
	echo "arg1 =  VCF file"
	echo "arg2 = output file name" 
	echo "arg3 = family pedigree file"
	echo "Example: popScore_analysis_pipeline.sh arg1 arg2 arg3"
	exit
fi

if [ $# == 3 ]
then  
# convert vcf file to annovar file
	echo "Converting VCF file to annovar input"
	perl ${ANNOVAR_PATH}/convert2annovar.pl -format vcf4old $1 -outfile ${2}.avinput -includeinfo

# write column names from vcf file to header file (will be used later)
	grep '#' $1 | tail -n 1 > ${2}.avinput.header

# if there is no annotated directory create annotated directory
	if [ $(ls -d $PWD/*/ | grep -c -w "annotated") == 0 ]
	then
		echo "Creating directory annotated/ to store annotation and analysis output"
		mkdir annotated
	fi

# check that all required ANNOVAR annotation files are in the humandb directory - if not download them
	MISSING=0
	for FILE in "hg19_ALL.sites.2014_09.txt" "hg19_cadd.txt" "hg19_esp6500si_all.txt" "hg19_snp137.txt" "hg19_wgEncodeGencodeBasicV19Mrna.fa" "hg19_wgEncodeGencodeBasicV19.txt" "hg19_mac63kFreq_ALL.txt"
	do	
		if [ ! -f ${ANNOVAR_PATH}/humandb/$FILE ]
		then
			MISSING=$(( $MISSING+1 ))
		fi
	done
	
	if [ $MISSING -gt 0 ]
	then
		echo "Missing required ANNOVAR annotations.  Downloading all necessary files and adding them to the ANNOVAR annotation directory."
		bash ${PSAP_PATH}/psap/get_annovar_annos.sh
	fi
	
# annotate with annovar
	echo "Annotating data with ANNOVAR"
	perl ${ANNOVAR_PATH}/table_annovar.pl ${2}.avinput -remove -outfile annotated/${2}.avinput ${ANNOVAR_PATH}/humandb/ -buildver hg19 -protocol wgEncodeGencodeBasicV19,mac63kFreq_ALL,esp6500si_all,1000g2014sep_all,snp137,cadd -operation g,f,f,f,f,f -nastring NA -otherinfo -argument -separate,,,,,-otherinfo

# annotate with popstat (requires ped file)
	bash ${PSAP_PATH}/psap/annotate_PSAP.sh ${2}.avinput $3 $PSAP_PATH

else
	echo $# "arguments provided"
	echo "Incorrect number of arguments"
	echo "Please provide VCF file, output file name, family pedigree file"
fi
