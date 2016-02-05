# $1 = annovar input file name (eg. FILE.avinput)
# $2 = PED file name (eg. FAMILY.ped)

PSAP_PATH=$3

if [ $# == 3 ]
then  
# annotate with popstat (requires ped file)
	echo "Annotating data with population model"
	Rscript ${PSAP_PATH}/psap/RScripts/generic_apply_popStat.R $1 ${2}.avinput $PSAP_PATH

# perform family analysis and validate
	echo "Performing family analysis and candidate validation"
	Rscript ${PSAP_PATH}/psap/RScripts/generic_candidate_analysis.R $1 ${2}.avinput $PSAP_PATH
	echo "Done.  All resutls can be found in the path ${PWD}/annotated/"
else
	echo $#
	echo "Incorrect number of arguments"
	echo "Please provide ANNOVAR input file name and family pedigree file"
fi
