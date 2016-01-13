# Initialization script that downloads all necessary annotation tables prior to running the PSAP pipeline for the first time
ANNOVAR_PATH= #ADD PATH TO ANNOVAR DIRECTORY HERE
PSAP_PATH= #ADD PATH TO PSAP DIRECTORY HERE

# Downloads Sep 2014 1000 Genomes allele frequencies from ANNOVAR
perl ${ANNOVAR_PATH}/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar 1000g2014sep ${ANNOVAR_PATH}/humandb/

# Downloads CADD V1.3 scores from University of Washington and formats the data into annovar format
wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz 
# Formats CADD data into annovar format
zcat whole_genome_SNVs.tsv.gz | awk -F"\t" '!/#/ {print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5","$6}' > hg19_cadd.txt
wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi > ${1}/humandb/hg19_cadd.txt.idx
rm whole_genome_SNVs.tsv.gz

# Downloads ESP 6500 allele frequencies from ANNOVAR
perl ${ANNOVAR_PATH}/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar esp6500si_all ${ANNOVAR_PATH}/humandb/

# Downloads dbSNP137 annotations from ANNOVAR
perl ${ANNOVAR_PATH}/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar snp137 ${ANNOVAR_PATH}/humandb/

# Downloads GencodeV19 gene map from UCSC
perl ${ANNOVAR_PATH}/annotate_variation.pl -downdb -buildver hg19 wgEncodeGencodeBasicV19 ${ANNOVAR_PATH}/humandb/

# Moves ExAC allele frequencies provided with the PSAP pipeline to the ANNOVAR annotation folder
mv ${PSAP_PATH}/lookups/hg19_mac63k_Freq_ALL.txt.gz ${ANNOVAR_PATH}/humandb/
gzip -d ${ANNOVAR_PATH}/humandb/hg19_mac63k_Freq_ALL.txt.gz
