# Initialization script that downloads all necessary annotation tables prior to running the PSAP pipeline for the first time
ANNOVAR_PATH=/lts/dclab/amy/annovar
PSAP_PATH=~/psap

cd $ANNOVAR_PATH
# Download Sep 2014 1000 Genomes allele frequencies from ANNOVAR
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar 1000g2014sep humandb/

# Download CADD V1.3 scores from University of Washington and formats the data into annovar format
wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz 
# Format CADD data into ANNOVAR format
zcat whole_genome_SNVs.tsv.gz | awk -F"\t" '!/#/ {print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5","$6}' > humandb/hg19_cadd.txt
wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi > humandb/hg19_cadd.txt.idx
rm whole_genome_SNVs.tsv.gz

# Download ESP 6500 allele frequencies from ANNOVAR
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar esp6500si_all humandb/

# Download dbSNP137 annotations from ANNOVAR
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar snp137 humandb/

# Download GencodeV19 gene map from UCSC
perl annotate_variation.pl -downdb -buildver hg19 wgEncodeGencodeBasicV19 humandb/

# Move ExAC allele frequencies provided with the PSAP pipeline to the ANNOVAR annotation folder and unzip
mv ${PSAP_PATH}/lookups/hg19_mac63kFreq_ALL.txt.gz humandb/
gzip -d humandb/hg19_mac63kFreq_ALL.txt.gz
