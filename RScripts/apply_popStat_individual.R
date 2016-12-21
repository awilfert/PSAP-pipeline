#!/bin/Rscript

print(getwd())
## NAME OF FILE FOR ANALYSIS, PROVIDED AS AN ARGUMENT WHEN CALLING THE RSCRIPT
args<-commandArgs(trailingOnly=T) ## args[1] = family - must be annovar annotated (avinput.hg19_multianno.txt) and have a separate header file for the vcf colums(header); args[2] = individual ID
dir = args[3]

fam.id<-strsplit(args[1],".avinput",fixed=T)
ped = read.table(args[4],stringsAsFactors=F)

## Individual ID - ASSUMES only one individual is being analyzed/annotated
indv.id = args[2]

## at some point this may be changed to an argument based system but for now it's hard coded
score = "CADD_Phred"
scale = seq(0,70,0.05)
lookup.genes = scan(file=paste(dir,"psap/lookups/lookup_genes.txt",sep=""),"character")

## READ IN AND FORMAT DATA
exome.raw<-read.table(file=paste("annotated/",fam.id,".avinput.hg19_multianno.txt",sep=""),sep="\t",stringsAsFactors=F,skip=1)
header<-read.table(file=paste("annotated/",fam.id,".avinput.hg19_multianno.txt",sep=""),sep="\t",stringsAsFactors=F,nrow=1)
vcf.header<-read.table(file=paste(fam.id,".avinput.header",sep=""),sep="\t",stringsAsFactors=F,comment.char="@")
stopifnot(indv.id %in% vcf.header) # CHECKS THAT SPECIFIED INDIVIDUAL IS IN THE DATA
n.annos=ncol(header)
names(exome.raw)=c(header[-n.annos],vcf.header)

if(any(grepl("cadd",names(exome.raw))) == T){
    exome.raw[,score] = as.numeric(sapply(exome.raw$cadd,function(x) unlist(strsplit(x,","))[2]))
}

stopifnot(any(grepl("CADD_Phred",names(exome.raw))))

# Extracts genotype info for specified individual
a1 = substr(exome.raw[,indv.id],1,1)
a2 = substr(exome.raw[,indv.id],3,3)
exome.raw[indv.id] = "NA"
if(length(which(a1 != a2)) > 0){
  exome.raw[which(a1 != a2),indv.id] = "het"
}
if(length(which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2)) > 0){
  exome.raw[which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2),indv.id] = "hom"
}
if(length(which(a1 == 0 & a2 == 0)) > 0){
  exome.raw[which(a1 == 0 & a2 == 0),indv.id] = "ref"
}

if(ped$V5[which(ped$V2 == indv.id)] == 1){
  exome.raw[which(exome.raw[,indv.id] %in% c("het","hom") & exome.raw[,"Chr"] == "X"),indv.id] = "hom"
  exome.raw[which(exome.raw[,indv.id] %in% c("het","hom") & exome.raw[,"Chr"] == "Y"),indv.id] = "hom"
}

if(ped$V5[which(ped$V2 == indv.id)] == 2){
  exome.raw[which(exome.raw[,"Chr"] == "Y"),indv.id] = "ref"
}
print(dim(exome.raw))

### CLEAN DATA: 1) REMOVE BLACKLIST GENES, 2) REMOVE VARIANTS WITH AF DISCREPANCIES, 3) REMOVE GENES NOT INCLUDED IN LOOKUP TABLES, 4) REMOVE REGIONS THAT ARE NOT COVERED IN ExAC, 5) REMOVE VARIANTS THAT DO NOT PASS TRANCHE FILTER 
# 1) REMOVE BLACKLISTED GENES-- I think it would be better to remove and output these genes to a separate file (like the missing data).  I also think this should include all low coverage genes because it's more generic.
bl<-scan(paste(dir,"psap/lookups/blacklist_122814.txt",sep=""),what="character")
bl.remove = unique(c(which(exome.raw$Gene.wgEncodeGencodeBasicV19 %in% bl),grep("^HLA", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^MUC", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^KRT", exome.raw$Gene.wgEncodeGencodeBasicV19),grep("^OR", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^TRBV", exome.raw$Gene.wgEncodeGencodeBasicV19)))

# 2) REMOVE AF DISCREPANCIES (ANYTHING THAT IS MISSING IN ExAC BUT PRESENT IN 1000GP OR ESP AT GREATER THAN 5% FREQUENCY
af.remove = which(is.na(exome.raw$mac63kFreq_ALL) == T & exome.raw[,"1000g2014sep_all"] > 0.05 | is.na(exome.raw$mac63kFreq_ALL) == T & exome.raw$esp6500si_all > 0.05)

# 3) REMOVE GENES NOT IN LOOKUP TABLES
lookup.remove = which(! exome.raw$Gene.wgEncodeGencodeBasicV19 %in% lookup.genes)

# 4) REMOVE LINES WHERE ALL AFs ARE MISSING
af = ped$V2[which(ped$V6 == 2)]
missing.remove = which(apply(exome.raw[af],1,function(row) return(sum(is.na(row)))) == length(af))

tmp.exome = exome.raw[-unique(c(bl.remove,af.remove,lookup.remove,missing.remove)),]

# 4a) REMOVE REGIONS NOT COVERED IN ExAC - JUST RETAINING PROTEIN CODING SITES AND SPLICE SITES
keep<-unique(c(grep("splic",tmp.exome$Func.wgEncodeGencodeBasicV19),which(is.na(tmp.exome$ExonicFunc.wgEncodeGencodeBasicV19)==FALSE)))
exome<-tmp.exome[keep,]

# 4b) SCORE INDELS
lookup.lof = read.table(file=paste(dir,"psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz",sep=""),stringsAsFactors=F)
indels = grep("^frameshift",exome$ExonicFunc.wgEncodeGencodeBasicV19)
gene.index = as.integer(factor(exome$Gene.wgEncodeGencodeBasicV19[indels],levels=lookup.lof[,1]))
exome[,score][indels] = lookup.lof[gene.index,2]

# 5) REMOVE VARIANTS THAT DO NOT PASS QUALITY FILTER OR HAVE MISSING pCADD SCORES
info<-exome[which(exome$FILTER=="PASS" & is.na(exome[,score]) == F | exome$FILTER=="." & is.na(exome[,score]) == F),c("Chr","Start","Ref","Gene.wgEncodeGencodeBasicV19","Func.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19","mac63kFreq_ALL","1000g2014sep_all","esp6500si_all","Alt",score,indv.id)]

## OUTPUT MISSING DATA/DATA NOT INCLUDED IN ANY OF THE ABOVE ANALYSES
id.raw = paste(exome.raw$Chr,exome.raw$Start,exome.raw$Ref,exome.raw$Alt,sep=":")
id.final = paste(info$Chr,as.numeric(info$Start),info$Ref,info$Alt,sep=":")
missing<-unique(exome.raw[which(! id.raw %in% id.final),c("Chr","Start","Ref","Gene.wgEncodeGencodeBasicV19","Func.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19","mac63kFreq_ALL","1000g2014sep_all","esp6500si_all","Alt",score,indv.id)])

rm(list=c("keep","exome","tmp.exome","exome.raw","af.remove","lookup.remove","bl.remove","bl","lookup.genes"))
class(info[,score]) = "numeric"

print("data cleaned, beginning annotation")
tmp = info[-which(names(info) == indv.id)]
tmp["Geno"] = info[,which(names(info) == indv.id)]
print(dim(tmp))
out = data.frame()

## SOURCES CODE THAT WILL FORMAT AND ANALYZE THE DATA FOR EACH MODE OF INHERITANCE (AD, AR, CHET or X) MODEL
## AD MODEL CODE
source(paste(dir,"psap/RScripts/apply_pop_stat_het.R",sep=""))
print("AD model complete")
  
## AR MODEL CODE
source(paste(dir,"psap/RScripts/apply_pop_stat_hom.R",sep=""))
print("AR-hom model complete")
  
## CHET MODEL CODE
source(paste(dir,"psap/RScripts/apply_pop_stat_chet_unphased.R",sep=""))
print("AR-het model complete")
print(dim(out))  
print("processing data")
out = out[which(!names(out) %in% c("i","j"))]
out[which(out$popScore == 0),"popScore"] = 1e-6

## WRITE OUTPUT FOR FAMILY
print("writing file")
write.out = out[which(is.na(out$popScore) == F),]
print(dim(write.out))
write.table(write.out,file=paste("annotated/",fam.id,"_",indv.id,"_popStat.txt",sep=""),col.names=T,row.names=F,sep="\t")
missing.out = missing[-which(missing$Func.wgEncodeGencodeBasicV19 == "exonic" | missing$Func.wgEncodeGencodeBasicV19 == "splicing;intronic" | missing$Func.wgEncodeGencodeBasicV19 =="splicing;exonic" | missing$Func.wgEncodeGencodeBasicV19 =="exonic;splicing"),]
write.table(missing.out,file=paste("annotated/",fam.id,"_",indv.id,"_missing_data.txt",sep=""),sep="\t",col.names=T,row.names=F)

