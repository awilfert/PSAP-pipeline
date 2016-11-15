#!/bin/Rscript

print(getwd())
## NAME OF FILE FOR ANALYSIS, PROVIDED AS AN ARGUMENT WHEN CALLING THE RSCRIPT
args<-commandArgs(trailingOnly=T) ## args[6] = family - must be annovar annotated (avinput.hg19_multianno.txt) and have a separate header file for the vcf colums(header); args[7] = family pedigree file

dir = args[3]

fam.id<-unlist(strsplit(args[1],".avinput",fixed=T))[1]

## PEDIGREE INFO
fam<-read.table(file=args[2],header=F,stringsAsFactors=F)
n.fam<-nrow(fam)
parents = c(unique(fam$V3[which(fam$V6==2 & fam$V3 !=0 )]),unique(fam$V4[which(fam$V6==2 & fam$V4 !=0)])) # ORDERED DAD THEN MOM
children = fam$V2[which(fam$V3 == parents[1] & fam$V4 == parents[2])]

## at some point this may be cahanged to an argument based system but for now it's hard coded
score = "scaled.cadd"
scale = seq(0,70,0.05)
lookup.genes = scan(file=paste(dir,"/psap/lookups/lookup_genes.txt",sep=""),"character")

## READ IN AND FORMAT DATA
exome.raw<-read.table(file=paste("annotated/",fam.id,".avinput.hg19_multianno.txt",sep=""),sep="\t",stringsAsFactors=F,skip=1)
header<-read.table(file=paste("annotated/",fam.id,".avinput.hg19_multianno.txt",sep=""),sep="\t",stringsAsFactors=F,nrow=1)
vcf.header<-read.table(file=paste(fam.id,".avinput.header",sep=""),sep="\t",stringsAsFactors=F,comment.char="@")
n.annos=ncol(header)
names(exome.raw)=c(header[-n.annos],vcf.header)
exome.raw$scaled.cadd = as.numeric(sapply(exome.raw$cadd,function(x) unlist(strsplit(x,","))[2]))

for(i in fam$V2){
  a1 = substr(exome.raw[,i],1,1)
  a2 = substr(exome.raw[,i],3,3)
  exome.raw[i] = NA
  if(length(which(a1 != a2)) > 0){
    exome.raw[which(a1 != a2),i] = "het"
  }
  if(length(which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2)) > 0){
    exome.raw[which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2),i] = "hom"
  }
  if(length(which(a1 == 0 & a2 == 0)) > 0){
    exome.raw[which(a1 == 0 & a2 == 0),i] = "ref"
  }
}

for(m in fam$V2[which(fam$V5 == 1)]){
  exome.raw[which(exome.raw[,m] %in% c("het","hom") & exome.raw[,"Chr"] == "X"),m] = "hom"
  exome.raw[which(exome.raw[,m] %in% c("het","hom") & exome.raw[,"Chr"] == "Y"),m] = "hom"
}

for(m in fam$V2[which(fam$V5 == 2)]){
  exome.raw[which(exome.raw[,"Chr"] == "Y"),m] = "ref"
}

### CLEAN DATA: 1) REMOVE BLACKLIST GENES, 2) REMOVE VARIANTS WITH AF DISCREPANCIES, 3) REMOVE GENES NOT INCLUDED IN LOOKUP TABLES, 4) REMOVE REGIONS THAT ARE NOT COVERED IN ExAC, 5) REMOVE MENDELIAN INCONSISTENCIES, 6) REMOVE VARIANTS THAT DO NOT PASS TRANCHE FILTER 
# 1) REMOVE BLACKLISTED GENES-- I think it would be better to remove and output these genes to a separate file (like the missing data).  I also think this should include all low coverage genes because it's more generic.
bl<-scan(paste(dir,"/psap/lookups/blacklist_122814.txt",sep=""),what="character")
bl.remove = unique(c(which(exome.raw$Gene.wgEncodeGencodeBasicV19 %in% bl),grep("^HLA", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^MUC", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^KRT", exome.raw$Gene.wgEncodeGencodeBasicV19),grep("^OR", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^TRBV", exome.raw$Gene.wgEncodeGencodeBasicV19)))

# 2) REMOVE AF DISCREPANCIES (ANYTHING THAT IS MISSING IN ExAC BUT PRESENT IN 1000GP OR ESP AT GREATER THAN 5% FREQUENCY
af.remove = which(is.na(exome.raw$mac63kFreq_ALL) == T & exome.raw[,"1000g2014sep_all"] > 0.05 | is.na(exome.raw$mac63kFreq_ALL) == T & exome.raw$esp6500si_all > 0.05)

# 3) REMOVE GENES NOT IN LOOKUP TABLES
lookup.remove = which(! exome.raw$Gene.wgEncodeGencodeBasicV19 %in% lookup.genes)

# 4) REMOVE LINES WHERE ALL AFs ARE MISSING
af = fam$V2[which(fam$V6 == 2)]
missing.remove = which(apply(exome.raw[af],1,function(row) return(sum(is.na(row)))) == length(af))

exome = exome.raw[-unique(c(bl.remove,af.remove,lookup.remove,missing.remove)),]

# 4) REMOVE MENDELIAN INCONSITENCIES - IF AFs HAVE BOTH PARENTS
autosome = data.frame(parents=c("ref,ref","ref,het","ref,hom","het,ref","het,het","het,hom","hom,ref","hom,het","hom,hom"),child=c("ref,het","ref,het","het","ref,het","ref,het,hom","het,hom","het","het,hom","hom"),stringsAsFactors=F)
male.x = data.frame(mom=c("ref","het","hom"),child=c("ref","ref,hom","hom"),stringsAsFactors=F)
female.x = data.frame(parents=c("ref,ref","ref,het","ref,hom","hom,ref","hom,het","hom,hom"),child=c("ref","ref,het","het","het","het,hom","hom"),stringsAsFactors=F)
male.y = data.frame(dad=c("ref","hom"),child=c("ref","hom"),stringsAsFactors=F)

inconsistencies = function(row){
  parent.genos = paste(row[parents[1]],row[parents[2]],sep=",")
  count = 0
  for(i in children){
    if(sum(is.na(row[fam$V2])) > 0){
      count = count + 1
    }else{
      if(row["Chr"] == "X" & fam$V5[which(fam$V2 == i)] == 1){ # MALE X
        genos = male.x[which(male.x$mom == row[parents[2]]),2] # ONLY LOOK AT MOM'S X
        if(nchar(genos) > 3){
          ok.genos = unlist(strsplit(genos,","))
        }else{
          ok.genos = genos
        }
      }else{
        if(row["Chr"] == "Y" & fam$V5[which(fam$V2 == i)] == 1){ # MALE Y
          genos = male.y[which(male.y$dad == row[parents[1]]),2] # ONLY LOOK AT DAD'S Y
          if(nchar(genos) > 3){
            ok.genos = unlist(strsplit(genos,","))
          }else{
            ok.genos = genos
          }
        }else{
          if(row["Chr"] == "X" & fam$V5[which(fam$V2 == i)] == 2){ # FEMALE X
            genos = female.x[which(female.x$parents == parent.genos),2] # LOOK AT BOTH PARENTS
            if(nchar(genos) > 3){
              ok.genos = unlist(strsplit(genos,","))
            }else{
              ok.genos = genos
            }
          }else{ # AUTOSOMES
            genos = autosome[which(autosome$parents == parent.genos),2] # LOOK AT BOTH PARENTS
            if(nchar(genos) > 3){
              ok.genos = unlist(strsplit(genos,","))
            }else{
              ok.genos = genos
            }
          }
        }
      }
      if(row[i] %in% ok.genos){
        count = count + 1
      }
    }
  }
  if(count == length(children)){
    return(matrix(row,nrow=1,byrow=T))
  }
}

if(length(which(parents != 0)) == 2){
  print("1 true")
  tmp.exome = data.frame(matrix(unlist(apply(exome,1,inconsistencies)),ncol=length(exome),byrow=T),stringsAsFactors=F)
  names(tmp.exome) = names(exome)
}else{
  print("2 true")
  tmp.exome = exome
}

# 5a) REMOVE REGIONS NOT COVERED IN ExAC - JUST RETAINING PROTEIN CODING SITES AND SPLICE SITES
keep<-unique(c(grep("splic",tmp.exome$Func.wgEncodeGencodeBasicV19),which(is.na(tmp.exome$ExonicFunc.wgEncodeGencodeBasicV19)==FALSE)))
exome<-tmp.exome[keep,]

# 5b) SCORE INDELS
lookup.lof = read.table(file=paste(dir,"/psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz",sep=""),stringsAsFactors=F)
indels = grep("^frameshift",exome$ExonicFunc.wgEncodeGencodeBasicV19)
gene.index = as.integer(factor(exome$Gene.wgEncodeGencodeBasicV19[indels],levels=lookup.lof[,1]))
exome$scaled.cadd[indels] = lookup.lof[gene.index,2]

# 6) REMOVE VARIANTS THAT DO NOT PASS QUALITY FILTER OR HAVE MISSING pCADD SCORES
info<-exome[which(exome$FILTER=="PASS" & is.na(exome[,score]) == F | exome$FILTER=="." & is.na(exome[,score]) == F),c("Chr","Start","Ref","Gene.wgEncodeGencodeBasicV19","Func.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19","mac63kFreq_ALL","1000g2014sep_all","esp6500si_all","Alt","cadd","scaled.cadd",fam$V2)]

## OUTPUT MISSING DATA/DATA NOT INCLUDED IN ANY OF THE ABOVE ANALYSES
id.raw = paste(exome.raw$Chr,exome.raw$Start,exome.raw$Ref,exome.raw$Alt,sep=":")
id.final = paste(info$Chr,as.numeric(info$Start),info$Ref,info$Alt,sep=":")
missing<-unique(exome.raw[which(! id.raw %in% id.final),c("Chr","Start","Ref","Gene.wgEncodeGencodeBasicV19","Func.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19","mac63kFreq_ALL","1000g2014sep_all","esp6500si_all","Alt","cadd","scaled.cadd",fam$V2)])

rm(list=c("keep","exome","tmp.exome","exome.raw","af.remove","lookup.remove","bl.remove","bl","lookup.genes"))
class(info[,score]) = "numeric"

indv.cols =c()
for(i in fam$V2) indv.cols = c(indv.cols,which(names(info) == i))

print("data cleaned, beginning loop")
final = info
for(m in indv.cols){
  print(m)
  tmp = info[-indv.cols]
  tmp["Geno"] = info[m]
  out = data.frame()

  print("beginning annotations")
  ## SOURCES CODE THAT WILL FORMAT AND ANALYZE THE DATA FOR EACH MODE OF INHERITANCE (AD, AR, CHET or X) MODEL
  ## AD MODEL CODE
  source(paste(dir,"/psap/RScripts/apply_pop_stat_het.R",sep=""))
  print("AD model complete")
  
  ## AR MODEL CODE
  source(paste(dir,"/psap/RScripts/apply_pop_stat_hom.R",sep=""))
  print("AR-hom model complete")
  
  ## CHET MODEL CODE
  source(paste(dir,"/psap/RScripts/apply_pop_stat_chet_unphased.R",sep=""))
  print("AR-het model complete")
  
  print("processing data")
  names(out)[which(names(out) == "popScore")] = paste("popScore",names(info)[m],sep=".")
  names(out)[which(names(out) == "Dz.Model")] = paste("Dz.Model",names(info)[m],sep=".")
  out[which(out[,paste("popScore.",names(info)[m],sep="")] == 0),paste("popScore.",names(info)[m],sep="")] = 1e-6
  out = out[which(!names(out) %in% c("i","j"))]
  names(out)[which(names(out) == "Geno")] = names(info)[m]
  final = merge(final,out,all.x=T,all.y=T,by=c("Chr","Start","Ref","Gene.wgEncodeGencodeBasicV19","Func.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19","mac63kFreq_ALL","1000g2014sep_all","esp6500si_all","Alt","cadd","scaled.cadd",names(info)[m]))
}

## WRITE OUTPUT FOR FAMILY
print("writing file")
write.out = final[-which(rowSums(is.na(final[grep("popScore.",names(final))])) == n.fam),]
write.table(write.out,file=paste("annotated/",fam.id,"_popStat.txt",sep=""),col.names=T,row.names=F,sep="\t")
missing.out = missing[-which(missing$Func.wgEncodeGencodeBasicV19 == "exonic" | missing$Func.wgEncodeGencodeBasicV19 == "splicing;intronic" | missing$Func.wgEncodeGencodeBasicV19 =="splicing;exonic" | missing$Func.wgEncodeGencodeBasicV19 =="exonic;splicing"),]
write.table(missing.out,file=paste("annotated/",fam.id,"_missing_data.txt",sep=""),sep="\t",col.names=T,row.names=F)

