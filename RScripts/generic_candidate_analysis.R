args = commandArgs(trailingOnly=T)
dir = args[3]
## IMPORTANT FAMILY INFO
fam.id<-strsplit(args[1],".avinput",fixed=T)
fam<-read.table(file=args[2],header=F,stringsAsFactors=F)
n.fam<-nrow(fam)
uf = fam$V2[which(fam$V6==1)]
n.uf = length(uf)
af<-fam$V2[which(fam$V6==2)]
n.af<-length(af)

id<-c(af,uf)
coverage.info = read.table(file=paste(dir,"/psap/lookups/gene_coverage_stats_final_12172014.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
low.coverage = coverage.info[which(coverage.info$Mean.Coverage < 10),"Gene"]
hgmd = read.table(file=paste(dir,"/psap/lookups/hgmd_pro_2013_4.12202014.annotated.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
hgmd.ad = unique(subset(hgmd,ModeInher == "AD")$Gene.wgEncodeGencodeBasicV19)
hgmd.ar = unique(subset(hgmd,ModeInher == "AR")$Gene.wgEncodeGencodeBasicV19)

models = c("DOM-het","REC-hom","REC-chet")
genos = c("het","hom")

# READ IN DATA
dat<-read.table(file=paste("annotated/",fam.id,"_popStat.txt",sep=""),sep="\t",header=T,stringsAsFactors=F,check.names=F)
dat$vid = paste(dat$Chr, dat$Start, dat$Alt,sep=":")

# ANALYSIS STEP 1: IDENTIFY SHARED VARIANTS AMONG AFFECTEDS
af.dat = list()

for(i in af){
  if(length(af.dat) == 0){
    for(m in 1:length(models)){
      tmp = dat[which(dat[,paste("Dz.Model.",i,sep="")] == models[m]),]
      if(models[m] == "REC-chet"){
        a1 = dat[which(dat[,paste("Dz.Model.",i,sep="")] == "DOM-het" & dat$Gene.wgEncodeGencodeBasicV19 %in% tmp$Gene.wgEncodeGencodeBasicV19),]
        a1[,paste("Dz.Model.",i,sep="")] = "REC-chet"
        tmp = rbind(tmp,dat[which(dat[,paste("Dz.Model.",i,sep="")] == "DOM-het"),])
      }
      af.dat[[m]] = tmp
    }
  }else{
    for(m in 1:length(models)){
      tmp = dat[which(dat[,paste("Dz.Model.",i,sep="")] == models[m]),]
      if(models[m] == "REC-het"){
        a1 = dat[which(dat[,paste("Dz.Model.",i,sep="")] == "DOM-het" & dat$Gene.wgEncodeGencodeBasicV19 %in% tmp$Gene.wgEncodeGencodeBasicV19),]
        a1[,paste("Dz.Model.",i,sep="")] = "REC-chet"
        tmp = rbind(tmp,dat[which(dat[,paste("Dz.Model.",i,sep="")] == "DOM-het"),])
      }
      af.dat[[m]] = af.dat[[m]][which(af.dat[[m]]$vid %in% tmp$vid),]
    }
  } 
}

# ANALYSIS STEP 2: VALIDATE VARIANT - if there are UFs
if(n.uf > 0){
  uf.dat = list()
  uf.ch = c()
  for(i in uf){
    if(length(uf.dat) == 0){
      for(g in 1:length(genos)){
        tmp = dat[which(dat[,i] == genos[g]),]
        uf.dat[[g]] = tmp
        if(genos[g] == "het"){
          gene.tab = table(af.dat[[3]]$Gene.wgEncodeGencodeBasicV19[which(af.dat[[3]][paste("Dz.Model.",af[1],sep="")] == "REC-chet" & af.dat[[3]]$vid %in% tmp$vid)])
          uf.ch = c(uf.ch,gene.tab[which(gene.tab > 1)])
        }
      }
    }else{
      for(g in 1:length(genos)){
        tmp = dat[which(dat[,i] == genos[g]),]
        uf.dat[[g]] = tmp
        if(genos[g] == "het"){
          gene.tab = table(af.dat[[3]]$Gene.wgEncodeGencodeBasicV19[which(af.dat[[3]][paste("Dz.Model.",af[1],sep="")] == "REC-chet" & af.dat[[3]]$vid %in% tmp$vid)])
          uf.ch = c(uf.ch,gene.tab[which(gene.tab > 1)])
        }
      }
    }
  }
  candidates = list()
  candidates[[1]] = af.dat[[1]][which(!af.dat[[1]]$vid %in% uf.dat[[1]]$vid & !af.dat[[1]]$vid %in% uf.dat[[2]]$vid),]
  candidates[[2]] = af.dat[[2]][which(!af.dat[[2]]$vid %in% uf.dat[[2]]$vid & !af.dat[[2]]$Gene.wgEncodeGencodeBasicV19 %in% uf.ch),]
  candidates[[3]] = af.dat[[3]][which(!af.dat[[3]]$Gene.wgEncodeGencodeBasicV19 %in% uf.ch & af.dat[[3]]$vid %in% uf.dat[[2]]$vid),]

  validated = data.frame()
  for(i in 1:length(candidates)) validated = rbind(validated,candidates[[i]])
}else{
  validated = data.frame()
  for(i in 1:length(af.dat)) validated = rbind(validated,af.dat[[i]])
}
        
## ANALYSIS STEP 3: FLAG SITES THAT HAVE LOW COVERAGE IN THE MAC61K DATA, MISSING ALLELE FREQUENCIES IN THE MAC61K DATA BUT NOT IN THE ESP AND 1000 GENOMES DATA, OR IS AN HGMD GENE WITH A MATCHING MODE OF INHERITANCE
validated$Flag = 0
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% low.coverage & is.na(validated$mac63kFreq_ALL)==F)] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% low.coverage & is.na(validated$mac63kFreq_ALL)==F)] + 100
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ad & validated[,paste("Dz.Model.",af[1],sep="")] == "DOM-het")] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ad & validated[,paste("Dz.Model.",af[1],sep="")] == "DOM-het")] + 20
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ar & validated[,paste("Dz.Model.",af[1],sep="")] %in% c("REC-hom","REC-chet"))] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ar & validated[,paste("Dz.Model.",af[1],sep="")] %in% c("REC-hom","REC-chet"))] + 3

## WRITE VALIDATED CANDIDATES TO FILE
write.table(validated[order(validated[,paste("popScore.",af[1],sep="")]),],file=paste("annotated/",fam.id,".report.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
