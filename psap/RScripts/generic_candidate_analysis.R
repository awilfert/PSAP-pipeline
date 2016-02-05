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
	uf.dat = data.frame()
	uf.dat = dat[unlist(lapply(1:nrow(dat),function(i) if(any(dat[i,uf] != "ref")) return(i))),]

  candiates = data.frame()
  for(i in 1:length(af.dat)) candidates= rbind(candidate,af.dat[[i]])
  
  candidates$validation = "ok"
  af.het.rows = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model",af,sep="")] == "DOM-het") & candidates$vid %in% uf.dat$vid) return(i))) # identifies het candidates in affected(s) that are also present in at least one unaffected
  candidates$validation[af.het.rows] = "violation" # flags the above identifeid variants as a violation of the the disease inheritance model
  uf.rec.rows = unlist(lapply(1:nrow(uf.dat), function(i) if(any(uf.dat[i,uf] != "DOM-het")) return(i))) # identifies rows in the unaffected data that are not DOM-het (all other hets are allowed)
  af.rec.rows = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model",af,sep="")] != "DOM-het") & candidates$vid %in% uf.dat$vid[uf.rec.rows]) return(i))) # identifies REC disease model candidates that are also present in at least one unaffected individual as a REC-hom or REC-chet
  candidates$validation[af.rec.rows] = "violation"
		
	# missing mate
	# 1) pull out chets tagged with violation
	af.chet.rows = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model",af,sep="")] != "REC-chet") & candidates$validation == "violation") return(i))) # identifies REC-chet rows in affected individuals that were present in unaffecteds as a REC-hom or REC-chet
	violations = candidates[af.chet.rows,]
	# 2) loop over violated chet genes
	for(g in violations$gene){
	#	3) look at subset of non-violated chets in genes with a chet violation
		ok.chets = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model",af,sep="")] == "REC-chet") & candidates$validation == "ok" & candidates$Gene.wgEncodeGencodeBasicV19 == g) return(i)))
	# 4) remove chets where mate vioated disease inheritance model
		candidates$violations = "mate violation"
	}
}else{
  candidates = data.frame()
  for(i in 1:length(af.dat)) candidates = rbind(candidates,af.dat[[i]])
  candidates$validation = "ok"
}

validated = candidates[which(candidates$validation == "ok"),]
    
## ANALYSIS STEP 3: FLAG SITES THAT HAVE LOW COVERAGE IN THE MAC61K DATA, MISSING ALLELE FREQUENCIES IN THE MAC61K DATA BUT NOT IN THE ESP AND 1000 GENOMES DATA, OR IS AN HGMD GENE WITH A MATCHING MODE OF INHERITANCE
validated$Flag = 0
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% low.coverage & is.na(validated$mac63kFreq_ALL)==F)] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% low.coverage & is.na(validated$mac63kFreq_ALL)==F)] + 100
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ad & validated[,paste("Dz.Model.",af[1],sep="")] == "DOM-het")] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ad & validated[,paste("Dz.Model.",af[1],sep="")] == "DOM-het")] + 20
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ar & validated[,paste("Dz.Model.",af[1],sep="")] %in% c("REC-hom","REC-chet"))] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ar & validated[,paste("Dz.Model.",af[1],sep="")] %in% c("REC-hom","REC-chet"))] + 3

## WRITE VALIDATED CANDIDATES TO FILE
write.table(validated[order(validated[,paste("popScore.",af[1],sep="")]),],file=paste("annotated/",fam.id,".report.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
