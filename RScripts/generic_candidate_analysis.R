args = commandArgs(trailingOnly=T)
dir = args[3]
## IMPORTANT FAMILY INFO
fam.id<-unlist(strsplit(args[1],".avinput",fixed=T))[1]
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
print("databases loaded")
# READ IN DATA
dat<-read.table(file=paste("annotated/",fam.id,"_popStat.txt",sep=""),sep="\t",header=T,stringsAsFactors=F,check.names=F)
dat$vid = paste(dat$Chr, dat$Start, dat$Alt,sep=":")

# ANALYSIS STEP 1: IDENTIFY SHARED VARIANTS AMONG AFFECTEDS
af.dat = list()

for(i in af){
    if(length(af.dat) == 0){
        for(m in 1:length(models)){
            tmp = dat[which(dat[,paste("Dz.Model.",i,sep="")] == models[m]),] # extract all CHETs from af
            if(models[m] == "REC-chet" & nrow(tmp) > 0){
                a1 = dat[which(dat[,paste("Dz.Model.",i,sep="")] == "DOM-het" & dat$Gene.wgEncodeGencodeBasicV19 %in% tmp$Gene.wgEncodeGencodeBasicV19),] # extract all CHET mates
                a1 = merge(a1[-which(names(a1) %in% c(paste("popScore.",i,sep=""),paste("Dz.Model.",i,sep="")))],tmp[c("Gene.wgEncodeGencodeBasicV19",paste("popScore.",i,sep=""),paste("Dz.Model.",i,sep=""))]) # rescore the DOM-het PSAPs to be equal to the REC-chet mates' PSAP and rename the DZ model for the mates
                tmp = rbind(tmp,a1)
            }
            af.dat[[m]] = tmp
        }
    }else{
        for(m in 1:length(models)){
            tmp = dat[which(dat[,paste("Dz.Model.",i,sep="")] == models[m]),]
            if(models[m] == "REC-het" & nrow(tmp) > 0 ){
                a1 = dat[which(dat[,paste("Dz.Model.",i,sep="")] == "DOM-het" & dat$Gene.wgEncodeGencodeBasicV19 %in% tmp$Gene.wgEncodeGencodeBasicV19),] # extract all CHET mates
                a1[,paste("Dz.Model.",i,sep="")] = "REC-chet" # rename the DZ model for the mates
                a1 = merge(a1[-which(names(a1) == paste("popScore.",i,sep=""))],tmp[c("Gene.wgEncodeGencodeBasicV19",paste("popScore.",i,sep=""))]) # rescore the DOM-het PSAPs to be equal to the REC-chet mates' PSAP
                tmp = rbind(tmp,a1)
            }
            af.dat[[m]] = af.dat[[m]][which(af.dat[[m]]$vid %in% tmp$vid),]
        }
    } 
}

# ANALYSIS STEP 2: VALIDATE VARIANT - if there are UFs
if(n.uf > 0){
    af.cols = c()
    for(i in af) af.cols = c(af.cols,grep(i,names(dat)))
    uf.dat = data.frame()
    uf.dat = dat[,-af.cols]

    candidates = data.frame()
    for(i in 1:length(af.dat)) candidates= rbind(candidates,af.dat[[i]])
    
    print("checking for disease model violations")
    candidates$validation = "ok"
    
    print("DOM-het disease model done")
    uf.het.rows = unlist(lapply(1:nrow(dat),function(i) if(any(!dat[i,uf] %in% c("ref",NA))) return(i))) # identify all rows where at least one UF is non-missing/ref
    af.het.rows = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model.",af,sep="")] == "DOM-het") & candidates$vid[i] %in% uf.dat$vid[uf.het.rows]) return(i))) # identifies het candidates in affected(s) that are also present in at least one unaffected
    candidates$validation[af.het.rows] = "violation" # flags the above identified variants as a violation of the the disease inheritance model
    
    print("REC-hom disease model done")
    uf.rec.rows = c()
    for(id in uf) uf.rec.rows = c(uf.rec.rows,grep("DOM-het",uf.dat[,paste("Dz.Model.",id,sep="")],invert=T)) # identifies rows in the unaffected data that are not DOM-het (all other genotypes and disease models are considered)
    uf.rec.rows.tab = table(uf.rec.rows)
    uf.rec.rows = as.numeric(names(uf.rec.rows.tab[which(uf.rec.rows.tab == length(uf))]))
    
    af.hom.rows = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model.",af,sep="")] == "REC-hom") & candidates$vid[i] %in% uf.dat$vid[uf.rec.rows]) return(i))) # identifies REC disease model candidates that are also present in at least one unaffected individual as a REC-hom or REC-chet
    candidates$validation[af.hom.rows] = "violation"
	
	print("REC-chet disease model done")
	af.chet.rows = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model.",af,sep="")] == "REC-chet") & candidates$vid[i] %in% uf.dat$vid[uf.rec.rows]) return(i))) # identifies REC disease model candidates that are also present in at least one unaffected individual as a REC-hom or REC-chet
    candidates$validation[af.chet.rows] = "violation"
	
	print("checking for chet mate violations")	
	# missing mate
	# 1) pull out chets tagged with violation
	af.chet.ok.rows = unlist(lapply(1:nrow(candidates), function(i) if(any(candidates[i,paste("Dz.Model.",af,sep="")] == "REC-chet") & candidates$validation[i] == "violation") return(i))) # identifies REC-chet rows in affected individuals that were present in unaffecteds as a REC-hom or REC-chet
	violations = candidates[af.chet.ok.rows,"Gene.wgEncodeGencodeBasicV19"]
	# 2) loop over violated chet genes
	missing_mate.chets = which(candidates[,paste("Dz.Model.",af,sep="")] == "REC-chet" & candidates$validation == "ok" & candidates$Gene.wgEncodeGencodeBasicV19 %in% violations)
	# 4) remove chets where mate vioated disease inheritance model
	candidates$validation[missing_mate.chets] = "mate violation"
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
