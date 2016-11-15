#EXTRACTING MAX pDa SCORES FOR MATERNAL AND PATERNAL HAPLOTYPE OF EACH GENE-- specifically, extracts the haplotype (genotype or diplotype)  with the max pDa score by gene based on the expected mode of inheritance.  one object for each mode of inheritance, will can be called in the spike in loop based on the genotpe from the HGMD variant included in the simulation
exome.AD<-subset(tmp,Geno=="het")
id<-paste(exome.AD$chr,exome.AD$Start)
exome.AD<-exome.AD[!duplicated(id),]

#just select genes with more than one variant
dups<-names(table(exome.AD$Gene.wgEncodeGencodeBasicV19)[which(table(exome.AD$Gene.wgEncodeGencodeBasicV19) >= 2)])
exome.AD<-exome.AD[exome.AD$Gene.wgEncodeGencodeBasicV19 %in% dups,]

if(nrow(exome.AD) > 0){
	popDat<-do.call(rbind,by(exome.AD,exome.AD$Gene.wgEncodeGencodeBasicV19, function(x) { x$scaled.cadd<-as.numeric(x$scaled.cadd); x<-x[order(x$scaled.cadd,decreasing=T,na.last=NA),]; return(x[2,]) }))
	popDat$Dz.Model = "REC-chet"
## take the min of the two variants (will require a "by" statement)
#cols = paste(score,c(".x",".y"),sep="")
#popDat$min.score<-as.numeric(apply(popDat[cols],1,min))
#popDat<-popDat[which(popDat$Start.x != popDat$Start.y & popDat$min.score>0),]

	lookup<-read.table(file=paste(dir,"/psap/lookups/full.chet.pCADD.gencodeV19.allsites.txt.gz",sep=""),stringsAsFactors=F,header=F)

### pop.score p(V|G)
## BEGIN popStat Calculations
### to get column info from lookup table
### NOTE THAT THIS SUPPOSES A PREDEFINED TABLE FORMAT, GRID IN UNITS OF .001.
### FIRST TWO COLUMNS ARE GENE NAME AND NUMBER OF VARIANTS IN THAT GENE, RESPECTIVELY
	popDat$j<-findInterval(popDat[,score],scale)+1

### to get row info from lookup table
	popDat$i<-as.integer(factor(popDat$Gene.wgEncodeGencodeBasicV19,levels=lookup[,1]))

### pop.score
	popDat$popScore<-unlist(apply(popDat,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
	out = rbind(out,popDat)
}
