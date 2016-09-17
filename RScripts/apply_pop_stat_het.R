#EXTRACTING MAX pDa SCORES FOR MATERNAL AND PATERNAL HAPLOTYPE OF EACH GENE-- specifically, extracts the haplotype (genotype or diplotype)  with the max pDa score by gene based on the expected mode of inheritance.  one object for each mode of inheritance, will can be called in the spike in loop based on the genotpe from the HGMD variant included in the simulation
exome.AD<-subset(tmp,Geno=="het")

print("calculating summary scores")

if(nrow(exome.AD) > 0){
	popDat<-do.call(rbind,by(exome.AD,exome.AD$Gene.wgEncodeGencodeBasicV19,function(x) x[which.max(as.numeric(x[,score])),]))
	popDat$Dz.Model = "DOM-het"
## for the het case: pull mac26k het look up table, make the "exome" of tom het hit in all geans in the healthy exome and the HGMD variant calculate het genotype probs and fill in missing data with singleton rate
## use the mac26k look up table that estimates probabilites based on the expected number of hets
	lookup<-read.table(file=paste(dir,"/lookups/full.het.pCADD.gencodeV19.allsites.txt.gz",sep=""),stringsAsFactors=F,header=F)

print("beginning popStat claculations")
## BEGIN popStat Calculations
### to get column info from lookup table
### NOTE THAT THIS SUPPOSES A PREDEFINED TABLE FORMAT, GRID IN UNITS OF .001.
### FIRST TWO COLUMNS ARE GENE NAME AND NUMBER OF VARIANTS IN THAT GENE, RESPECTIVELY
	popDat$j<-findInterval(popDat[,score],scale)+1

### to get row info from lookup table
	popDat$i<-as.integer(factor(popDat$Gene.wgEncodeGencodeBasicV19,levels=lookup[,1]))

### pop.score p(C > c)
	popDat$popScore<-apply(popDat,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup)
	out = rbind(out,popDat)
}

