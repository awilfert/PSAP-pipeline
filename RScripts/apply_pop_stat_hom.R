#EXTRACTING MAX pDa SCORES FOR MATERNAL AND PATERNAL HAPLOTYPE OF EACH GENE-- specifically, extracts the haplotype (genotype or diplotype)  with the max pDa score by gene based on the expected mode of inheritance.  one object for each mode of inheritance, will can be called in the spike in loop based on the genotpe from the HGMD variant included in the simulation
exome.AR<-subset(tmp,Geno=="hom")

if(nrow(exome.AR) > 0){
	popDat<-do.call(rbind,by(exome.AR,exome.AR$Gene.wgEncodeGencodeBasicV19,function(x) x[which.max(x[,score]),]))
	popDat$Dz.Model = "REC-hom"
## for the hom case: mac26k hom look up table, make the "exome" of top hom hit in all genes in the healthy exome and the HGMD variant, calcualte hom genotype probs and fill in missing data with singleton rate
## use the mac26k look up table that estimated probabilities based on the expected number of homs
	lookup<-read.table(file=paste(dir,"/lookups/full.hom.pCADD.gencodeV19.allsites.txt.gz",sep=""),stringsAsFactors=F,header=F)

## BEGIN popStat Calculations
### to get column info from lookup table
### NOTE THAT THIS SUPPOSES A PREDEFINED TABLE FORMAT, GRID IN UNITS OF .001.
### FIRST TWO COLUMNS ARE GENE NAME AND NUMBER OF VARIANTS IN THAT GENE, RESPECTIVELY
	popDat$j<-findInterval(popDat[,score],scale)+1

## to get row info from lookup table
	popDat$i<-as.integer(factor(popDat$Gene.wgEncodeGencodeBasicV19,levels=lookup[,1]))

	popDat$popScore<-unlist(apply(popDat,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
	out = rbind(out,popDat)
}
