#!/usr/bin/env Rscript

library("optparse")
# get options

option_list = list(
  make_option(c("-v", "--VCF"), type="character", default=NULL, help="multisample VCF file from mutect2 [default= %default]", metavar="character"),
  make_option(c("-C", "--CNV_tumors"), type="character", default=NULL, help="CNV file from hatchet for tumors [default= %default]", metavar="character"),
  make_option(c("-n", "--normal_id"), type="character", default=NULL, help="ID of the normal sample in VCF [default= %default]", metavar="character"),
  make_option(c("-d", "--min_dp"), type="numeric", default=20, help="minimum dp to consider a mutation [default= %default]", metavar="character"),
  make_option(c("-c", "--min_cnv_length"), type="numeric", default=1, help="minimum length in bp to consider a copy number [default= %default]", metavar="numeric"),
  make_option(c("-s", "--chunk_size"), type="numeric", default=100000, help="size of VCF chunks for reading [default= %default]", metavar="numeric"),
  make_option(c("-o", "--output_folder"), type="character", default=".", help="folder to save the results [default= %default]", metavar="character"),
  make_option(c("-G", "--kept_genes"), type="character", default=NULL, help="Gene spans file to compute FPKMs [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

print(args)

library(VariantAnnotation)

if(is.null(args$output_folder)) {output_folder="."} else {output_folder=args$output_folder}
if(is.null(args$chunk_size)) {chunk_size = 100000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$min_dp)) {min_dp = 20} else {min_dp = as.numeric(args$min_dp)}
if(is.null(args$min_cnv_length)) {min_cnv_length = 1} else {min_cnv_length = as.numeric(args$min_cnv_length)}
if(is.null(args$kept_genes)) {kept_genes = c()} else {kept_genes = unlist(strsplit(args$kept_genes,","))}

normal_id = args$normal_id

print( args$CNV_tumors )

CNV_tumors = read.table(args$CNV_tumors ,comment.char = "", header = T)
colnames( CNV_tumors ) = gsub("X.","",colnames(CNV_tumors))

# only consider CNV with a sufficient length to de-noise the data
CNV_tumors$seglength = CNV_tumors$END - CNV_tumors$START + 1
CNV_tumors = CNV_tumors[which(CNV_tumors$seglength > min_cnv_length),]

if(!file.exists(paste(args$VCF,".tbi",sep=""))){
  system(paste("tabix -p vcf ", args$VCF, sep=""))
}

vcf <- open(VcfFile(args$VCF,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf, "hg38")

tumors_reformated4pyclone = vector("list",ncol(vcf_chunk)-1)
for(i in 1:(ncol(vcf_chunk)-1) ) tumors_reformated4pyclone[[i]] = data.frame(mutation_id=NULL, ref_counts=NULL, var_counts=NULL, normal_cn=NULL, minor_cn=NULL, major_cn=NULL, gene=NULL)
funcl = c()
exonicfuncl = c()
while(dim(vcf_chunk)[1] != 0) {
  #muts = unlist(lapply(rownames(vcf_chunk), function(s) unlist(strsplit(s,"_"))[[2]]))
  vcf_chunk = vcf_chunk#[which(nchar(muts)==3),]
  vcf_chunk = vcf_chunk[-as.numeric( grep("HLA|chrUn|_decoy|_alt|_random", rowRanges(vcf_chunk)@seqnames) ),] # remove mutations on weird chromosomes
  AD_matrix = geno(vcf_chunk,"AD")
  GT_matrix = geno(vcf_chunk,"GT")
  print(colnames(GT_matrix))
  DP_matrix = geno(vcf_chunk,"DP")
  gene = unlist(info(vcf_chunk)$Gene.refGene) # possibly change to Gene.ensGene if ensembl ref used by annovar
  tumors_mutid = as.numeric(which( (GT_matrix[,normal_id]=="0/0" | GT_matrix[,normal_id]=="0|0") &
                                                ( rowMeans(DP_matrix[,colnames(DP_matrix)!=normal_id]>min_dp) ==1
                                                  | gene %in% kept_genes)))

  mutid=rownames(AD_matrix)
  #print(head(GT_matrix))
  #print(head(DP_matrix))
  print(tumors_mutid)
  #print( colnames(DP_matrix)!=normal_id )
  #print(head(GT_matrix[,normal_id]))
  # for the tumor 1 sample
  #print( head( AD_matrix[tumors_mutid,colnames(AD_matrix)!=normal_id] ) )

  rc = sapply(1:(ncol(AD_matrix)-1), function(i) matrix(unlist( AD_matrix[tumors_mutid,colnames(AD_matrix)!=normal_id][,i] ),ncol=2,byrow=T)[,1] ) #unlist(lapply( tumors_mutid, function(i) AD_matrix[i,colnames(DP_matrix)!=normal_id] ))
  vc = sapply(1:(ncol(AD_matrix)-1), function(i) matrix(unlist( AD_matrix[tumors_mutid,colnames(AD_matrix)!=normal_id][,i] ),ncol=2,byrow=T)[,2] )
  gene = unlist(info(vcf_chunk)$Gene.refGene)[tumors_mutid] #possibly change to Gene.ensGene if ensembl ref used by annovar
  gene = unlist(lapply(gene, function(g) unlist(strsplit(g,"[\\]"))[1]))
  
  #funcl     = c(funcl, unlist(info(vcf_chunk)$Func.refGene)[tumors_mutid] )
  func       = unlist(info(vcf_chunk)$Func.refGene)[tumors_mutid] #possibly change to Func.ensGene if ensembl ref used by annovar
  #exonicfuncl = c(exonicfuncl, unlist(info(vcf_chunk)$ExonicFunc.refGene)[tumors_mutid] )
  exonicfunc = unlist(info(vcf_chunk)$ExonicFunc.refGene)[tumors_mutid] #possibly change to ExonicFunc.ensGene if ensembl ref used by annovar
  for(i in 1:ncol(rc)) tumors_reformated4pyclone[[i]] = rbind(tumors_reformated4pyclone[[i]],
           data.frame(mutation_id=mutid[tumors_mutid], ref_counts=rc[,i], var_counts=vc[,i], normal_cn=NA, minor_cn=NA, major_cn=NA, gene=gene,gene.func = func,gene.exonicfunc=exonicfunc))
  vcf_chunk = readVcf(vcf, "hg38")
}
tumor_ids = colnames(AD_matrix)[colnames(AD_matrix)!=normal_id]

#get_CN <-function(mutid, CNV_data, type="tumor"){
mut_starts = t( sapply( strsplit( as.character(tumors_reformated4pyclone[[1]][,"mutation_id"]) , ":|_",perl=T) , function(x) x[1:2]) )
mut_starts = data.frame(chr=as.character(mut_starts[,1]) , start=as.numeric(mut_starts[,2]), stringsAsFactors = F )
for(j in 1:nrow(CNV_tumors) ){
  ids = which(mut_starts$chr==CNV_tumors$CHR[j] & mut_starts$start >= CNV_tumors$START[j]  & mut_starts$start <= CNV_tumors$END[j] )
  if(length(ids)>0){
    cn_tmp = as.numeric(strsplit( as.character(CNV_tumors$cn_clone1[j]), split="\\|")[[1]])
    #if(max(cn_tmp)==0)  cn_tmp = as.numeric(strsplit( as.character(CNV_tumors$cn_clone2[j]), split="\\|")[[1]])
    for(k in 1:length(tumors_reformated4pyclone) ) tumors_reformated4pyclone[[k]][ids,c("normal_cn" ,"minor_cn" , "major_cn")] = matrix(c(2,min(cn_tmp),max(cn_tmp)) , ncol=3,
                                                                                                                                      nrow=length(ids),byrow = T)
  }
}

print(dim(tumors_reformated4pyclone[[1]]))

## NAs to check
# check random one
NAmut.maj  = which( is.na( tumors_reformated4pyclone[[1]]$major_cn ) )
NAmut.norm = which( is.na( tumors_reformated4pyclone[[1]]$normal_cn ) ) # should be the same as NAmut.maj
muttest = mut_starts[sample(NAmut.maj,1),] 
CNV_tumors[which(CNV_tumors$CHR==muttest$chr & CNV_tumors$START <= muttest$start & CNV_tumors$END >= muttest$start  ),]

# remove NAs
if(length(NAmut.norm)>0){ for(i in 1:ncol(rc)) tumors_reformated4pyclone[[i]] = tumors_reformated4pyclone[[i]][-NAmut.norm,] }

# remove major copy = 0
zeromut.maj  = which( tumors_reformated4pyclone[[1]]$major_cn==0 )
if(length(zeromut.maj)>0){ for(i in 1:ncol(rc)) tumors_reformated4pyclone[[i]] = tumors_reformated4pyclone[[i]][-zeromut.maj,] }

for(i in 1:length(tumors_reformated4pyclone)) write.table(tumors_reformated4pyclone[[i]], file=paste(output_folder,"/",tumor_ids[i],"_reformated4pyclone.csv",sep=""), quote=F, sep="\t", row.names = F)
# with DP > 100 in all samples
DPmat = c()
for(i in 1:length(tumors_reformated4pyclone)) DPmat = cbind(DPmat , rowSums(tumors_reformated4pyclone[[i]][,2:3]) )

DPthres = min_dp
for(i in 1:length(tumors_reformated4pyclone)) write.table(tumors_reformated4pyclone[[i]][which( rowSums(DPmat>DPthres)==ncol(DPmat) ),], file=paste(output_folder,"/",tumor_ids[i],"_minDP",DPthres,"_reformated4pyclone.csv",sep=""), quote=F, sep="\t", row.names = F)

# from list of genes
select.genes = c("TP53", "RB1", "STK11", "KEAP1", "MEN1" , "EIF1AX" , "ARID1A", "KRAS","NRAS","HRAS",
                 "MYC","MYCL1","MYCN","SOX2","CCNE1","RICTOR","NKX2-1","FGFR1","CDKN2A",
                 "NOTCH1","NOTCH2","NOTCH3","NOTCH4",
                 "PTEN","PI3KCA","ADAMTS12","ADAMTS2","GAS7","NTM","SMARCA2","NTRK2","NTRK3","LAMA1","PCLO","MEGF8")


library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)

for(i in 1:(length(tumors_reformated4pyclone)-1) ){
  for(j in (i+1):length(tumors_reformated4pyclone) ){
  dat <- tibble(T1 = as.numeric(tumors_reformated4pyclone[[i]]$var_counts)/(as.numeric(tumors_reformated4pyclone[[i]]$ref_counts) +  as.numeric(tumors_reformated4pyclone[[i]]$var_counts) ),
              T2 = as.numeric(tumors_reformated4pyclone[[j]]$var_counts)/(as.numeric(tumors_reformated4pyclone[[j]]$ref_counts) +  as.numeric(tumors_reformated4pyclone[[j]]$var_counts) ) )

  plot_ITH = ggplot(data = dat, mapping = aes(x = T1, y = T2)) + geom_pointdensity() +  scale_color_viridis() + #coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    theme_minimal() + labs(x=paste0("Sample ",i," (",tumor_ids[i],")"),y=paste0("Sample ",j," (",tumor_ids[j],")"))+ coord_fixed(xlim = c(0,1),ylim=c(0,1)) + 
    geom_point(data = dat[tumors_reformated4pyclone[[1]]$gene%in%select.genes & tumors_reformated4pyclone[[1]]$gene.func=="exonic",],size=2, color="red", shape=16) +
    geom_text(data = dat[tumors_reformated4pyclone[[1]]$gene%in%select.genes & tumors_reformated4pyclone[[1]]$gene.func=="exonic",] , color = "red", vjust = "inward", hjust = "inward",
              label = tumors_reformated4pyclone[[1]]$gene[tumors_reformated4pyclone[[1]]$gene%in%select.genes & tumors_reformated4pyclone[[1]]$gene.func=="exonic"] )

  svg(paste0(output_folder,"/",tumor_ids[i],"_",tumor_ids[j],"_scatterplot_smallVariants.svg"),h=4,w=4)
  print( plot_ITH  )
  dev.off()
  }
}
  
print(dim(tumors_reformated4pyclone[[1]]))

#for(i in 1:length(tumors_reformated4pyclone)) write.table(tumors_reformated4pyclone[[i]], 
#                                                          file=paste(output_folder,"/",tumor_ids[i],"_all_reformated4pyclone.csv",sep=""), quote=F, sep="\t", row.names = F)
# exonic
for(i in 1:length(tumors_reformated4pyclone)) write.table(tumors_reformated4pyclone[[i]][which( tumors_reformated4pyclone[[i]]$gene.exonicfunc%in%c("frameshift_deletion", "nonframeshift_deletion", "frameshift_substitution","nonframeshift_substitution", "nonsynonymous_SNV" , "stopgain","stoploss") ),], 
                                                          file=paste(output_folder,"/",tumor_ids[i],"_exonic_reformated4pyclone.csv",sep=""), quote=F, sep="\t", row.names = F)
