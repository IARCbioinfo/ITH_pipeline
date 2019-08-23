args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  return(res)
}

argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$help)) {help=FALSE} else {help=TRUE}

if(is.null(args$VCF) | is.null(args$CNV_tumor1) | is.null(args$CNV_tumor2) | is.null(args$normal_id) | is.null(args$tumor1_id)
   | is.null(args$tumor2_id)| help) {
  cat("

      reformat4pyclone.r: build a tsv file to input to pyclone from one multisample VCF from mutect2 and one CNV.txt file from facets 

      Mandatory arguments:
      --VCF=file_name             - multisample VCF file from mutect2
      --CNV_tumors=file_name      - CNV file from hatchet for tumors
      --normal_id=id              - ID of the normal sample in VCF
      --tumor1_id=id              - ID of the tumor 1 sample in VCF
      --tumor2_id=id              - ID of the tumor 2 sample in VCF

      Optional arguments:
      --CNV_normal=file_name      - CNV file from facets tool for normal
      --min_dp=value              - minimum dp to consider a mutation (default: 20)
      --min_cnv_length=value      - minimum length in pb to consider a copy number (default:1)
      --output_folder=path        - folder to save the models (default: .)
      --chunk_size=value          - size of reading (default: 1000)
      --kept_genes=list           - list of kept genes separated by commas
      --help                      - print this text \n\n")
  q(save="no")
}

library(VariantAnnotation)

if(is.null(args$output_folder)) {output_folder="."} else {output_folder=args$output_folder}
if(is.null(args$chunk_size)) {chunk_size = 100000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$min_dp)) {min_dp = 20} else {min_dp = as.numeric(args$min_dp)}
if(is.null(args$min_cnv_length)) {min_cnv_length = 1} else {min_cnv_length = as.numeric(args$min_cnv_length)}
if(is.null(args$kept_genes)) {kept_genes = c()} else {kept_genes = unlist(strsplit(args$kept_genes,","))}

normal_id = args$normal_id


CNV_tumors = read.table(args$CNV_tumors ,comment.char = "", header = T)

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
  muts = unlist(lapply(rownames(vcf_chunk), function(s) unlist(strsplit(s,"_"))[[2]]))
  vcf_chunk = vcf_chunk[which(nchar(muts)==3),]
  AD_matrix = geno(vcf_chunk,"AD")
  GT_matrix = geno(vcf_chunk,"GT")
  DP_matrix = geno(vcf_chunk,"DP")
  gene = unlist(info(vcf_chunk)$Gene.refGene)
  tumors_mutid = as.numeric(which( (GT_matrix[,normal_id]=="0/0" | GT_matrix[,normal_id]=="0|0") &
                                                ( rowMeans(DP_matrix[,colnames(DP_matrix)!=normal_id]>min_dp) ==1
                                                  | gene %in% kept_genes)))

  mutid=rownames(AD_matrix)
  # for the tumor 1 sample
  rc = sapply(1:(ncol(AD_matrix)-1), function(i) matrix(unlist( AD_matrix[tumors_mutid,colnames(AD_matrix)!=normal_id][,i] ),ncol=ncol(AD_matrix)-1,byrow=T)[,1] ) #unlist(lapply( tumors_mutid, function(i) AD_matrix[i,colnames(DP_matrix)!=normal_id] ))
  vc = sapply(1:(ncol(AD_matrix)-1), function(i) matrix(unlist( AD_matrix[tumors_mutid,colnames(AD_matrix)!=normal_id][,i] ),ncol=ncol(AD_matrix)-1,byrow=T)[,2] )
  gene = unlist(info(vcf_chunk)$Gene.refGene)[tumors_mutid]
  gene = unlist(lapply(gene, function(g) unlist(strsplit(g,"[\\]"))[1]))
  
  funcl     = c(funcl, unlist(info(vcf_chunk)$Func.refGene)[tumors_mutid] )
  exonicfuncl = c(exonicfuncl, unlist(info(vcf_chunk)$ExonicFunc.refGene)[tumors_mutid] )
  for(i in 1:ncol(rc)) tumors_reformated4pyclone[[i]] = rbind(tumors_reformated4pyclone[[i]],
           data.frame(mutation_id=mutid[tumors_mutid], ref_counts=rc[,i], var_counts=vc[,i], normal_cn=NA, minor_cn=NA, major_cn=NA, gene=gene))
  vcf_chunk = readVcf(vcf, "hg38")
}
tumor_ids = colnames(AD_matrix)[colnames(AD_matrix)!=normal_id]

#get_CN <-function(mutid, CNV_data, type="tumor"){
mut_starts = t( sapply( strsplit( as.character(tumors_reformated4pyclone[[1]][,"mutation_id"]) , ":|_",perl=T) , function(x) x[1:2]) )
mut_starts = data.frame(chr=as.character(mut_starts[,1]) , start=as.numeric(mut_starts[,2]), stringsAsFactors = F )
for(j in 1:nrow(CNV_tumors) ){
  ids = which(mut_starts$chr==CNV_data$CHR[j] & mut_starts$start >= CNV_data$START[j]  & mut_starts$start <= CNV_data$END[j] )
  if(length(ids)>0){
    cn_tmp = as.numeric(strsplit( as.character(CNV_tumors$cn_clone1[j]), split="\\|")[[1]])
    #if(max(cn_tmp)==0)  cn_tmp = as.numeric(strsplit( as.character(CNV_tumors$cn_clone2[j]), split="\\|")[[1]])
    for(k in 1:length(tumors_reformated4pyclone) ) tumors_reformated4pyclone[[k]][ids,c("normal_cn" ,"minor_cn" , "major_cn")] = matrix(c(2,min(cn_tmp),max(cn_tmp)) , ncol=3,
                                                                                                                                      nrow=length(ids),byrow = T)
  }
}
## NAs to check
# check random one
NAmut.maj  = which( is.na( tumors_reformated4pyclone[[1]]$major_cn ) )
NAmut.norm = which( is.na( tumors_reformated4pyclone[[1]]$normal_cn ) ) # should be the same as NAmut.maj
muttest = mut_starts[sample(NAmut.maj,1),] 
which(CNV_data$CHR==muttest$chr & CNV_data$START <= muttest$start & CNV_data$END >= muttest$start  )

# remove NAs
for(i in 1:ncol(rc)) tumors_reformated4pyclone[[i]] = tumors_reformated4pyclone[[i]][-NAmut.norm,]

# remove major copy = 0
zeromut.maj  = which( tumors_reformated4pyclone[[1]]$major_cn==0 )
for(i in 1:ncol(rc)) tumors_reformated4pyclone[[i]] = tumors_reformated4pyclone[[i]][-zeromut.maj,]

for(i in 1:length(tumors_reformated4pyclone)) write.table(tumors_reformated4pyclone[[i]], file=paste(output_folder,"/",tumor_ids[i],"_reformated4pyclone.csv",sep=""), quote=F, sep="\t", row.names = F)
