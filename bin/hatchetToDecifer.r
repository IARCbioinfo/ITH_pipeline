args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  return(res)
}

argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$help)) {help=FALSE} else {help=TRUE}

if(is.null(args$VCF) | is.null(args$CNV) | is.null(args$normal_id) | is.null(args$tumors_id) | help) {
  cat("

      hatchetToDecifer.r: build a tsv file to input to decifer from one multisample VCF from mutect2 and one bbc.ucn file from hatchet 

      Mandatory arguments:
      --VCF=file_name             - multisample VCF file from mutect2
      --CNV=file_name             - CNV file from hatchet (best.bbc.ucn file)
      --normal_id=id              - ID of the normal sample in VCF
      --tumors_id=id              - ID of the tumor samples in VCF, separated with comma (e.g. --tumors_id=T1,T2)

      Optional arguments:
      --min_dp=value              - minimum dp to consider a mutation (default: 20)
      --output_folder=path        - folder to save the results (default: .)
      --output_file=file          - output file (default: NORMAL_ID_reformated4decifer.tsv)
      --chunk_size=value          - size of reading (default: 10000)
      --help                      - print this text \n\n")
  q(save="no")
}

library(VariantAnnotation)

if(is.null(args$output_folder)) {output_folder="."} else {output_folder=args$output_folder}
if(is.null(args$chunk_size)) {chunk_size = 10000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$min_dp)) {min_dp = 20} else {min_dp = as.numeric(args$min_dp)}

normal_id = args$normal_id
tumors_id = unlist(strsplit(args$tumors_id,","))
if(is.null(args$output_file)) {output_file=paste(output_folder,"/",normal_id,"_reformated4decifer.tsv",sep="")} else {
  output_file=paste(output_folder,args$output_file,sep="/")}

CNV = read.table(args$CNV, quote="\"", stringsAsFactors=F, sep="\t", header=T, comment.char = "")
if(!file.exists(paste(args$VCF,".tbi",sep=""))) stop("Please index you input VCF file with tabix (tabix -p vcf file.vcf.gz)")

vcf <- open(VcfFile(args$VCF,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf, "hg19")

while(dim(vcf_chunk)[1] != 0) {
  AD_matrix = geno(vcf_chunk,"AD")
  GT_matrix = geno(vcf_chunk,"GT")
  DP_matrix = geno(vcf_chunk,"DP")
  tumors_mutid = as.numeric(which( (GT_matrix[,normal_id]=="0/0" | GT_matrix[,normal_id]=="0|0") &
                                     apply(DP_matrix[,tumors_id],1,function(x){min(x)})>min_dp &
                                     grepl("chr", rownames(vcf_chunk))))
  mutid=rownames(AD_matrix)
  # for the tumor 1 sample
  for(t in tumors_id){
    rc = unlist(lapply(AD_matrix[tumors_mutid,t], function(x) x[1])) # reference counts
    vc = unlist(lapply(AD_matrix[tumors_mutid,t], function(x) x[2])) # variant counts
    if(!exists("reformated4decifer")){
      reformated4decifer = data.frame(sample_index=match(t,tumors_id)-1, sample_label=t, character_label=mutid[tumors_mutid], 
                                      ref=rc, var=vc, row.names = NULL)
    } else {
      reformated4decifer = rbind(reformated4decifer,
                  data.frame(sample_index=match(t,tumors_id)-1, sample_label=t, character_label=mutid[tumors_mutid],
                             ref=rc, var=vc, row.names = NULL))
    }
  }
  vcf_chunk = readVcf(vcf, "hg19")
}
reformated4decifer$character_index = match(reformated4decifer$character_label, unique(reformated4decifer$character_label))
reformated4decifer = reformated4decifer[c("sample_index","sample_label","character_index","character_label","ref","var")]
x = colnames(CNV)[length(colnames(CNV))] ; n_clone_hatchet = as.numeric(substr(x,nchar(x),nchar(x)))

write.table(paste(length(tumors_id)," #samples",sep=""), file=output_file, append=T, quote = F, col.names = F, row.names = F, sep = "\t")
write.table(paste(length(unique(reformated4decifer$character_label))," #characters",sep=""), file=output_file, append=T,
            quote = F, col.names = F, row.names = F, sep = "\t")
write.table("#sample_index\tsample_label\tcharacter_index\tcharacter_label\tref\tvar", 
            file=output_file, append=T, quote = F, col.names = F, row.names = F, sep = "\t")

for(i in 1:nrow(reformated4decifer)){
  mutid = as.character(reformated4decifer[i,"character_label"])
  mut_start = as.numeric(strsplit(strsplit(mutid, ":")[[1]][2],"_")[[1]][1])
  mut_chr = strsplit(mutid, ":")[[1]][1]
  mut_tumor = as.character(reformated4decifer[i,"sample_label"])
  id = which(CNV[,"X.CHR"] == mut_chr & CNV[,"START"] <= mut_start & CNV[,"END"] >= mut_start & CNV$SAMPLE == mut_tumor)
  if(length(id)==0){
    dat = cbind(reformated4decifer[i,], 1, 1, 1)
  } else {
    dat = reformated4decifer[i,]
    for(clone in 1:n_clone_hatchet){
      prop = CNV[id,paste("u_clone",clone,sep="")] * (1/(1-CNV[id,"u_normal"]))
      mat = as.numeric(unlist(strsplit(CNV[id,paste("cn_clone",clone,sep="")],"|"))[1])
      pat = as.numeric(unlist(strsplit(CNV[id,paste("cn_clone",clone,sep="")],"|"))[3])
      dat = cbind(dat, mat, pat, prop)
    }
  }
  write.table(dat, file=output_file, append=T, quote = F, col.names = F, row.names = F, sep = "\t")
}

