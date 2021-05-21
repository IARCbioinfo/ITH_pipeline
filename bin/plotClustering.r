####################################################
##             Plot clusterings                   ##
####################################################

library(ggplot2)
library(readr)
library(ggpubr)
library(patchwork)
library(tibble)
library(stringr)
library(dplyr)
library(tidyr)

# get command line arguments
args = commandArgs(trailingOnly=TRUE) #segmentation file

# read bbc file
bbc = read_tsv(args[1])
ucn.files = grep(grep(list.files(".",pattern="bbc.ucn"),pattern = "best|chosen|n2|n3",value = T),pattern = ".png",invert = T,value = T)
ucn = lapply( ucn.files,read_tsv)
for(i in 1:length(ucn)){
  clonecols = which( str_detect(colnames(ucn[[i]]),"cn_clone") )
  for(j in clonecols){ #split CN column into numeric columns for major and minor CN
    res = as_tibble(t( sapply( ucn[[i]] %>% pull(colnames(ucn[[i]])[j]) , function(x) str_split(x,"\\|")[[1]] ) ))
    colnames(res) = paste0(colnames(ucn[[i]])[j],c(".major",".minor"))
    res = res %>% mutate_if(is.character,as.numeric)
    ucn[[i]] = bind_cols(ucn[[i]],res)
  }
  clonecols2  = which( str_detect(colnames(ucn[[i]]),"cn_clone[0-9]+.m") )
  uclonecols2 = which( str_detect(colnames(ucn[[i]]),"u_clone[0-9]+") )
  ucn[[i]]$RDR_predicted = (rowSums( (ucn[[i]][,clonecols2[seq(1,length(clonecols2),2)]]+
                                       ucn[[i]][,clonecols2[seq(2,length(clonecols2),2)]])*ucn[[i]][,uclonecols2] ) + 2*ucn[[i]]$u_normal)/2
  
  #put in long format for easier ggplot use
  ucn[[i]] = ucn[[i]]%>% pivot_longer(cols=colnames(ucn[[i]])[clonecols2],names_to = "allele",values_to = "cn_clone_numeric")
  ucn[[i]] = ucn[[i]] %>% mutate(clone=str_remove(str_remove(allele,".major|.minor"),"cn_") , allele=str_remove(allele,"cn_clone[0-9]+."))
}

# chromosome offsets
offset = c(0,  248956422,  491149951,  689445510,  879660065, 1061198324, 1232004303, 1391350276, 1536488912, 1674883629, 1808681051, 1943767673, 
           2077042982, 2191407310, 2298451028, 2400442217, 2490780562, 2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522, 3031042417,
           3088269832, 3088286401)

# plot
theme_hatchet <- theme_classic() 
gg_clust = ggplot(bbc,aes(x=0.5-BAF,y=RD,col=factor(CLUSTER) )) + geom_point() + guides(col=F) + theme_hatchet + grids() +facet_wrap(.~SAMPLE) 
gg_RDR = ggplot(bbc, aes(x=START+offset[as.numeric(factor(bbc$`#CHR`,levels=paste0("chr",1:22)))],
                         xend=END+offset[as.numeric(factor(bbc$`#CHR`,levels=paste0("chr",1:22)))],y=RD,yend=RD,col=factor(CLUSTER))) + geom_segment(size =2) + 
  geom_vline(xintercept = offset[1:22],col=rgb(0.5,0.5,0.5)) + 
  geom_text(data=tibble(chr=paste0("chr",1:22),pos=0.5*(offset[1:22]+offset[2:23]),SAMPLE=bbc$SAMPLE[1] ), 
            mapping = aes(x=pos,label=chr),inherit.aes = F, y=max(bbc$RD)*0.92,angle=90,col=rgb(0.5,0.5,0.5), size=2.5 ) + 
  facet_grid(SAMPLE~.) + coord_cartesian(expand = F) +
  guides(col=F)  + theme_hatchet + xlab("Genomic location") + ylab("Read Depth Ratio") + grids(axis = "y") + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
gg_BAF = ggplot(bbc, aes(x=START+offset[as.numeric(factor(bbc$`#CHR`,levels=paste0("chr",1:22)))],
                         xend=END+offset[as.numeric(factor(bbc$`#CHR`,levels=paste0("chr",1:22)))],y=BAF,yend=BAF,col=factor(CLUSTER))) + geom_segment(size =2) + 
  geom_vline(xintercept = offset[1:22],col=rgb(0.5,0.5,0.5)) + 
  geom_text(data=tibble(chr=paste0("chr",1:22),pos=0.5*(offset[1:22]+offset[2:23]),SAMPLE=bbc$SAMPLE[2] ), 
            mapping = aes(x=pos,label=chr),inherit.aes = F, y=min(bbc$BAF)+0.008,angle=90,col=rgb(0.5,0.5,0.5), size=2.5 ) + 
  facet_grid(SAMPLE~.) + coord_cartesian(expand = F) +
  guides(col=F)  + theme_hatchet + xlab("Genomic location") + ylab("B-allele frequency") + grids(axis = "y") + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

# plot clones
for(i in 1:length(ucn)){
  gg_profile <- ggplot(ucn[[i]], aes(x=START+offset[as.numeric(factor(ucn[[i]]$`#CHR`,levels=paste0("chr",1:22)))],
                xend=END+offset[as.numeric(factor(ucn[[i]]$`#CHR`,levels=paste0("chr",1:22)))],y=cn_clone_numeric,yend=cn_clone_numeric,col=factor(CLUSTER))) + 
    geom_segment(size =2) + 
    geom_vline(xintercept = offset[1:22],col=rgb(0.5,0.5,0.5)) + 
    geom_text(data=tibble(chr=paste0("chr",1:22),pos=0.5*(offset[1:22]+offset[2:23]),clone="clone1" ), 
            mapping = aes(x=pos,label=chr),inherit.aes = F, y=max(ucn[[i]]$cn_clone_numeric)*0.85,angle=90,col=rgb(0.5,0.5,0.5), size=2.5 ) + 
    facet_grid(clone~allele) + coord_cartesian(expand = F) +
    guides(col=F)  + theme_hatchet + xlab("Genomic location") + ylab("Predicted CN") + grids(axis = "y") + 
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  
  gg_RDR_pred = ggplot(ucn[[i]], aes(x=START+offset[as.numeric(factor(ucn[[i]]$`#CHR`,levels=paste0("chr",1:22)))],
                           xend=END+offset[as.numeric(factor(ucn[[i]]$`#CHR`,levels=paste0("chr",1:22)))],y=RDR_predicted,yend=RDR_predicted,col=factor(CLUSTER))) + geom_segment(size =2) + 
    geom_vline(xintercept = offset[1:22],col=rgb(0.5,0.5,0.5)) + 
    geom_text(data=tibble(chr=paste0("chr",1:22),pos=0.5*(offset[1:22]+offset[2:23]),SAMPLE=ucn[[i]]$SAMPLE[1] ), 
              mapping = aes(x=pos,label=chr),inherit.aes = F, y=max(ucn[[i]]$RDR_predicted)*0.93,angle=90,col=rgb(0.5,0.5,0.5), size=2.5 ) + 
    facet_grid(SAMPLE~.) + coord_cartesian(expand = F) +
    guides(col=F)  + theme_hatchet + xlab("Genomic location") + ylab("Predicted RDR") + grids(axis = "y") + 
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
 
  # get clone proportion table
  proptab = ucn[[i]] %>% dplyr::select(grep(colnames(ucn[[i]]),pattern="SAMPLE|u_",value = T)) %>% distinct()
  proptab = proptab %>% pivot_longer(2:ncol(proptab))
  proptab$value[proptab$SAMPLE==ucn[[i]]$SAMPLE[2]] = -proptab$value[proptab$SAMPLE==ucn[[i]]$SAMPLE[2]]
  
  gg_clones <- ggplot(proptab,aes(y=name, x=value,fill=name)) + geom_col(position=position_stack() ) + 
    theme_hatchet +  theme_void() + coord_cartesian(ylim=c(-0.25,length( unique(proptab$name) )), xlim=c(-1.1,1) )+
    geom_segment(data = tibble(start=c(0,0),end=c(-1,1)), mapping = aes(x=start,xend=end,y=0,yend=0), arrow=arrow(),inherit.aes = F) +
    geom_text(data = tibble(x=c(-0.5,0.5),y=c(-0.5,-0.5),label=unique(proptab$SAMPLE) ),aes(x=x,y=y,label=label),inherit.aes = F ,size=3) +
    geom_text(data = tibble(x=proptab$value[!duplicated(proptab$name)]-0.19,y=unique(proptab$name), label=str_remove(unique(proptab$name),"u_") ),
              aes(x=x,y=y,label=label,col=y),inherit.aes = F ,size=3) +
    grids(axis = "xy") + geom_vline(xintercept = 0,linetype="solid") + 
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(length(unique(proptab$name))-1,"Set1"),"lightgray")) + 
    scale_color_manual(values = c(RColorBrewer::brewer.pal(length(unique(proptab$name))-1,"Set1"),"lightgray")) + guides(fill=F,col=F)
   
  ggsave(paste0(args[2],ucn.files[i],"_profiles.png"), (gg_clust+gg_clones) / gg_BAF / gg_RDR / gg_RDR_pred / gg_profile +  
           plot_layout(heights = c(1, 2,2,2,3) , widths = c(1,2,2,2,2)) + plot_annotation(tag_levels = c('A')), height = 4*3.2,width = 4*2 )
}

# save to file
ggsave(paste0(args[2],args[1],"_QC.png"), gg_clust / gg_RDR / gg_BAF +  plot_layout(heights = c(1, 2,2), widths = c(1,2,2)) , height = 4*2.1,width = 4*1.5 )
ggsave(paste0(args[2],args[1],"_QC.pdf"), gg_clust / gg_RDR / gg_BAF +  plot_layout(heights = c(1, 2,2), widths = c(1,2,2)) , height = 4*2.1,width = 4*1.5 )
