#set up  R options to increase the max heap size of JVM  
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
#  restartSession()
#rm(list=ls())
# 
# Target human Population from : 
# 1. Normal Level Natural Radiation Area (NLNRA) (<1.5 mGy/yr)
# 2. High Level Natural Radiation Area (HLNRA) (>1.5mGy/yr)
# 
# Blood samples (Peripheral Blood Mononuclear Cells) were taken from 6 NLNRA individuals and 6  HLNRA individuals. 
# Thsi study aims to analyse the RNA sequencing data obained from these samples for differential expression among groups :
#   1. NLNRA (+Gy) vs NLNRA (basal) : comparison_type:normal 
# 2. HLNRA (+Gy) vs HLNRA (basal) : comparison_type:high
# 3. HLNRA vs NLNRA : comparison_type:basal
# 4. HLNRA (+Gy) vs NLNRA (+Gy) : comparison_type:exposed
# 
# The RNA seq-data was adapter trimmed,checked for quality, aligned to the reference genome Hg38 using HISAT2(a splice aware aligner). The alignment files were mapped to genes. featurecounts,a R program using a Hg38 GTF file gives counts of reads to gene regions identified by the NCBI Entrez ID for each of the sample. The code starts by inputing the sample files( bam format) to featurecounts .


comparison_type<<-"exposed" # one of the four comparison types enlisted above
fc<<- as.double(2.0) # fold change at which differential expression of genes is observed at
lnfc=log2(fc)
setwd(paste0("/media/pavitra/Backup Plus3/final_deseq2_entrez_ids/FC_",fc,"/",comparison_type,"/"))

dir()
#required_packages
libraries<- function(){
  #library(BiocManager)
  library(BiocVersion)
  library(DESeq2)
  library(ggplot2)
  library(plotly)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(tidyverse)
  library("BiocParallel")
  register(MulticoreParam(8))
  library(Rsubread)
  library(xlsx)
  library(readr)
  library(rstudioapi) #for restarting R session command 
}


#comparison_type= one of normal,high, basal , exposed
call_title <-function(comparison_type){
  if(comparison_type=="normal"){
    return("NLNRA2GY_vs_NLNRA")
  } else if (comparison_type=="high"){
    return("HLNRA2GY_vs_HLNRA")
  }
  else if (comparison_type=="basal"){
    return("HLNRA_vs_NLNRA")
  } else if (comparison_type=="exposed"){
    return ("HLNRA2GY_vs_NLNRA2GY")
  } else {}
  
}

gp_num <-function(comparison_type){
  if(comparison_type=="normal"){
    return("I")
  } else if (comparison_type=="high"){
    return("II")
  }
  else if (comparison_type=="basal"){
    return("III")
  } else if (comparison_type=="exposed"){
    return ("IV")
  } else {}
  
}

featurecounts <-function(sampleFiles){
  fc=featureCounts(files=sampleFiles,annot.inbuilt ="hg38",GTF.featureType = "gene_id",chrAliases = "/media/pavitra/Backup Plus/all_bam_files/hg38/refseq_to_chr.txt",isPairedEnd = TRUE,nthreads = 8)
  return(fc)
}
create_sample_table<-function(comparison_type){
  #comparison type is one of normal ; high ; basal ; exposed
  sample_table_file_name=paste0(comparison_type,"_sample_table_hg38_2.tab")
  sample_table=read.table(sample_table_file_name)
  return(sample_table)
}

count_matrix <- function(comparison_type) {
  
  #creating the count matrix to feed to deseq2
  featurecounts_outfile = paste0("counts_",comparison_type,"_hg38.txt")
  a=read.delim(featurecounts_outfile,sep=" ")
  counts=a[,(3:length(a))]
  head(counts)
  return(counts)
  
}
run_deseq <-function(count_matrix, sampleinfo_table){
  dds=DESeqDataSetFromMatrix(countData = count_matrix,colData = sampleinfo_table,design= ~ condition)
  #for getting the normalized counts
  d<- estimateSizeFactors(dds)
  s=sizeFactors(d)
  normlzd_dds= as.data.frame(counts(d, normalized=T))
  colnames(normlzd_dds)=table$name
  
  outfile=paste0(call_title(comparison_type),"_normalized_counts.tsv")
  write.table(normlzd_dds,outfile,sep = "\t",row.names = TRUE,col.names = NA)
  outfile=paste0(call_title(comparison_type),"_normalized_counts.xlsx")
  write.xlsx2(normlzd_dds, file = outfile, sheetName = "Sheet1", col.names = TRUE, row.names = TRUE,append = FALSE)
  
  dds<- DESeq(dds)
  return(dds)
  
}
deseq2_results <-function(dds,comparison_type){
  if(comparison_type=="normal"||comparison_type=="high"){
    res=results(dds,contrast = c("condition","exposed","basal"))
    print(summary(res))
    # df=write_deseq2_results(res,comparison_type)
    return(res)
    
  } else {
    res=results(dds,contrast = c("condition","high","normal"))
    print(summary(res))
    #   df=write_deseq2_results(res,comparison_type)
    return(res)
  }
}
write_deseq2_results<-function(res,comparison_type){
  #res is results from deseq
  df=data.frame(res)
  df=rownames_to_column(df,var='EntrezGeneID')
  
  outfile=paste0(call_title(comparison_type),"_deseq2_results.tsv")
  write_tsv(df,outfile)
  outfile=paste0(call_title(comparison_type),"_deseq2_results.xlsx")
  write.xlsx2(df, file = outfile, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
  return(df)
}
thresholds<-function(res_df,res,padj,lnfc,comparison_type){
  print(head(res_df))
  filter1=filter(res_df,complete.cases(res))
  print("after removing low counts and outliers")
  print(dim(filter1))
  
  print("after putting filter for padjusted value<0.05")
  filter2_padj=filter(data.frame(filter1),padj<0.05)
  write_tsv(filter2_padj,paste0(call_title(comparison_type),"_padj_lt0.05.tsv"))
  write.xlsx2(filter2_padj,paste0(call_title(comparison_type),"_padj_lt0.05.xlsx"), col.names = TRUE, row.names = FALSE, append = FALSE)
  
  print(dim(filter2_padj))
  
  print(paste0("after putting filter for fc >"),fc)
  filter3_lnfc=filter(filter2_padj,abs(log2FoldChange)>lnfc)
  print(dim(filter3_lnfc))
  
  return(filter3_lnfc)
}

up_regulated_genes <-function(df,fc,comparison_type){
  #upregulated genes 
  lnfc=log2(fc)
  up=filter(df,log2FoldChange>lnfc)
  print(count(up))
  write_tsv(up,paste0("upregulated_",comparison_type,"_fc_",fc,".tsv"))
  write.xlsx2(up, file = paste0("upregulated_",comparison_type,"_fc_",fc,".xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
  return(up)
  
}

#downregulated genes
down_regulated_genes <-function(df,fc,comparison_type){
  #upregulated genes 
  lnfc=log2(fc)
  down=filter(df,log2FoldChange < -lnfc)
  print(count(down))
  write_tsv(down,paste0("down_",comparison_type,"_fc_",fc,".tsv"))
  write.xlsx2(down, file = paste0("down_",comparison_type,"_fc_",fc,".xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
  return(down)
}

#annotation of a list of entrez ids
Annotation<-function(df){
  a=AnnotationDbi::select(org.Hs.eg.db,keys=df$EntrezGeneID,col = c("SYMBOL", "GENENAME"),keytype = "ENTREZID")
  annotated_df=left_join(df,a,by=c('EntrezGeneID'='ENTREZID'))
  write_tsv(annotated_df,paste0("annotated_",comparison_type,"_deseq2_results_unfiltered_fc_1.3.tsv"))
  write.xlsx2(annotated_df,paste0("annotated_",comparison_type,"_deseq2_results_unfiltered_fc_1.3.xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
  return(annotated_df)
}

#PCA and Dispersion Plot
pca <-function(des){
  print(comparison_type)
  vst=varianceStabilizingTransformation(des)
  p <-plotPCA(vst,intgroup=('condition'),)
  # p <-p + geom_text(aes(label = des$name),position = position_nudge(y=1))
  t=call_title(comparison_type)
  p <- p + ggtitle(t)
  p<-p + theme(plot.title = element_text(color="black", size=18, face="bold.italic"),axis.text.x = element_text(face="bold", size=15),axis.text.y = element_text(face="bold", size=15),axis.title.x=element_text(face="bold", size=15),axis.title.y=element_text(face="bold",size=15),legend.title = element_text(face="bold",size=15),legend.text = element_text(face="bold",size=15))+geom_point(size=5)
  ggsave(paste0("pca_",comparison_type,".png"),device="png",width = 15,height = 8)
  
  flname=paste0("dispersion_",comparison_type,".png")
  png(flname,height = 548,width = 978)
  plotDispEsts(des,cex.axis=1.5,cex.lab=1.5,cex=1.3)
  dev.off()
  
}

#Volcano Plot: log transformed adjusted p-values plotted on the y-axis and log2 fold change values on the x-axis. 
vp<-function(annotated_whole,padj,fc,n_genes){
  lnfc=log2(fc)
  print(lnfc)
  #annotated_whole$test = annotated_whole$padj<padj & abs(annotated_whole$log2FoldChange)>lnfc
  t=call_title(comparison_type)
  t=paste0(t,"\npadj< ",padj,", fc>",fc,"\nNumber of DE genes: ",n_genes)
  annotated_whole$diffexpressed <-"NO"
  annotated_whole$diffexpressed[annotated_whole$log2FoldChange>lnfc & annotated_whole$padj<padj]<-"UP"
  annotated_whole$diffexpressed[annotated_whole$log2FoldChange < -lnfc & annotated_whole$padj<padj]<-"DOWN"
  print(head(annotated_whole$log2FoldChange,50))
  print(head(annotated_whole$padj,50))
  print(head(annotated_whole$diffexpressed,50))
  g=ggplot(annotated_whole,aes(x=log2FoldChange,y=-log10(padj),col=diffexpressed)) + 
    #geom_point(aes(color=test),size=1.5,alpha=1.0) + 
    geom_point()+
    scale_color_manual(values=c('red','black','blue'))+
    geom_vline(xintercept=lnfc,color='darkgreen',linetype=1,aes(label="1"))+
    geom_vline(xintercept=-lnfc,color='darkgreen',linetype=1)+
    geom_hline(yintercept=-log10(padj),color='blue',linetype=1)+
    theme_bw()+
    theme(plot.title = element_text(color="black", size=18, face="bold.italic"),axis.text.x = element_text(face="bold", size=16),axis.text.y = element_text(face="bold", size=16),axis.title.x=element_text(face="bold", size=16),axis.title.y=element_text(face="bold",size=16),legend.title=element_text(face="bold",size=14),legend.text=element_text(face="bold",size=13))
  g <- g + ggtitle(t)
  
  
  ggplotly(g)  
  ggsave(paste0("volcanoplot_",comparison_type,"_",fc,".png"), device="png",width = 15,height = 8)
}

goplots <-function(go_res,ont,orientation,comparison_type,fc){
  if(ont=='all'){
    print("yes the ont is all")
    ont="GO_terms"
  }
  d=enrichplot::dotplot(go_res,showCategory=25,title=paste0("\n",ont," enriched for ",orientation," genes at fold change ",fc," in Group ",gp_num(comparison_type),":",call_title(comparison_type),"\n\n"))
  d=d+scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  d=d+theme(plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.9,vjust = 0.9),axis.text.y= element_text(size=16,face="bold"),axis.text.x= element_text(size=16,face="bold"),axis.title.x = element_text(size=16,face="bold"),legend.title=element_text(face="bold",size=14),legend.text=element_text(face="bold",size=13))
  ggsave(paste0("dotplots/",comparison_type,"_",orientation,"_",ont,"_","fc_",fc,"_dotplot.png"), device="png",width = 18,height =15)
  enrichplot::heatplot(go_res)
  g=goplot(go_res,title=paste0(call_title(comparison_type)," ",orientation," ",ont," fc_",fc))
  g=g+ggtitle(paste0(call_title(comparison_type)," ",orientation," ",ont))
  # ggsave(paste0("goplots/",comparison_type,"_",orientation,"_",ont,"_","fc_",fc,"_goplot.png"), device="png",width = 18,height = 6,title=paste0(call_title(comparison_type)," ",orientation," ",ont))
}

cnetplots <- function(go_res,ont,orientation,comparison_type,fc){
  fold_changes=res$log2FoldChange
  # dev.off()
  #png(paste0(comparison_type,"_",orientation,"_",ont,"_cnetplot.png"))
  d= cnetplot(go_res,showCategory=5,foldChange=fold_changes,colorEdge=TRUE,cex_label_category=2,cex_label_gene=0.8)
  d=d+ggtitle(paste0(comparison_type,"_",orientation,"_",ont,"_fc_",fc,"cnetplot (top5)"))
  d=d+theme(legend.position = "none", legend.title=element_text(size=16, face = "bold"), plot.title = element_text(size=20, face = "bold"))
  ggsave(paste0("cnetplots/",comparison_type,"_",orientation,"_",ont,"_fc_",fc,"_cnetplot_5.png"), device="png",width = 30,height = 8,limitsize = FALSE)
  
  d= cnetplot(go_res,showCategory=10,foldChange=fold_changes,colorEdge=TRUE,cex_label_category=2.5,cex_label_gene=2,cex_category =3,cex_gene=2)
  d=d+ggtitle(paste0(comparison_type,"_",orientation,"_",ont,"_fc_",fc,"cnetplot (top 10) "))
  d=d+theme(legend.position = "none", plot.title = element_text(size=40, face = "bold"))
  ggsave(paste0("cnetplots/",comparison_type,"_",orientation,"_",ont,"_fc_",fc,"_cnetplot_10.png"), device="png",width =40 ,height = 20,limitsize = FALSE)
  
  d= cnetplot(go_res,showCategory=25,foldChange=fold_changes,colorEdge=TRUE,cex_label_category=5,cex_label_gene=3.5,cex_category =4.5,cex_gene=3)
  # d=d+scale......(labels=function(x) str_wrap(x, width=45))
  #scale_colour_discrete(labels = function(x) str_wrap(x, width = 5))
  d=d+ggtitle(paste0(comparison_type,"_",orientation,"_",ont,"_fc_",fc,"cnetplot (top 25) "))
  d=d+theme(legend.position = "none", plot.title = element_text(size=30, face = "bold"))
  
  ggsave(paste0("cnetplots/",comparison_type,"_",orientation,"_",ont,"_fc_",fc,"_cnetplot_25.png"), device="png",width =60 ,height = 50,limitsize = FALSE)
  
}

#Gene Ontology analysis using enrichGO
ego <-function(df,uni_df,ont,orientation,comparison_type){
  
  print(count(df))
  print(count(uni_df))
  o=clusterProfiler::enrichGO(gene = df$EntrezGeneID,
                              OrgDb = org.Hs.eg.db,
                              keyType = "ENTREZID",
                              ont =ont,
                              universe=uni_df$EntrezGeneID,
                              pvalueCutoff = 0.1,
                              readable=TRUE
  )
  ont_result=o@result[, c(1,6)]
  # name=paste0(df)
  write_tsv(ont_result,paste0("revigo/tsv/rvgo_",comparison_type,"_",orientation,"_",ont,"_fc_",fc,".tsv"))
  write.xlsx2(ont_result, file = paste0("revigo/xlsx/rvgo_",comparison_type,"_",orientation,"_",ont,"_fc_",fc,".xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
  
  write_tsv(o@result,paste0("clusterprofiler/tsv/clusterprofiler@enrichgo_result_",comparison_type,"_",orientation,"_",ont,"_fc_",fc,".tsv"))
  write.xlsx2(o@result, file = paste0("clusterprofiler/xlsx/clusterprofiler@enrichgo_result_",comparison_type,"_",orientation,"_",ont,"_fc_",fc,".xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
  
  goplots(o,ont,orientation,comparison_type,fc)
  print("doing cnetplots")
  cnetplots(o,ont,orientation,comparison_type,fc)
  return(o)
}


libraries()

# get the meta information about the samples
# this could be:
#normal:"NLNRA2GY_vs_NLNRA"
#high: "HLNRA2GY_vs_HLNRA"
#basal: "HLNRA_vs_NLNRA"
#exposed: "HLNRA2GY_vs_NLNRA2GY"

# table will have the order in which the samples occur
table=create_sample_table(comparison_type)
table
#count matrix to be used for deseq
cm=count_matrix(comparison_type)

padj=0.05

c=cm
c=rownames_to_column(c,var='EntrezGeneID')
write_tsv(c,paste0("countmatrix_",comparison_type,".tsv"))
write.xlsx2(c, file = paste0("countmatrix_",comparison_type,".xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
#run DESeq
# defaults for the DESeq
des=run_deseq(cm,table)

dir()

res=deseq2_results(des,comparison_type)
#resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
DESeq2::summary(res)
unfiltered_df=write_deseq2_results(res,comparison_type)

n=call_title(comparison_type)
n
lnfc=log2(fc)
lnfc
#annotation of a list of entrez ids

annotated_whole=Annotation(unfiltered_df)

filtered=thresholds(annotated_whole,res,padj,lnfc,comparison_type)
sorted_padj=filtered[order(filtered$padj),]
#upregulated genes
upreg=up_regulated_genes(df = sorted_padj,comparison_type =comparison_type,fc = fc)

#downregulated genes
downreg=down_regulated_genes(df=sorted_padj,comparison_type =comparison_type,fc = fc)



FLname=paste0(n,"_padj_",padj,"_fc_",fc,".tsv")
write_tsv(sorted_padj,file = FLname)
FLname=paste0(n,"_padj_",padj,"_fc_",fc,".xlsx")
write.xlsx2(sorted_padj, file = FLname, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)

dir()
# check analysis

#principal component analysis
pca(des)
n_genes=nrow(sorted_padj)
n_genes
vp(unfiltered_df,padj,fc,n_genes)


# Gene Ontology Analysis

ont_list<- list("BP","MF","CC","all")


dir.create(paste0("clusterprofiler/tsv/"),recursive=TRUE)
dir.create(paste0("clusterprofiler/xlsx/"),recursive=TRUE)
dir.create(paste0("cnetplots"))
ego(upreg,uni_df = unfiltered_df,ont = "BP",orientation = "upregulated",comparison_type = comparison_type)
ego(downreg,uni_df = unfiltered_df,ont = "BP",orientation = "downregulated",comparison_type = comparison_type)
ego(sorted_padj,uni_df = unfiltered_df,ont = "BP",orientation ="differentially expressed",comparison_type =comparison_type)

ego(upreg,uni_df = unfiltered_df,ont = "MF",orientation = "upregulated",comparison_type = comparison_type)
ego(downreg,uni_df = unfiltered_df,ont = "MF",orientation = "downregulated",comparison_type = comparison_type)
ego(sorted_padj,uni_df = unfiltered_df,ont ="MF",orientation = "differentially expressed",comparison_type =comparison_type)


ego(upreg,uni_df = unfiltered_df,ont = "CC",orientation = "upregulated",comparison_type = comparison_type)
ego(downreg,uni_df = unfiltered_df,ont = "CC",orientation = "downregulated",comparison_type = comparison_type)
ego(sorted_padj,uni_df = unfiltered_df,ont = "CC",orientation = "differentially expressed",comparison_type =comparison_type)

ego(upreg,uni_df = unfiltered_df,ont = "all",orientation = "upregulated",comparison_type = comparison_type)
ego(downreg,uni_df = unfiltered_df,ont = "all",orientation = "downregulated",comparison_type = comparison_type)
ego(sorted_padj,uni_df = unfiltered_df,ont = "all",orientation = "differentially expressed",comparison_type =comparison_type)
#rm(list=ls())
#restartSession()



