



library(ontologyIndex)
library(openxlsx)
library(ComplexHeatmap)
library(MatrixGenerics)
library(circlize)
#########load go.obo file

GO_infor<-get_OBO(
  "/go.obo",###The file could be downloaded from https://geneontology.org/docs/download-ontology/
  #propagate_relationships = "is_a",
  extract_tags = "everything",
  merge_equivalent_terms = TRUE
)


##########load enriched results
tmp_enrichment.Workbook<-loadWorkbook("/862+13803.VS.NC.5.enrichment.xlsx")
prefix="862+13803.VS.NC.5.enrichment"
tmp_enrichment<-list()
for(i in names(tmp_enrichment.Workbook)){
  tmp_enrichment[[i]]=read.xlsx(tmp_enrichment.Workbook,sheet=i)
  rownames(tmp_enrichment[[i]])<-tmp_enrichment[[i]]$ID
}
tmp_gene_set_name=c()
for(k in names(tmp_enrichment)){
  
  tmp_enrichment[[k]][grep("kill",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("_NK_",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("_nk_",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("_NKCELL_",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("_nkcell_",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("natural_killer",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("NATURAL_KILL",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("natural_kill",tmp_enrichment[[k]]$Description),"NK_related"]="Y"
  tmp_enrichment[[k]][grep("_T_cell",tmp_enrichment[[k]]$Description),"T_cell_related"]="Y"
  tmp_enrichment[[k]][grep("_t_cell",tmp_enrichment[[k]]$Description),"T_cell_related"]="Y"
  tmp_enrichment[[k]][grep("tcell",tmp_enrichment[[k]]$Description),"T_cell_related"]="Y"
  tmp_enrichment[[k]][grep("Tcell",tmp_enrichment[[k]]$Description),"T_cell_related"]="Y"
  tmp_enrichment[[k]]<-tmp_enrichment[[k]][,-3]
  
  tmp_results<-tmp_enrichment[[k]][tmp_enrichment[[k]]$p.adjust<0.05,]
  tmp_gene_set_name<-unique(c(tmp_gene_set_name,rownames(tmp_results[tmp_results$gs_subcat%in%"GO:BP",])))
  
}


##########get P value of all terms in each cluster
p.value.df<-data.frame(matrix(data=NA,nrow=length(tmp_gene_set_name),ncol = 5))
rownames(p.value.df)<-tmp_gene_set_name
colnames(p.value.df)<-c("Cluster.1","Cluster.2","Cluster.3","Cluster.4","Cluster.5"   )
p.adj.df<-p.value.df
q.value.df<-p.value.df
for(k in names(tmp_enrichment)){
  p.value.df[,k]<-tmp_enrichment[[k]][rownames(p.value.df),"pvalue"]
  p.adj.df[,k]<-tmp_enrichment[[k]][rownames(p.adj.df), "p.adjust" ]
  q.value.df[,k]<-tmp_enrichment[[k]][rownames(q.value.df),"qvalue"]
}
p.value.df[is.na(p.value.df)]=1
p.adj.df[is.na(p.adj.df)]=1
q.value.df[is.na(q.value.df)]=1
tmp_p.value.df<-t(scale(t(log(p.value.df))))

############conducting term p-value based clustring to find out cluster-specific go terms
col_fun = colorRamp2(c(max(tmp_p.value.df),max(tmp_p.value.df)*3/5,max(tmp_p.value.df)*1/4,0,1/4*min(tmp_p.value.df),3/5*min(tmp_p.value.df),min(tmp_p.value.df)), c("#0088C3","#6FBCDB","#BCE5F1","#f3f3f3","#FFC482","#FF7B4C","#FF5831") )

enrichment_heatmap<-Heatmap(tmp_p.value.df,
                            name = "Scaled(log Pvalue)",
                            show_column_names = T,
                            col = col_fun,
                            show_row_names = F,
                            column_names_rot = 90,
                            row_km = 10,
                            #left_annotation   = ha,
                            cluster_rows = T,
                            cluster_columns = F,raster_device ="tiff",use_raster =T)
pdf(paste0(prefix,".heatmap.pdf"),height = 5,width = 3)

print(enrichment_heatmap)
dev.off()



############get clustered terms and their ancestor terms

enrichment_heatmap1 <- draw(enrichment_heatmap)

name_of_each_cluster_enrichment1<-row_order(enrichment_heatmap)
name_of_each_cluster_enrichment_pvalue1<-row_order(enrichment_heatmap)
name_of_each_cluster_enrichment=name_of_each_cluster_enrichment_pvalue=list()
for(tmp_cluster1 in names(name_of_each_cluster_enrichment1)){
  tmp_name<-paste0("Enrichment_cluster.",tmp_cluster1)
  name_of_each_cluster_enrichment[[tmp_name]]<-tmp_p.value.df[name_of_each_cluster_enrichment1[[tmp_cluster1]],]
  
  name_of_each_cluster_enrichment_pvalue[[tmp_name]]<-p.value.df[rownames(name_of_each_cluster_enrichment[[tmp_name]]),]
}




tmp_GO_list<-list()
tmp_GO_pvalue_list<-list()
tmp_GO_ancestors_list<-list()
tmp_GO_ancestors_summary_list<-list()
for(tmp_cluster in names(name_of_each_cluster_enrichment)){
  tmp_terms<-rownames(name_of_each_cluster_enrichment[[tmp_cluster]])
  tmp_terms<-gsub("GOBP_","",tmp_terms)
  tmp_terms<-gsub("_"," ",tmp_terms)
  tmp_terms<-tolower(tmp_terms)
  tmp_GO_list[[tmp_cluster]]<-data.frame("GO_ID"=names(GO_infor$name[GO_infor$name%in%tmp_terms]),
                                         "GO_name"=GO_infor$name[names(GO_infor$name[GO_infor$name%in%tmp_terms])])
  tmp_GO_list[[tmp_cluster]]$terms<-paste0("GOBP_",tmp_GO_list[[tmp_cluster]]$GO_name)
  tmp_GO_list[[tmp_cluster]]$terms<-toupper(tmp_GO_list[[tmp_cluster]]$terms)
  tmp_GO_list[[tmp_cluster]]$terms<-gsub(" ","_",tmp_GO_list[[tmp_cluster]]$terms)
  tmp_GO_pvalue_list[[tmp_cluster]]=tmp_GO_list[[tmp_cluster]]
  
  
  tmp_GO_list[[tmp_cluster]]<-cbind(tmp_GO_list[[tmp_cluster]],name_of_each_cluster_enrichment[[tmp_cluster]][tmp_GO_list[[tmp_cluster]]$terms,])
  tmp_GO_pvalue_list[[tmp_cluster]]<-cbind(tmp_GO_pvalue_list[[tmp_cluster]],name_of_each_cluster_enrichment_pvalue[[tmp_cluster]][tmp_GO_pvalue_list[[tmp_cluster]]$terms,])
  
  
  tmp_ancestors<-get_ancestors(GO_infor, terms=tmp_GO_list[[tmp_cluster]]$GO_ID)
  tmp_ancestors.df<-data.frame(matrix(data=NA,nrow=length(tmp_GO_list[[tmp_cluster]]$GO_ID),ncol=length(tmp_ancestors)))
  rownames(tmp_ancestors.df)<-tmp_GO_list[[tmp_cluster]]$GO_ID
  colnames(tmp_ancestors.df)<-tmp_ancestors
  for(tmp_go_terms in tmp_GO_list[[tmp_cluster]]$GO_ID){
    tmp_ancestors.df[tmp_go_terms,get_ancestors(GO_infor, terms=tmp_go_terms)]=1
  }
  tmp_GO_ancestors_list[[tmp_cluster]]<-tmp_ancestors.df
  
  sorted_GO_ID<-names(sort(colSums(tmp_ancestors.df,na.rm = T),decreasing = T))
  
  tmp_summary.df<-data.frame("GO_ID"=sorted_GO_ID,
                             "GO_counts"=sort(colSums(tmp_ancestors.df,na.rm = T),decreasing = T),
                             "GO_names"=GO_infor$name[sorted_GO_ID])
  tmp_GO_ancestors_summary_list[[tmp_cluster]]<-tmp_summary.df
  
  
  
  
  #######get ancestor terms of first level in GO tree
  first_level<-c("GO:0044848","GO:0044419","GO:0051703","GO:0065007","GO:0009987","GO:0098754","GO:0032502","GO:0040007","GO:0042592","GO:0002376","GO:0051179","GO:0040011","GO:0008152","GO:0032501","GO:0043473","GO:0000003","GO:0022414","GO:0050896","GO:0048511","GO:0023052","GO:0016032")
  first_level<-intersect(first_level,colnames(tmp_ancestors.df))
  plot.scaled_p.first.df<-data.frame(matrix(data=NA,nrow=0,ncol=3))
  colnames(plot.scaled_p.first.df)<-c("GO_id","Ancestors_GO_id","Ancestors_term")
  for(tmp_ancestor1 in first_level){
    tmp_summary<-tmp_ancestors.df[,tmp_ancestor1]
    names(tmp_summary)<-rownames(tmp_ancestors.df)
    tmp_summary<-na.omit(tmp_summary)
    
    tmp.data.frame<-data.frame("GO_id"=names(tmp_summary),
                               "Ancestors_GO_id"=rep(tmp_ancestor1,length(names(tmp_summary))),
                               "Ancestors_term"=rep(GO_infor$name[tmp_ancestor1],length(names(tmp_summary))))
    
    plot.scaled_p.first.df<-rbind(plot.scaled_p.first.df,tmp.data.frame)
    
    
  }
  
  plot.scaled_p.first.df[,colnames(tmp_GO_list[[tmp_cluster]][plot.scaled_p.first.df$GO_id,])]=tmp_GO_list[[tmp_cluster]][plot.scaled_p.first.df$GO_id,]
  
  
  plot.p.first.df<-plot.scaled_p.first.df[,1:3]
  plot.p.first.df<-cbind(plot.p.first.df,tmp_GO_pvalue_list[[tmp_cluster]][rownames(plot.p.first.df),])
  
  plot.p.first.df[,colnames(tmp_GO_pvalue_list[[tmp_cluster]][plot.p.first.df$GO_id,])]=tmp_GO_pvalue_list[[tmp_cluster]][plot.p.first.df$GO_id,]
  
  tmp_go_names<-names(table(plot.p.first.df$Ancestors_GO_id)[table(plot.p.first.df$Ancestors_GO_id)>4])
  p.list<-list()
  for(tmp_cluster1 in c("Cluster.1","Cluster.2","Cluster.3","Cluster.4","Cluster.5")){
    
    tmp_df1<-plot.p.first.df[plot.p.first.df$Ancestors_GO_id%in%tmp_go_names,]
    tmp_df1<-tmp_df1[,c("GO_id", "Ancestors_GO_id","Ancestors_term","GO_ID","GO_name","terms",tmp_cluster1)]
    colnames(tmp_df1)[7]<-"pvalue"
    plot.p.dot_plot<-ggplot(tmp_df1,aes(x=-log(pvalue) ,y=Ancestors_term,fill=Ancestors_term))+
      geom_jitter(aes(color=Ancestors_term),width = 0.1,height =0.1,stroke = 0.3)+
      stat_summary(color="black")+
      theme_bw()+
      theme(legend.position="none",text = element_text(size = 15))+
      xlab("-log(p value)")+
      ggtitle(paste0("DEG.",tmp_cluster1,"'s pvalue"))+ylab("Ancestors term(first level)")
    p.list[[tmp_cluster1]]=plot.p.dot_plot
    
    
  }
  pdf(paste0(prefix,".cluster.",tmp_cluster,".first_level.dotplot.pvalue_from.DEG_cluster.",".pdf"),width = 10,height = 4)
  print(p.list)
  dev.off()
  
  write.xlsx(plot.p.first.df,paste0(prefix,".cluster.",tmp_cluster,".first_level.dotplot.DEG_cluster.",".xlsx"))
  
  tmp_heatmap.list<-list()
  for(tmp_first_level_GO in c("GO:0002376")){
    
    tmp_heatmap_df<-plot.p.first.df[plot.p.first.df$Ancestors_GO_id%in%tmp_first_level_GO,]
    tmp_row_title=unique(tmp_heatmap_df$Ancestors_term)
    rownames(tmp_heatmap_df)<-tmp_heatmap_df$GO_name
    tmp_heatmap_df<-tmp_heatmap_df[,c("Cluster.1",    "Cluster.2","Cluster.3",    "Cluster.4",  "Cluster.5")]
    
    
    col_fun = colorRamp2(c(min(-log(tmp_heatmap_df)),1/4*max(-log(tmp_heatmap_df)),max(-log(tmp_heatmap_df))), c("#ffffff","#f1cbcd","#e00d17") )
    
    
    enrichment_heatmap<-Heatmap(-log(tmp_heatmap_df[order(rowMins(as.matrix(tmp_heatmap_df))),]),
                                name = "-log Pvalue",
                                show_column_names = T,
                                col = col_fun,
                                show_row_names = T,
                                column_names_rot = 90,
                                #left_annotation   = ha,
                                cluster_rows = F,
                                cluster_columns = F,
                                border=1,row_title=tmp_row_title)
    tmp_heatmap.list[[tmp_first_level_GO]]=enrichment_heatmap
  }
  
  pdf(paste0(prefix,".cluster.",tmp_cluster,".first_level.GO.pvalue.heatmap.pdf"),width = 6,height = 10)
  print(tmp_heatmap.list)
  dev.off()
  
  
  
}
