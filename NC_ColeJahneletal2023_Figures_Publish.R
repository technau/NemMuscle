### script for generating figures in NC-2023 Cole, Jahnel et al.
# generate datasets with NC_ColeJahneletal_DataProcessing.R first
# Setup ----
## load libraries 
library(easypackages)
libraries("Seurat", "Matrix", "readxl",# "dplyr", 
          "RColorBrewer", 'ggplot2','patchwork','pals')

## set color palettes
LibCP=c(brewer.paired(7),'black') #library
gene.cp=c('lightgrey',rev(brewer.pal(11 , "Spectral" ))) #gene rainbow
clust.cp.separate = unique (c(cols25(25),alphabet2(26),glasbey(32),alphabet(26))) #distinct color palette
muscle.cp=c("steelblue","steelblue3","#960096","darkgreen", "palegreen3",'tan3',
                        "black","grey60","grey60","grey60","grey80","grey80",
                        "grey80","grey80","grey80","grey80")#long palette for fine clusters
muscle.cp.coarse=c("steelblue3","#960096","darkgreen", "palegreen3",'tan3',"black","grey50","grey80","grey90") #coarse cluster palette
                        
## select whether or not you save output files: 
save.files = T
                        
## load genes / annotations: this is Sup.Data1.S1
load ('Genes.RData')

# Select parts to run:
Figure_1 = F
Figure_2 = F
Figure_3 = T
Figure_4 = F
Figure_5 = F

## load datasets ----
 # load the data objects
 load (file = '~/Robj/Tissues.Robj')
 load (file = '~/Robj/muscle.Robj')
 if (Figure_4 == T)
 {
 load (file = '~/RObj/nem64.Robj')
 load (file = '~/RObj/nem64.gastrodermis.RObj')
 load (file = '~/RObj/nem64.muscle7.Robj')
 }
 ### Set up the indices 
  #CHECK THAT IT IS WHAT YOU EXPECT :
  t=DimPlot(Tissues7,label = T, repel=T, cols=clust.cp.separate)+NoLegend()
  m=DimPlot(Muscle7,label = T,reduction = 'umap', cols=clust.cp.separate)+NoLegend()
  Muscle7 = SetIdent(Muscle7,value = 'IDs.fine')
  pie(table(Idents(Muscle7)), col = c(muscle.cp))
  t+m
  
  Muscle7 = SetIdent(Muscle7,value='IDs')
  
  # ID the muscle cells in the Tissue library:
  muscle.cells = Muscle7@assays$RNA@data@Dimnames[[2]]
  # ID the non muscle cells in the full dataset
  namesT5=levels(Tissues7)
  nonmuscle <- WhichCells(Tissues7,idents =namesT5[3:8])
  # ID the various muscle populations
  T_rm = WhichCells(Muscle7, idents = c('TR'))
  M_rm = WhichCells(Muscle7, idents = 'MR')
  P_em = WhichCells(Muscle7, idents = 'PM')
  R_em = WhichCells(Muscle7, idents = 'CM')
  allMuscle = c(T_rm,M_rm,P_em,R_em)

## load BULK data ----

source ('~/Robj/NC_ColeJahnel_etal_BULKdata.R')
ann.counts.filtered <- na.omit(ann.counts)
ann.mean.tpm.filtered <- na.omit(ann.mean.tpm)
DEBulk.filteredUP <- ann.counts.filtered[(ann.counts.filtered$logFC <= -1.5 & ann.counts.filtered$FDR < .001),]
DEBulk.filteredDOWN <- ann.counts.filtered[(ann.counts.filtered$logFC >= 1.5 & ann.counts.filtered$FDR < .001),]
#clean up the workspace:
rm(ann.mean.tpm.m, gene.info, gene.info.length, gene.info.length.de, nves, 
   cols, genes.to.plot,lmean.tpm, mean.tpm)  
# generate a list of TFs:
t<-nchar(DEBulk.filteredUP$deM_TF)
TFBulkNVE = DEBulk.filteredUP$GeneId[which(t>0)]
TFBulk <- genes$gene_short_name[match(TFBulkNVE,genes$NVE)]
  # names = levels(Tissues7)
  # levels(Tissues7@active.ident) <- c(1:length(names))

 
  # make the figures/data output ----
### Figure_1 ----
  if (Figure_1)
  {
    #### Figure 1C ----
    GOI = c(DEBulk.filteredDOWN$GeneId,DEBulk.filteredUP$GeneId)
    plot.these <- match(GOI,ann.mean.tpm$GeneId)
    cp <-colorRampPalette(c("white","light grey","dark red"))
    F1C=pheatmap::pheatmap(tpm[rev(plot.these),c(3:4,1:2)],cluster_cols = F, cluster_rows = T,scale = 'row',
                           labels_row = rev(GOI),
                           show_colnames = T,color = cp(50),show_rownames = F)
    source.fig1c=tpm[rev(plot.these),c(3:4,1:2)]
    rownames(source.fig1c)=rev(GOI)
    write.csv(source.fig1c,file='~/raw.source.data/source.fig1c.csv')
    # save list of upregulated genes from bulkdataset | extended data 1
    DEBulk.filtered = rbind(DEBulk.filteredUP,DEBulk.filteredDOWN)
    DEBulk.filtered$gene <- 'NA'
    DEBulk.filtered$go.annotation <- 'NA'
    
    ind = NULL
    for ( i in 1:length(DEBulk.filtered$GeneId))
      ind[i]=match(DEBulk.filtered$GeneId[i],annotations$NVE) 
    
    DEBulk.filtered$gene = annotations$gene_short_name[ind]
    DEBulk.filtered$go.annotation = annotations$gene_ontology_pfam[ind]
    write.csv(DEBulk.filtered,"~/raw.source.data/DataS2_BULK.csv")
    
    #### Figure 1D ----     
    # calculate cluster-specific genes for Tissues:  
    all.markers_Tissues <- FindAllMarkers(Tissues7,
                                          logfc.threshold = 1,
                                          only.pos = TRUE,
                                          return.thresh = 0.00001)# 
    
    # add GO terms and NVEs associated with this list:
    all.markers_Tissues$go.annotation <- 'NA'
    all.markers_Tissues$NVE <- 'NA'
    for (i in 1:length(levels(Tissues7@active.ident))) # 
    {
      x=all.markers_Tissues[as.numeric(all.markers_Tissues$cluster)==i,][1:length(which(as.numeric(all.markers_Tissues$cluster)==i)),7]
      anInd = match(genes[match(x,genes$gene_short_name),1],annotations$NVE)
      print(unique(annotations$gene_ontology_pfam[anInd]))
      all.markers_Tissues[as.numeric(all.markers_Tissues$cluster)==i,][1:length(which(as.numeric(all.markers_Tissues$cluster)==i)),8]<-annotations$gene_ontology_pfam[anInd]
      all.markers_Tissues[as.numeric(all.markers_Tissues$cluster)==i,][1:length(which(as.numeric(all.markers_Tissues$cluster)==i)),9]<-annotations$NVE[anInd]
    }  
    
    write.csv(all.markers_Tissues,file = "~/raw.source.data/DataS3_Tissues.csv")
    
    all.markers_variable=all.markers_Tissues
    list = NULL
    for (i in 1:length(levels(Tissues7@active.ident))) # 
    {
      x=all.markers_variable[as.numeric(all.markers_variable$cluster)==i,][1:min(10,length(which(as.numeric(all.markers_variable$cluster)==i))),7]
      if (is.na (x) ==F)
        list=c(list,x)  }
    list = unique(list)
    
    # generate population average dataset:
    AveTissues <- AverageExpression (Tissues7, return.seurat = T, use.scale = F)
    
    
    F1D=DoHeatmap(AveTissues, features = list, group.bar = T, draw.lines = F, label = F, size = 3, 
                  slot = 'scale.data')   + scale_fill_gradient2(low='white', mid="lightgrey", high="darkred")  +
      theme(legend.position = 'none',axis.text.y = element_blank())+
      labs(title = "Cell clustering and analysis",subtitle = "Single cell transcriptomes")
    source.fig1d=F1D$data
    write.csv(source.fig1d,file='~/raw.source.data/source.fig1d.csv')
    
    # theme(axis.text.x =
    layoutF1CDE=c(area(1,1),area(1,2),area(1,3))#,area(2,2))
    plot(layoutF1CDE)
    # generate plotting list for figure      
    all.markers_variable=all.markers_Tissues# 
    clusterLabels = NULL
    for (i in 1:length(levels(Tissues7@active.ident))) # 
    {
      x=all.markers_variable[as.numeric(all.markers_variable$cluster)==i,][1:min(3,length(which(as.numeric(all.markers_variable$cluster)==i))),7]
      if (is.na (x) ==F)
        clusterLabels=c(clusterLabels,x)  
    }
    clusterLabels = unique(clusterLabels)
    
    #### Figure 1D.2 ----
    Tissues7@active.assay='MAGIC_RNA'
    F1D.2=    FeaturePlot(Tissues7,features = c('Nve-MyHC-st'),
                          # min.cutoff = 0.5, #1.5 no inference
                          max.cutoff = 4, # this sets the limit for showing expression
                          cols = c("lightgrey","darkred"),# reduction = "umap",
                          order = T, label = T, repel = T, label.size = 4)& NoAxes()#&NoLegend()
    Tissues7@active.assay='RNA'
    
    #### Figure 1E.2 ----
    Muscle7@active.assay='MAGIC_RNA'
    F1E.2=
      FeaturePlot(Muscle7,features = c('Nve-MELC4','MYPH-like8','Nve-MELC3'),#c('CD151-like3','SEGN-like3','TBA1-like7','PRD14-like3'),#
                  # min.cutoff = 1, max.cutoff = 3, # this sets the limit for showing expression
                  cols = gene.cp,# c("lightgrey","darkred"), #reduction = "umap",
                  ncol=1, order = T, label = F, repel = T, label.size = 4)& NoAxes()&NoLegend()
    Muscle7@active.assay='RNA'
    
    #### Figure 1E ----
        F1E=DimPlot(Muscle7,reduction = 'umap',label = F, pt.size = 0.25,repel=T,
                    group.by = 'IDs.fine',
                cols = c('steelblue3',muscle.cp ))+ NoLegend()+NoAxes()+
      labs(title = 'Cell plot',subtitle = 'Muscle Subset')
    
    F1c=ggplot()+labs(title = "Differential gene expression",
                      subtitle = 'Bulk transcriptomes')+ theme_classic() #F1C is not a plot...
    
    ### Plot Figure 1C:E ----
    F1C
    F1c+F1D+F1E+plot_layout(ncol = 3)
    
    #### Figure S1 ----
    FS1g=VlnPlot(Tissues7,features = c('nCount_RNA','nFeature_RNA'),
                 group.by = 'orig.ident' ,cols = LibCP,pt.size = 0.1)
    source.figS1g=FS1g$data
    write.csv(source.figS1g,file='~/raw.source.data/source.figS1g.csv')
    

    #### Figure S2b ----   
    S2b=DotPlot(Tissues7, features = (clusterLabels),
                cols = c('grey','darkred')) + RotatedAxis ()+
      # NoLegend()+
      labs(title = 'Top 3 Cluster Marker Genes')
    source.figS2b=S2b$data
    write.csv(source.figS2b,file='~/raw.source.data/source.figS2b.csv')
    
    # Figure 1E | UMAP of MUSCLE dataset in clusters
    Muscle7<-SetIdent(Muscle7,value = 'IDs.fine')
    
    
    
    #### Figure S2c ----  
    #plot the Bulk genes on the full dataset:
    plot.ind <- match(DEBulk.filteredUP$GeneId, genes$NVE,nomatch = '0')
    plot.ind=plot.ind[plot.ind>=1]
    plot.these <- (genes$gene_short_name[plot.ind])
    # filter
    ind2 = NULL
    ind2 = which(rowSums(as.matrix(Tissues7@assays$RNA@counts[plot.these,T_rm]))>50) 
    plot.these = plot.these[ind2]
    S2c=DotPlot(Tissues7, features = rev(plot.these), dot.min = 0.01,col.min = 0,
                  scale.by = 'size',
                  cols =c("lightgrey","darkred")) +
      RotatedAxis()+theme(axis.text.x = element_blank())+NoLegend()+
      labs(title = 'Bulk DEG in scDataset')
    Tissues.CP = cols25(length(levels(Tissues7)))
    source.figS2c=S2c$data
    write.csv(source.figS2c,file='~/raw.source.data/source.figS2c.csv')
    
    #### Figure S2a ----
    S2a.1=DimPlot(Tissues7,label = F,repel = T,cols = Tissues.CP)+
      NoAxes()+labs(title='Tissues scDataset | Clusters')
    #### Figure S2b ----
    S2a.2=DimPlot(Tissues7,group.by = 'orig.ident',order = levels(Tissues7$orig.ident)[c(1,2,5,3,4)],
                  cols = LibCP)+NoAxes()+labs(title='Tissue Origin')
    
    layoutS2 = c(area(1,1,2,1),area(1,2,2,2),area(3,1,3,2),area(4,1,4,2))
    plot(layoutS2)
    
    ## plot Figure S2 ----
    S2a.1+S2a.2+S2b+S2c+plot_layout(design = layoutS2)
  }

  # put Muscle7 into Tissues7
  Muscle7@active.ident=Muscle7$IDs
  coi=Muscle7@assays$RNA@counts@Dimnames[[2]]
  Tissues7 <- SetIdent(Tissues7, cells = coi, value = Muscle7@active.ident[coi])
  levels(Tissues7)
  # re-order for plotting:
  Tissues7@active.ident = factor(Tissues7@active.ident ,
                                 levels(Tissues7@active.ident )[c(1:9,12,15,13,14,11,10)])
## Figure_2 ----
  if (Figure_2)
  {
  
    
    #supplement:
    #### Figure S3a ----
    FS3a = DimPlot(Muscle7, reduction = 'umap', group.by = 'orig.ident',cols=LibCP)&NoAxes()
    
    Muscle7$orig.ident = droplevels(as.factor(Muscle7$orig.ident))
    ids.cluster.library = as.data.frame(table(Idents(Muscle7), Muscle7@meta.data$orig.ident))
    colnames(ids.cluster.library) = c('ID','Library','CellCount')
    
    #### Figure S3b ----
    FS3b =  ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                            x=Library)) +
      geom_bar(mapping =aes(fill=ID, y= (CellCount),
                            x=(Library)),
               position="fill", stat="identity", width = 0.5)+
      scale_fill_manual(values = muscle.cp.coarse)+
      theme(axis.text.x = element_text(#face="bold", color="#993333", 
        size=8, angle=-45,hjust=0,vjust = 0.5))+
      geom_area(mapping =aes(fill=ID, y= (CellCount),
                             x=as.integer(Library)),
                position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
      geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                            x=(Library)),
               position="fill", stat="identity", width = 0.5)+
      ggtitle("Distribution of cell types in time and space")
    
    Muscle7@active.assay='MAGIC_RNA'
    #### Figure S3c ----
    
    FS3c = FeaturePlot(Muscle7,'Nve-MyHC-st',order=T,cols=gene.cp,min.cutoff = 0.2)&NoAxes()
    Muscle7@active.assay='RNA'
    
    
    FS3a+FS3b+FS3c
    #### Figure S3e ----
    #1) filter out the TFs:
    Muscle7<- SetIdent(Muscle7,value = 'IDs')
    Muscle7@active.assay='RNA'
    GOI = setdiff(Muscle7@assays$RNA@counts@Dimnames[[1]],TF_list$gene_short_name)
    #remove ribosomal:
    GOI = setdiff(GOI,ribosomal)
    #generate DEG list:
    all.markers_Muscle <- FindAllMarkers(Muscle7,
                                         features = GOI,
                                         logfc.threshold = 0.6,
                                         only.pos = T,
                                         min.pct = 0.1,
                                         return.thresh = 0.0001)# ,
    
    # add GO terms and NVEs associated with this list:
    all.markers_Muscle$go.annotation <- 'NA'
    all.markers_Muscle$NVE <- 'NA'
    for (i in 1:length(levels(Muscle7@active.ident))) # 
    {
      x=all.markers_Muscle[as.numeric(all.markers_Muscle$cluster)==i,][1:length(which(as.numeric(all.markers_Muscle$cluster)==i)),7]
      anInd = match(genes[match(x,genes$gene_short_name),1],annotations$NVE)
      print(unique(annotations$gene_ontology_pfam[anInd]))
      all.markers_Muscle[as.numeric(all.markers_Muscle$cluster)==i,][1:length(which(as.numeric(all.markers_Muscle$cluster)==i)),8]<-annotations$gene_ontology_pfam[anInd]
      all.markers_Muscle[as.numeric(all.markers_Muscle$cluster)==i,][1:length(which(as.numeric(all.markers_Muscle$cluster)==i)),9]<-annotations$NVE[anInd]
    }  
    # save the list of structural DEGs from the 'muscle set'
    if (save.files){
    write.csv(all.markers_Muscle,"~/raw.source.data/DataS4_MuscleDEGs.csv")
    }
    data1=Muscle7
    
    # generate plotting list for figure 
    all.markers=all.markers_Muscle
    list = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x=all.markers[as.numeric(all.markers$cluster)==i,][1:min(10,length(which(as.numeric(all.markers$cluster)==i))),7]
      if (is.na (x) ==F)
        list=c(list,x)
    }
    
    #Image the list
    FS3e=DotPlot(data1, features = unique(c(list)), 
                   scale.by='size' , col.min = 0, col.max = 3,  
                   cols = c('lightgrey','darkred')) + FontSize(6,8)+
      RotatedAxis()  +#coord_flip()+
      # NoLegend()+
      labs(title = 'Top 10 DEGs')
    FS3e
    source.figS3e=FS3e$data
    write.csv(source.figS3e,file='~/raw.source.data/source.figS3e.csv')
    
    all.markers_variable=all.markers_Muscle# 
    listM.1 = NULL
    for (i in 1:4)#6:11
    {
      x=all.markers_variable[as.numeric(all.markers_variable$cluster)==i,][1:min(10,length(which(as.numeric(all.markers_variable$cluster)==i))),7]
      if (is.na (x) ==F)
        listM.1=c(listM.1,x)  
    }
    listM.1 = unique(listM.1)
    
    #### Figure S3d ----
    #run GO analysis #note: this is not as is in the supplement; also the DEGs are slightly different
    #also need to update the GO terms in the genes file.
    FigS3d = F
    if (FigS3d)
    { 
      library(topGO)
      library(tidyr)
      #get the GO terms into TOPGO friendly format
      Unfold <- genes %>% 
        mutate(`GO:terms` = strsplit(as.character(`GO:terms`), ",")) %>% 
        unnest(`GO:terms`) 
      geneID2GO <- Unfold %>% split(x = .$`GO:terms`, f = .$gene_short_name)
      geneNames <- names(geneID2GO)
      
      #set dataset:
      all.markers=all.markers_Muscle
      
      #process data
      GOdata = NULL
      # this takes a really long time:
      for(i in 1:length(unique(all.markers$cluster)))
      {
        
        #identify the GOI from the DE list
        myInterestingGenes <- all.markers$gene[all.markers$cluster==unique(all.markers$cluster)[i]] #list of genes you want to perform GO enrichment for
        
        geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
        names(geneList) <- geneNames
        GOdata[[i]] <- new("topGOdata",ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
      }
      
      # filter and generate figures.
      p=NULL
      resultFis=NULL
      for(i in 1:length(unique(all.markers$cluster)))
      {
        #run the test...
        resultFis[[i]] <- runTest(GOdata[[i]], algorithm = "classic", statistic = "fisher") 
        pvalFis <- score(resultFis[[i]])
        #filter for only >0.05
        pvalFis = pvalFis[pvalFis>=0.05]
        allRes_intgenes<- GenTable(GOdata[[i]], pvalues = resultFis[[i]], orderBy = "pvalues", topNodes=20)
        allRes_intgenes$pvalues<-as.numeric(allRes_intgenes$pvalues)
        #convert NA values to zero; only if p value is so small it cannot be displayed in r
        allRes_intgenes[is.na(allRes_intgenes)]<-0.00000001
        
        #plot GOenrichment 
        
        david_col=c("darkred")
        p[[i]]=
          ggplot2::ggplot(allRes_intgenes, ggplot2::aes(x=(reorder(Term,(-log10(pvalues)))), y=(-log10(pvalues)))) +
          ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
                                col=david_col,fill=david_col) +
          ggplot2::coord_flip()+
          ggplot2::xlab("Biological Process") +
          ggplot2::ylab("Enrichment -log10 p-value") +
          ggplot2::labs(title="GO enrichment",subtitle=levels(all.markers$cluster)[i])
        
        # p2[[i]]=
        #   ggplot2::ggplot(allRes_intgenes, ggplot2::aes(x=fct_inorder(Term), y=(-log10(pvalues)))) +
        #   ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
        #                         col=david_col,fill=david_col) +
        #   ggplot2::coord_flip()+
        #   ggplot2::xlab("Biological Process") +
        #   ggplot2::ylab("Enrichment -log10 p-value") +
        #   ggplot2::labs(title="GO enrichment",subtitle=levels(all.markers$cluster)[i])
      }
      
      #print figures
      p[[1]]+p[[2]]+p[[3]]+p[[4]]+p[[5]]+
        p[[6]]+p[[7]]+p[[8]]+
        p[[9]]
    }
    
    layoutF2 = c(area(1,1,4,1),area(1,3,4,3),area(4,2,4,2))
    plot(layoutF2)
    #### Figure 2a ----
    
    Fig2A.raw=DotPlot(Muscle7, features = rev(listM.1), 
                      idents = levels(Muscle7),
                      dot.min = 0.01,col.min = 0,
                      scale.by = 'size',cols =c("lightgrey","darkred")) +RotatedAxis()+coord_flip()
    source.fig2a=Fig2A.raw$data
    write.csv(source.fig2a,file='~/raw.source.data/source.fig2a.csv')
    
    NVEGenes = genes$NVE [match(rev(listM.1),genes$gene_short_name)]# filterGOIs
    
    plot.these <- match(NVEGenes,ann.mean.tpm$GeneId)
    cp <-colorRampPalette(c("white","lightgrey","darkred"))
    F2aBulkRaw=pheatmap::pheatmap(tpm[rev(plot.these),],cluster_cols = F, cluster_rows = F,scale = 'row',
                                  labels_row = rev(NVEGenes),
                                  show_colnames = T,color = cp(50))
    
    
    #### Figure 2c ----
    # generate Plots of paralogs:
    StructuralParalogs<- read_excel("D:/Dropbox/Nematostella/documents/Genelists/NVeannotations_NoDash.xls",
                                    sheet = 'musclegeneparalogs')
    Sabrina =StructuralParalogs
    # indexing NVE list
    ind = NULL
    for (i in 1:length(Sabrina$NVE))
      ind[i] = match(Sabrina$NVE[i],genes$NVE,nomatch = 0)
    ind = unique(ind[ind>0])
    
    #  # filter step: 
    ind2 = NULL
    ind2 = which(rowSums(as.matrix(Muscle7@assays$RNA@counts[ind,c(allMuscle)]))>250)
    indExp = ind[ind2]
    p=genes[indExp,]
    Paralogs = unique( c('Nve-MyHC-nm',genes$gene_short_name[indExp]))
    
    Fig2C.raw=DotPlot(Muscle7, features = (Paralogs), 
                      # idents = levels(Muscle7)[2:5], #better without because otherwise the CM is really low in the muscle genes compared to the rest!
                      col.min = 0, col.max = 1, dot.min = 0.01,
                      scale.by = 'size',cols =c("lightgrey","black"))+RotatedAxis() + 
      FontSize(8)+coord_flip()
    source.fig2c=Fig2C.raw$data
    write.csv(source.fig2c,file='~/raw.source.data/source.fig2c.csv')
    
    ### Plot Figure 2
    Fig2A.raw+Fig2C.raw
    NVEGenes = unique( c('NVE4823',genes$NVE[indExp]))
    GenesNames = unique( c('Nve-MyHC-nm',genes$gene_short_name[indExp]))
    plot.these <- match(NVEGenes,ann.mean.tpm$GeneId)
    plot.these <- rev(plot.these)
    cp <-colorRampPalette(c('white',"light grey", "black"))
    #Fig2C.bulk
    F2cBulkRaw=pheatmap::pheatmap(tpm[plot.these,],cluster_cols = F, cluster_rows = F,scale = 'row',labels_row = rev(NVEGenes),
                                  show_colnames = T,color = cp(50), fontsize = 8)# +scale_colour_gradient(low ='light blue',high ='purple') # this is more promissing...
    #export as 2x5.5
    F2aBulkRaw
    F2cBulkRaw
    
    #### Figure 2d ----
    # speed plot:
    speed = read_excel("~/SupplementaryData2.xlsx")
    speed$Group = as.factor((speed$Group))
    speed$log = log10(speed$Speed)
    speed$colour = as.factor(as.numeric(speed$Group))
    

    library(ggplot2)
    Fig2d=ggplot(speed, aes(x=Group, y=Speed,colour=colour,notch=F)) + 
      geom_boxplot(fill = c("steelblue3", "darkred", "darkgreen")) + 
      scale_color_manual(values = c("black", "black", "black"))+ theme_classic()+
      theme(legend.position="none")+
      scale_x_discrete(limits=c("Tentacle.Retraction", "Body.Retraction", "Body.Parastalsis"))+
      scale_y_continuous(trans='log10')
    source.fig2d=Fig2d$data
    write.csv(source.fig2d,file='~/raw.source.data/source.fig2d.csv')
    
    
    Fig2A.raw+NoLegend()+Fig2C.raw+NoLegend()+Fig2d+plot_layout(design=layoutF2)
    #### Figure S4 ----
    #supplemental:
    
    # generate Plots of in situ genes:
    inSituPlot<- read_excel("D:/Dropbox/Nematostella/documents/Genelists/NVeannotations_NoDash.xls",
                            sheet = 'insitus')
    
    FigS4=DotPlot(Tissues7, features = c(inSituPlot$gene_short_name), 
                    col.min = 0, col.max = 1, #dot.min = 0.01,
                    scale.by = 'size',
                    cols =c("lightgrey","darkred"))+
      RotatedAxis() + FontSize(8)
    
    source.FigS4=FigS4$data
    write.csv(source.FigS4,file='~/raw.source.data/source.FigS4.csv')
    
    # lists of 'structural genes' filtered for muscle:
    #   detectable expression
    data1 = Muscle7
    # filter the list according to lack of expression outside muscle cells:
    f=20
    M_RM_ST = genes$gene_short_name[which(rowSums(as.matrix(data1@assays$RNA@counts[genes$gene_short_name,M_rm]))>(f))]
    M_RM_ST = setdiff(M_RM_ST,TF_list$gene_short_name)# remove TFs
    M_RM_ST = setdiff(M_RM_ST,ribosomal)# remove ribosomal proteins
    M_RM_ST <- M_RM_ST[which(rowSums(as.matrix(Tissues7@assays$RNA@counts[M_RM_ST,nonmuscle]))<200)]# filter for specificity
    T_rm_ST = genes$gene_short_name[which(rowSums(as.matrix(data1@assays$RNA@counts[genes$gene_short_name,T_rm]))>f)]
    T_rm_ST = setdiff(T_rm_ST,TF_list$gene_short_name)# remove TFs
    T_rm_ST = setdiff(T_rm_ST,ribosomal)# remove ribosomal proteins
    T_rm_ST <- T_rm_ST[which(rowSums(as.matrix(Tissues7@assays$RNA@counts[T_rm_ST,nonmuscle]))<200)]# filter for specificity
    
    P_em_ST = genes$gene_short_name[which(rowSums(as.matrix(data1@assays$RNA@counts[genes$gene_short_name,P_em]))>f)]
    P_em_ST = setdiff(P_em_ST,TF_list$gene_short_name)# remove TFs
    P_em_ST = setdiff(P_em_ST,ribosomal)# remove ribosomal proteins
    P_em_ST <- P_em_ST[which(rowSums(as.matrix(Tissues7@assays$RNA@counts[P_em_ST,nonmuscle]))<200)]# filter for specificity
    
    R_em_ST = genes$gene_short_name[which(rowSums(as.matrix(data1@assays$RNA@counts[genes$gene_short_name,R_em]))>f)]
    R_em_ST = setdiff(R_em_ST,TF_list$gene_short_name)# remove TFs
    R_em_ST = setdiff(R_em_ST,ribosomal)# remove ribosomal proteins
    R_em_ST <- R_em_ST[which(rowSums(as.matrix(Tissues7@assays$RNA@counts[R_em_ST,nonmuscle]))<200)]# filter for specificity
    
    # set up a Venn Diagram:
    library (VennDiagram)
    a1=(M_RM_ST)
    a2=(T_rm_ST)
    a3=(P_em_ST)
    a4=(R_em_ST)
    plot.new()
    # FigS2.2 center
    ST.venn.plot <- draw.quad.venn(length(a1),length(a2),length(a3),length(a4),
                                   length(intersect(a1,a2)),length(intersect(a1,a3)),
                                   length(intersect(a1,a4)),length(intersect(a2,a3)), 
                                   length(intersect(a2,a4)), length(intersect(a3,a4)), 
                                   length(intersect(a1,(intersect(a2,a3)))),
                                   length(intersect(a1,(intersect(a2,a4)))),length(intersect(a1,(intersect(a3,a4)))),
                                   length(intersect(a2,(intersect(a3,a4)))),
                                   length(intersect(a1,(intersect(a2,intersect(a3,a4))))),
                                   category = c('Mesentery','Tentacle','Parietal','Ring'),
                                   fill =c("#960096","steelblue3","darkgreen","palegreen3"),
                                   cex=2)
    
    StructuralGeneSet <- unique(c(a1,a2,a3,a4))
    StructuralGeneSet <- genes[match(StructuralGeneSet, genes$gene_short_name),]
    write.csv(StructuralGeneSet,"D:/Alison/manuscripts/MusclePaper//NatureCom.2023/accepted/raw.source.data/DataS5_MuscleAllStructural.csv")
    
    
    FastCommon <-  setdiff(setdiff(intersect(a1,a2),intersect(a3,a4)),a3)
    FastCommon <- genes[match(FastCommon, genes$gene_short_name),]
    
    
    FS5fast=DotPlot(Tissues7, features = rev(FastCommon$gene_short_name), dot.min = 0.01,col.min = 0,
                      scale.by = 'size',cols =c("lightgrey","darkred")) +
      RotatedAxis()+FontSize(5,)
    source.FS5fast=FS5fast$data
    write.csv(source.FS5fast,file='~/raw.source.data/source.FS5fast.csv')
    
    SlowCommon <- setdiff(setdiff(setdiff(intersect(a3,a4),intersect(a1,a2)),a1),a2)
    FS5slow=DotPlot(Tissues7, features = c(SlowCommon), dot.min = 0.01,col.min = 0,
                      scale.by = 'size',cols =c("lightgrey","darkred")) +RotatedAxis()
    source.FS5slow=FS5slow$data
    write.csv(source.FS5slow,file='~/raw.source.data/source.FS5slow.csv')
    
    # AllMuscle <-  intersect(a1,(intersect(a2,intersect(a3,a4))))
    # FS5all=DotPlot(Tissues7, features = c(AllMuscle), dot.min = 0.01,col.min = 0,
    #                  scale.by = 'size',cols =c("lightgrey","darkred")) +
    #   RotatedAxis()+NoLegend()
    # source.FS5all=FS5all$data
    # write.csv(source.FS5all,file='~/raw.source.data/source.FS5all.csv')
    # 
    layout.FigS5=c(area(1,1,1,4),area(2,1),area(2,2,2,4))
    plot(layout.FigS5)
    FS5slow+ggplot()+FS5fast+plot_layout(design = layout.FigS5)
    
  }
  
## Figure_3 ----
  if (Figure_3)
  {
    SelectedTFs<- c("NVE3960","NVE19510","NVE17916","NVE24920","NVE23578","NVE19720","NVE8199",
    "NVE7115","NVE10444","NVE1324","NVE24655","NVE8968","NVE19721","NVE18116","NVE18117",
    "NVE2043","NVE1803","NVE4006","NVE10560","NVE4528","NVE9897","NVE24170","NVE4363","NVE4362","NVE8569","NVE8570")
    SelectedTFs = genes$gene_short_name[match(SelectedTFs,genes$NVE)]
      # read_excel('D:/Alison/manuscripts/MusclePaper/Archive/AnalysisSCRNAseq/SelectedTFs.xlsx')
    
    Fig3A=DotPlot(Tissues7, features = (SelectedTFs), col.min = 0, col.max = 1, 
                  scale.by = 'size',cols =c("lightgrey","darkred"),dot.min = 0.01)+RotatedAxis()+FontSize(8)# 
    Fig3A#+coord_flip()
    source.fig3a=Fig3A$data
    write.csv(source.fig3a,file='~/raw.source.data/source.fig3a.csv')
    
    #supplementary materials:

    # DEG algorythm
    all.markers_TF_Muscle <- FindAllMarkers(Muscle7, 
                                            features = intersect(TF_list$gene_short_name,rownames(Muscle7)),
                                            logfc.threshold = 0.1,
                                            only.pos = T,
                                            min.pct = 0.01,
                                            return.thresh = 0.0001)# 
    # add GO terms associated with this list:
    all.markers_TF_Muscle$NVE <- 'NA'
    all.markers_TF_Muscle$go.annotation <- 'NA'
    
    for (i in 1:length(levels(Muscle7@active.ident))) # 
    {
      x=all.markers_TF_Muscle[as.numeric(all.markers_TF_Muscle$cluster)==i,][1:length(which(as.numeric(all.markers_TF_Muscle$cluster)==i)),7]
      anInd = match(genes[match(x,genes$gene_short_name),1],annotations$NVE)
      print(unique(annotations$gene_ontology_pfam[anInd]))
      all.markers_TF_Muscle[as.numeric(all.markers_TF_Muscle$cluster)==i,][1:length(which(as.numeric(all.markers_TF_Muscle$cluster)==i)),8]<-annotations$NVE[anInd]
      all.markers_TF_Muscle[as.numeric(all.markers_TF_Muscle$cluster)==i,][1:length(which(as.numeric(all.markers_TF_Muscle$cluster)==i)),9]<-annotations$gene_ontology_pfam[anInd]
    }  
    
    # save the list of DETFs from the 'muscle set'
    if (save.files){
    write.csv(all.markers_TF_Muscle,"~/raw.source.data/DataS6_MuscleDETFs.csv")
    }
    all.markers_variable=all.markers_TF_Muscle
    list_TF = NULL
    for (i in 1:9)#(length(levels(Muscle7@active.ident))))
    {
      x=all.markers_variable[as.numeric(all.markers_variable$cluster)==i,][1:min(10,length(which(as.numeric(all.markers_variable$cluster)==i))),7]
      if (is.na (x) ==F)
        list_TF=c(list_TF,x)
    }
    list_TF = unique(list_TF)
    
    FS3f=DotPlot(Muscle7, features = (list_TF), col.min = 0, col.max = 3, 
                   scale.by = 'size',cols =c("lightgrey","darkred"))+
      RotatedAxis()+FontSize(6,8)+
      labs(title = 'Top 10 DETFs')# 
    source.FS3f=FS3f$data
    write.csv(source.FS3f,file='~/raw.source.data/source.FS3f.csv')
    
    FS8a=DotPlot(Tissues7, features = (list_TF), col.min = 0, col.max = 3, 
                 scale.by = 'size',cols =c("lightgrey","darkred"))+
      RotatedAxis()+FontSize(6,8)+
      labs(title = 'Top 10 DETFs')# 
    source.FS8a=FS8a$data
    write.csv(source.FS8a,file='~/raw.source.data/source.FS8a.csv')
    
    if (Figure_2 && Figure_3){FS3e+FS3f+plot_layout(ncol=1)}
    
    TFOI = TF_list
    ind=match(genes$gene_short_name,TFOI$gene_short_name,nomatch = 0)
    TFOI = TFOI[ind[ind>0],]
    data1 = Muscle7
    # must have detection above f reads to count...
    f=5
    M_RM_TF = TFOI$gene_short_name[which(rowSums(as.matrix(data1@assays$RNA@counts[TFOI$gene_short_name,M_rm]))>=f)]
    M_RM_TFf <- M_RM_TF[which(rowSums(as.matrix(Tissues7@assays$RNA@data[M_RM_TF,nonmuscle]))<50)]# 75 filter for specificity
    T_rm_TF = TFOI$gene_short_name[which(rowSums(as.matrix(Muscle7@assays$RNA@counts[TFOI$gene_short_name,T_rm]))>=f)]
    T_rm_TFf <- T_rm_TF[which(rowSums(as.matrix(Tissues7@assays$RNA@data[T_rm_TF,nonmuscle]))<50)]# filter for specificity
    
    P_em_TF = TFOI$gene_short_name[which(rowSums(as.matrix(data1@assays$RNA@counts[TFOI$gene_short_name,P_em]))>=f)]
    P_em_TFf <- P_em_TF[which(rowSums(as.matrix(Tissues7@assays$RNA@data[P_em_TF,nonmuscle]))<50)]# filter for specificity
    
    R_em_TF = TFOI$gene_short_name[which(rowSums(as.matrix(data1@assays$RNA@counts[TFOI$gene_short_name,R_em]))>=f)]
    R_em_TFf <- R_em_TF[which(rowSums(as.matrix(Tissues7@assays$RNA@data[R_em_TF,nonmuscle]))<50)]# filter for specificity
    
    # set up a Venn Diagram:
    
    a1=(M_RM_TF)
    a2=(T_rm_TF)
    a3=(P_em_TF)
    a4=(R_em_TF)
    a1.f=(M_RM_TFf)
    a2.f=(T_rm_TFf)
    a3.f=(P_em_TFf)
    a4.f=(R_em_TFf)
    dev.off()
    #FigS4.1B
    TF.venn.plot <- VennDiagram::draw.quad.venn(length(a1),length(a2),length(a3),length(a4),
                                                length(intersect(a1,a2)),length(intersect(a1,a3)),
                                                length(intersect(a1,a4)),length(intersect(a2,a3)), 
                                                length(intersect(a2,a4)), length(intersect(a3,a4)), 
                                                length(intersect(a1,(intersect(a2,a3)))),
                                                length(intersect(a1,(intersect(a2,a4)))),length(intersect(a1,(intersect(a3,a4)))),
                                                length(intersect(a2,(intersect(a3,a4)))),
                                                length(intersect(a1,(intersect(a2,intersect(a3,a4))))),
                                                category = c('Mesentery','Tentacle','Parietal','Ring'), 
                                                fill =c("#960096","steelblue3","dark green","palegreen3"),,
                                                cex=2)
    list.S8.2C=c(a2.f,a1.f,a3.f,a4.f)
    FS8C=DotPlot(Tissues7, features = unique(rev(list.S8.2C)), col.min = 0, col.max = 1, # dot.min = 0.2,
                   scale.by = 'size',cols =c("lightgrey","darkred"),dot.min = 0.0)+RotatedAxis()+FontSize(8)# 
    source.FS8C=FS8C$data
    write.csv(source.FS8C,file='~/raw.source.data/source.FS8C.csv')
    
    Muscl.TFoi <- genes[match(list.S8.2C, genes$gene_short_name),]
    
    write.csv(Muscl.TFoi,"~/raw.source.data/DataS7_MuscleAllTF.csv")
    
    # candidate gene family plotting   #FigS4.1D
    TF_Classes = c("zf-C2H2","T-box","bZIP","HLH","Fork_head","Homeobox","HMG_box",
                   "GATA","Ets","SRF","zf-C4") # this is in TF_list column 6
    
    HLH = TF_list[which(TF_list$TFclass == 4),2]
    SRF = TF_list[which(TF_list$TFclass == 10),2]
    GATA = TF_list[which(TF_list$TFclass == 8),2]
    
    FKH = TF_list[which(TF_list$TFclass == 5),2]
    TBX = TF_list[which(TF_list$TFclass == 2),2]
    HMG = TF_list[which(TF_list$TFclass == 7),2]
    Nk <- c('NVE10217','NVE10555','NVE10557','NVE10559','NVE10560','NVE11289','NVE11830','NVE18161',
            'NVE20899','NVE2477','NVE2478','NVE2480','NVE2485','NVE24919','NVE24920',
            'NVE6215','NVE7180')
    Nk <- genes[match(Nk,genes$NVE),2]
    
    SOX <- c('NVE13400','NVE15777','NVE16845','NVE17897','NVE2102','NVE426',
             'NVE23841','NVE3409','NVE20190','NVE24655','NVE23709','NVE12762','estExt_fgenesh1_pm.C_1290009',
             'NVE6454')
    SOX<- genes[match(SOX,genes$NVE),2]
    Homeobox = setdiff(TF_list[which(TF_list[,7] == 6),2],Nk)
    
    GOIs = rev(c(SRF,'Nve-MyoC',Nk,FKH,GATA,SOX,TBX, Homeobox, HLH))
    
    GOIs = intersect(GOIs,genes$gene_short_name)
    GOIsf <- (GOIs[which(rowSums(as.matrix(Muscle7@assays$RNA@counts[GOIs,c(T_rm,M_rm,P_em,R_em)]))>=5)]) 
    
    FS8D=DotPlot(Tissues7, features = unique(GOIsf), col.min = 0, col.max = 1, # dot.min = 0.2,
                   scale.by = 'size',cols =c("lightgrey","darkred"),dot.min = 0.0)+RotatedAxis()+FontSize(6)# 
    source.FS8D=FS8D$data
    write.csv(source.FS8D,file='~/raw.source.data/source.FS8D.csv')
    
    layout.FS8 =c(area(1,1,1,3),area(2,2,2,3),area(3,1,3,3))
    plot(layout.FS8)
    FS8a+FS8C+FS8D+plot_layout(design = layout.FS8)
  }
  # 
## Figure_4 ----
  if (Figure_4)
  {
    Tissues7@active.assay='MAGIC_RNA'
    Fig4d.1=FeaturePlot(Tissues7,features = c('Nve-MELC4'),
                cols = gene.cp,
                order = T, label = F, repel = T, label.size = 4)& 
      NoAxes()&NoLegend()
    Tissues7@active.assay='RNA'
    
    Fig4d.2=FeaturePlot(nem64,features = c('Nve-MELC4'),
                        cols = gene.cp,
                        order = T, label = F, repel = T, label.size = 4)& 
      NoAxes()&NoLegend()

    F4e=DimPlot(nem64.Muscle7,cells.highlight = dimnames(nem64.endo)[[2]],
                cols.highlight = 'black',
                label = F,repel = T)&NoAxes()
    
    FS4c=VlnPlot(nem64.Muscle7,c('NvSoxB.2','NvNem64','CALM-like3','Nve-MELC4','NVE8057','NvNem24','DMBX1-like'),
                 idents = levels(data1)[1:2],cols=c(muscle.cp),
                 group.by = 'orig.ident',
                 split.by = 'IDs',stack=T)
    
    #combine comparable cells: nem was only the head 
    nem=WhichCells(SetIdent(nem64.Muscle7,value='orig.ident'),idents = 'nem64')
    head=WhichCells(SetIdent(nem64.Muscle7,value='orig.ident'),idents = c('pharynx','tentacle','pharynx2'))
    tissue=setdiff(dimnames(nem64.Muscle7@assays$RNA@counts)[[2]],c(nem,head))
    nem64.Muscle7$genotype = NA
    nem64.Muscle7$genotype[nem]='nem64mutant.head'
    nem64.Muscle7$genotype[head]='wildtype.head'
    nem64.Muscle7$genotype[tissue]='wildtype.body'
    nem64.Muscle7$genotype = as.factor(nem64.Muscle7$genotype)
    nem64.Muscle7$genotype = factor(nem64.Muscle7$genotype,levels(nem64.Muscle7$genotype)[c(2,3,1)])
     clust.cp=c(muscle.cp.coarse)
    {
      
      ids.cluster.library = as.data.frame(table(Idents(nem64.Muscle7), nem64.Muscle7@meta.data$genotype))
      colnames(ids.cluster.library) = c('ID','Genotype','CellCount')
      
      genotype.plot = DimPlot(nem64.Muscle7,group.by = 'genotype',
                              order = rev(levels(nem64.Muscle7$genotype)),
                              cols = c('grey50','grey80','orange'))+NoAxes()+DarkTheme()
      labs(title = 'Time | Genotype origin')#+NoLegend()
      
      mutant.plot =  DimPlot(nem64.Muscle7,cells.highlight = dimnames(nem64.endo)[[2]],
                             cols.highlight = 'black',
                             label = F,repel = T)&NoAxes()&NoLegend()
      
      cluster.plot =DimPlot(nem64.Muscle7, label = F,label.size = 4, 
                            repel = T,order=rev(levels(nem64.Muscle7@active.ident)),#ordering is optional
                            cols = clust.cp,split.by = 'orig.ident')+NoAxes()+
        labs(title = 'Clusters | ID')+NoLegend()
      
      dist.clust2=
        ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                        x=Genotype)) +
        geom_bar(mapping =aes(fill=ID, y= (CellCount),
                              x=(Genotype)),
                 position="fill", stat="identity", width = 0.5)+
        scale_fill_manual(values = c(muscle.cp.coarse,muscle.cp.coarse[2:9]))+
        theme(axis.text.x = element_text(#face="bold", color="#993333", 
          size=8, angle=-45,hjust=0,vjust = 0.5))+
        geom_area(mapping =aes(fill=ID, y= (CellCount),
                               x=as.integer(Genotype)),
                  position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
        geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                              x=(Genotype)),
                 position="fill", stat="identity", width = 0.5)+
        ggtitle("Distribution of cell types in time and space")
      
     }
     Fig4f = dist.clust2
     
     Fig4d.1+F4e+Fig4d.2+Fig4f+plot_layout(design = c(area(1,1),area(1,2),area(2,1),area(2,2)))
  } 


## Figure_5 ----
  if (Figure_5)
{ 
  data1=Tissues7 # make sure the Muscle subset IDs were imported. ln 254 after Figure 1
  data1[["RNA"]]@meta.features <- data.frame(row.names = rownames(data1[["RNA"]]))

  #first generate list of all differentially expressed TFs 
  all.markers_TF <- FindAllMarkers(data1, only.pos = F, min.pct = 0.001,logfc.threshold = 0.4,
                                   features = intersect(TF_list$gene_short_name,rownames(data1)),
                                   return.thresh = 0.00001)#
  TFOI=unique(all.markers_TF$gene)

  #scale these; data scaling corrects for background noise from popped gastrodermal cells
  data1=ScaleData(data1,split.by = 'orig.ident',features = TFOI)
# tree reconstruction TFs:
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = F,
                            features = TFOI,slot='scale.data')
  
  ape::plot.phylo(data1@tools$BuildClusterTree)+title('Variable TF subset')

# tree reconstruction DE expressed genes without TFs
  all.markers <- FindAllMarkers(data1, only.pos = F, min.pct = 0.001,logfc.threshold = 1,
                                   return.thresh = 0.00001)#
  goi=unique(all.markers$gene)
  
  goi=setdiff(goi,TF_list$gene_short_name)
  data1=ScaleData(data1,split.by = 'orig.ident',features = goi)
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = F,slot='scale.data',
                            features = goi)
  ape::plot.phylo(data1@tools$BuildClusterTree,direction = 'right')+title('Variable Structural subset')
}