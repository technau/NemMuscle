### script for importing and analyzing data in NC-2023 Cole, Jahnel et al.
# Setup ----
## load libraries 
library(easypackages)
libraries("Seurat",  "readxl", #"Matrix","dplyr",
          "RColorBrewer", 'ggplot2','patchwork','pals')

## set color palettes
LibCP=c(brewer.paired(7),'black') #library
gene.cp=c('lightgrey',rev(brewer.pal(11 , "Spectral" ))) #gene rainbow
clust.cp.separate = unique (c(cols25(25),alphabet2(26),glasbey(32),alphabet(26))) #distinct color palette
muscle.cp=c("steelblue3","#960096","darkgreen", "palegreen3",'tan3',
                        "black","grey60","grey60","grey60","grey80","grey80",
                        "grey80","grey80","grey80","grey80")#long palette for fine clusters
muscle.cp.coarse=c("steelblue3","#960096","darkgreen", "palegreen3",'tan3',"black","grey50","grey80","grey90") #coarse cluster palette

## select whether or not you save output files: 
save.files = T

## load genes / annotations: this is Sup.Data1.S1
load ('Genes.RData')

# generate the files and analysis ----
  # Select parts to run:
  LoadData_raw = F
  Tissues7_analysis = T
  Muscle7_analysis = T
  Mutant = T
  
  if (LoadData_raw)
  {
    #  import the data matrices from the 10x Cellranger pipeline 'out' directory
    raw.datat <- Read10X(data.dir = "D:/Dropbox/Nematostella/data/Import4R/tissues2/tentacles")
    raw.datam <- Read10X(data.dir = "D:/Dropbox/Nematostella/data/Import4R/tissues2/mesentery")
    raw.datap <- Read10X(data.dir = "D:/Dropbox/Nematostella/data/Import4R/tissues2/pharynx")
    raw.datab <- Read10X(data.dir = "D:/Dropbox/Nematostella/data/Import4R/tissues2/bodywall")
    raw.data2 <- Read10X(data.dir = 'Z:/sequencing/Alison/Nematostella/NVE_mapping/TwistSequencing/TwistMutantTissue_12000x2')
    raw.datam2 <- Read10X(data.dir = 'Z:/sequencing/Alison/Nematostella/NVE_mapping/MesenteryFemale')
    raw.dataipBW<- Read10X(data.dir = 'Z:/sequencing/Alison/Nematostella/NVE_mapping/BW.interparietal')
    
    #  set gene model names to the annotations
    rownames(raw.datat) <- genes[, 2]
    rownames(raw.datam) <- genes[, 2]
    rownames(raw.datap) <- genes[, 2]
    rownames(raw.datab) <- genes[, 2]
    rownames(raw.data2) <- genes[, 2]
    rownames(raw.datam2) <- genes[, 2]
    rownames(raw.dataipBW) <- genes[, 2]
    
    #  generate individual Seurat objects:
    tentacle  <- CreateSeuratObject(counts = raw.datat, project = "tentacle")
    mesentery  <- CreateSeuratObject(counts = raw.datam, project = "mesentery")
    pharynx  <- CreateSeuratObject(counts = raw.datap, project = "pharynx")
    bodywall  <- CreateSeuratObject(counts = raw.datab, project = "bodywall")
    pharynx2  <- CreateSeuratObject(counts = raw.data2, project = "pharynx2")
    mesenteryF2  <- CreateSeuratObject(counts = raw.datam2, project = "mesenteryF")
    BW.ip  <- CreateSeuratObject(counts = raw.dataipBW, project = "BW.IP")    
    
    #  filter each library for outliers
    VlnPlot(tentacle, features = c('nFeature_RNA','nCount_RNA')) #Visual guide
    tentacle <- subset(x = tentacle, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)
    VlnPlot(mesentery, features = c('nFeature_RNA','nCount_RNA')) #Visual guide
    mesentery <- subset(x = mesentery, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)
    VlnPlot(pharynx, features = c('nFeature_RNA','nCount_RNA')) #Visual guide
    pharynx <- subset(x = pharynx, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)
    VlnPlot(bodywall, features = c('nFeature_RNA','nCount_RNA')) #Visual guide
    bodywall <- subset(x = bodywall, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)
    VlnPlot(pharynx2, features = c('nFeature_RNA','nCount_RNA')) #Visual guide
    pharynx2 <- subset(x = pharynx2, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)    
    VlnPlot(mesenteryF2, features = c('nFeature_RNA','nCount_RNA')) #Visual guide
    mesenteryF2 <- subset(x = mesenteryF2, subset = nFeature_RNA > 300 & nFeature_RNA < 1500)
    VlnPlot(BW.ip, features = c('nFeature_RNA','nCount_RNA')) #Visual guide
    BW.ip <- subset(x = BW.ip, subset = nFeature_RNA > 200)
    
    #  clean up the workspace
    rm (raw.datam, raw.datab,raw.datap, raw.datat, raw.data2, raw.datam2,raw.dataipBW)
    
    #  add library identifier to each sample in each library, and scale all genes
    tentacle <- RenameCells(tentacle, add.cell.id = "tentacle")
    mesentery <- RenameCells(mesentery, add.cell.id = "mesentery")
    mesenteryF2 <- RenameCells(mesenteryF2, add.cell.id = "F-mesentery2")
    pharynx <- RenameCells(pharynx, add.cell.id = "pharynx")
    bodywall <- RenameCells(bodywall, add.cell.id = "bodywall")
    pharynx2 <- RenameCells(pharynx2, add.cell.id = "pharynx2")
    BW.ip <- RenameCells(BW.ip, add.cell.id = "BWip")
    
    # generate the combined dataset:
    data1 <- merge (mesenteryF2, c(mesentery,bodywall,tentacle,pharynx,pharynx2,BW.ip))
    data1@assays$RNA@meta.features$NVE = genes$NVE
    rm (tentacle, mesentery, pharynx, bodywall, pharynx2,mesenteryF2,BW.ip)
    Tissues7.Raw=data1
    
    #save raw data for later:
    if (save.files){
      save (Tissues7.Raw,file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/Robj/Tissues.raw.Robj')
    }
  } 
## Seurat Analysis of full Tissue dataset ----
  if (Tissues7_analysis) 
  { 
    load (file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/Robj/Tissues.raw.Robj')
    
    Tissues7.Raw$orig.ident = as.factor(Tissues7.Raw$orig.ident)
    Tissues7 = Tissues7.Raw
    #  image the data in terms of gene content
     VlnPlot(object = Tissues7, features = c("nFeature_RNA", "nCount_RNA"), 
            group.by = 'orig.ident',ncol = 2,cols=LibCP)
    
    #  Filter dataset to include only the highest information containing cells:
    Tissues7 <- subset(x = Tissues7, subset =  nCount_RNA < 5000)#
    VlnPlot(object = Tissues7, features = c("nFeature_RNA", "nCount_RNA"), 
            ncol = 2,group.by = 'orig.ident',cols=LibCP)
    
    # Normalize the dataset and calculate variable features:  
    Tissues7 <- NormalizeData(object = Tissues7, scale.factor = 5000)
    
    # Select variable genes from each dataset separately then merge the list
    {
      list=  NULL
      vargenelist <- SplitObject(Tissues7, split.by = "orig.ident")
      for (i in 1:length(vargenelist)) {
        vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]], 
                                                selection.method = "vst",
                                                nfeatures = 500, verbose = FALSE)
      }
      
      for (i in 1:length(vargenelist)) {
        x <- vargenelist[[i]]@assays$RNA@var.features
        list=c(list,x)}
      list=unique(list)
      length(list)
      
    # set var.features and save copy of the list for recovery:
      Tissues7@assays$RNA@var.features = list
      Tissues7@misc$variable.genes.all = list
      
    } 
    
    # Scale the data
    Tissues7<- ScaleData(Tissues7,split.by = 'orig.ident')
    
    #  run different reduction algorithms:
    Tissues7 <- RunPCA(Tissues7, pcs.compute = 50, do.print = F) 
    #  choose number of PC dimensions from standard deviation:
    d=as.integer(which(Tissues7@reductions$pca@stdev>2))
    length(d)
    PCAPlot(Tissues7,group.by='orig.ident')
    
    #  UMAP  
    Tissues7 <- RunUMAP(Tissues7,  n.neighbors = 30L, 
                     min.dist = 0.3,
                     local.connectivity = 30,
                     spread = 0.6,reduction = 'pca', dims = d,
                     reduction.name = 'umap')
    
    # check output:
    
    DimPlot(Tissues7, reduction = 'umap', group.by = 'orig.ident',cols=LibCP)&
      FeaturePlot(Tissues7,c('NvSoxC','Nve-MyHC-st','NvGATA','NvFoxA'),
                  order = T,cols=gene.cp,label=T,reduction = 'umap')&
      NoAxes()&NoLegend() 
    
    #  cluster the cells    
    Tissues7 <- FindNeighbors(object = Tissues7,reduction ="pca",dims = d,
                              nn.method = 'annoy', annoy.metric = 'cosine',
                              k.param = 30)
    Tissues7 <- FindClusters(object = Tissues7,resolution = 0.1,random.seed = 0)#
    
    # re-order clusters
    Tissues7 <- BuildClusterTree(object = Tissues7, dims = c(d),reorder = TRUE, reorder.numeric = TRUE)
    
    #check output
    DimPlot(object = Tissues7,reduction = 'umap', pt.size = 0.5,
            cols =clust.cp.separate, label = T, label.size = 6)+NoAxes()
    
    # assign cluster ID
    clusterNames<- read_excel("D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/SupplementaryData1.xlsx", 
                              sheet = 'DataS8') 
    goi = clusterNames$Marker
    DotPlot(Tissues7,features = goi)+RotatedAxis()
    cl <-length(levels(Tissues7@active.ident))
    C.suffix <-seq(1:cl)
    
    # scale genes of interest:
    Tissues7<- ScaleData(Tissues7,split.by = 'orig.ident',model.use = 'linear', 
                         use.umi = F, features= goi)
    
    # make average matrix for clusters and genes
    
    g=length(goi)
    clName = vector()
    m=matrix(0L,g,cl)
    for (j in 1:cl)
    {
      for (i in 1:g)
        m[i,j]=mean(Tissues7@assays$RNA@scale.data[goi[i],WhichCells(Tissues7,idents = C.suffix[j])])
      clName[j]=as.integer(which.max(m[,j]))
    }
    sort(clName)
    clusterNames$ID[clName]
    
    DimPlot(Tissues7,cols=clust.cp.separate)&NoAxes()
    Tissues7@active.ident = factor(Tissues7@active.ident,
        levels(Tissues7@active.ident)[order(clName)])
    levels(Tissues7@active.ident) = clusterNames$ID[clName][order(clName)]
    
    #save the IDs in metadata:
    Tissues7@meta.data$IDs = Tissues7@active.ident
    
    DimPlot(Tissues7,label = T, cols = clust.cp.separate, repel = T,
            reduction='umap')+NoLegend()+NoAxes() 
    
    
    Tissues7 = Tissues7     #  save the dataset: 
    
    DimPlot(Tissues7,label = T,cols=clust.cp.separate)&NoAxes()
    
    goi=c('MYPH-like8','Nve-MyHC-st',
          'NvGATA','TCF24-like','Nve-MELC3','Nve-MELC5',
          'Nve-MRLC2','MTPN-like1','MP20-like3',
          'COA-like1','GELS1-like2','Nve-MELC4','NvVAX-EMX-like','NvFoxA'
    )
    Tissues7=Rmagic::magic(Tissues7,genes=goi) 
    
    if (save.files){
      save (Tissues7, file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/Robj/Tissues.Robj')
    }
  } 
  
## Seurat Analysis of "Muscle/Gastroderm" dataset ----
  if (Muscle7_analysis)
  {
    load (file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/Robj/Tissues.Robj')
    #select cells:
    muscle.cells = WhichCells(Tissues7,idents =  c('gastrodermis','retractor muscle'))
    #generate object:
    Muscle7 = CreateSeuratObject(Tissues7@assays$RNA@counts[,muscle.cells],)
    
    #keep previous data normalization and scaling
    Muscle7@assays$RNA@data=Tissues7@assays$RNA@data[,muscle.cells]
    Muscle7@assays$RNA@scale.data=Tissues7@assays$RNA@scale.data[,muscle.cells]
    Muscle7$orig.ident = Tissues7$orig.ident[muscle.cells]

    data1=Muscle7
    
    data1 <- FindVariableFeatures(data1)
    dim(data1)
    data1@misc$variable.genes.all=data1@assays$RNA@var.features
    
    
    # scale genes according to the full dataset:
    {
      #set the two dataset:
      data1.subset = data1
      data.full = Tissues7#AllData_Tissues#
      
      #Scale data in full from variable genes of subset
      coi = data1.subset@assays$RNA@counts@Dimnames[[2]] #pull out cells of interest
      t=ScaleData(data.full,split.by = 'orig.ident',
                  features = data1.subset@assays$RNA@var.features)
      # #import scaling to subset:
      data1.subset@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
      
      # re-set your working dataset
      data1=data1.subset
      
      #clean up the workspace
      rm(t, data1.subset,data.full)
    }
    
    data1 <- RunPCA(data1, pcs.compute = 50, do.print = F)
    ElbowPlot(object = data1, ndims = 50)#
    d = c(1:23)
    
    data1 <- RunUMAP(data1,  n.neighbors = 10L,spread =0.6,
                     reduction = 'pca',
                     dims = d,reduction.name ='umap',min.dist = 0.4,
                     local.connectivity = 10,return.model = T)# 
    DimPlot(data1, reduction = 'umap',cols=clust.cp.separate)+NoAxes()&
      
      DimPlot(data1, reduction = 'umap', group.by = 'orig.ident',cols=LibCP,
              order = (levels(data1@meta.data$orig.ident)))+NoAxes()+labs(title = 'tissue origin')
    
    data1 <- FindNeighbors(object = data1,reduction = "pca",k.param = 10,
   dims = as.integer(which(data1@reductions$pca@stdev > 2)),  
   nn.method = 'annoy',
   annoy.metric = 'cosine'
    )#   
    data1 <- FindClusters(object = data1,resolution = 0.3,random.seed = 0) 
    data1 <- BuildClusterTree(object = data1, reorder = TRUE, 
      reorder.numeric = TRUE, dims = d)
    
    DimPlot(data1, reduction = 'umap',cols=clust.cp.separate)+NoAxes()
    FeaturePlot(data1,c('MYPH-like8'),split.by = 'orig.ident',
                order = T,cols=gene.cp,label=T,reduction = 'umap'
    )&NoAxes()&NoLegend()
    
    # ID the clusters
    clusterNames<- read_excel("D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/SupplementaryData1.xlsx", 
      sheet = 'DataS9') 
    
    clusterNames=clusterNames
    goi = clusterNames$gene_short_name
    # assign cluster ID to the individual libraries
    data1<-ScaleData(data1,split.by = 'orig.ident',features=  goi)
    
    cl <-length(levels(data1@active.ident))
    C.suffix <-seq(1:cl)
    
    g=length(goi)
    clName = vector()
    m=matrix(0L,g,cl)
    for (j in 1:cl)
    {
      for (i in 1:g)
        m[i,j]=mean(data1@assays$RNA@scale.data[goi[i],WhichCells(data1,idents = C.suffix[j])])
      clName[j]=as.integer(which.max(m[,j]))
    }
    sort(clName)
    clusterNames$label[(clName)]
    clusterNames$label_merge[(clName)]
    DimPlot(data1, label = T, reduction = 'umap',cols=clust.cp.separate)+NoLegend()
    # arrange the clusters with the differentiated muscle populations first
    data1@active.ident = factor(data1@active.ident,
        levels(data1@active.ident)[order(clName)])
    #update names
    levels(data1@active.ident) = clusterNames$label[clName][order(clName)]
    data1$IDs.fine = data1@active.ident
    
    data1 <- SetIdent(data1, value = 'seurat_clusters')
    data1 <- BuildClusterTree(data1,reorder = T,reorder.numeric = T, dims= c(d))
    data1@active.ident = factor(data1@active.ident,
        levels(data1@active.ident)[order(clName)])
    
    levels(data1@active.ident) = clusterNames$label_merge[clName][order(clName)]
    DimPlot(object = data1,reduction = 'umap', pt.size = 1,
            cols =clust.cp.separate, label = T, label.size = 5)  +NoLegend() +
      DimPlot(object = data1,reduction = 'umap', pt.size = 1,group.by = 'IDs.fine',
              cols =clust.cp.separate, label = T, label.size = 5)  +NoLegend() 
    
    #save IDs in a metadata slot:
    data1$IDs = data1@active.ident
    # save the subset:
    Muscle7 = data1
    
    #' secondarily split BW muscle:
    coi=WhichCells(Muscle7,idents = 'endomes.BW.muscle')
    data1 = CreateSeuratObject(Muscle7@assays$RNA@counts[,coi],)
    
    data1 <- FindVariableFeatures(data1)
    dim(data1)
    data1@misc$variable.genes.all=data1@assays$RNA@var.features
    
    # scale genes according to the full dataset:
    {
      #set the two dataset:
      data1.subset = data1
      data.full = Tissues7#AllData_Tissues#
      
      #Scale data in full from variable genes of subset
      coi = data1.subset@assays$RNA@counts@Dimnames[[2]] #pull out cells of interest
      t=ScaleData(data.full,split.by = 'orig.ident',
                  features = data1.subset@assays$RNA@var.features)
      # #import scaling to subset:
      data1.subset@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
      
      # re-set your working dataset
      data1=data1.subset
      
      #clean up the workspace
      rm(t, data1.subset,data.full)
    }
    
    data1 <- RunPCA(data1, pcs.compute = 50, do.print = F)
    ElbowPlot(object = data1, ndims = 50)#
    d=as.integer(which(data1@reductions$pca@stdev>2))
    
    data1 <- RunUMAP(data1,  n.neighbors = 10L,spread =0.6,
                     reduction = 'pca',
                     dims = d,reduction.name ='umap',min.dist = 0.4,
                     local.connectivity = 1)# 
    DimPlot(data1, reduction = 'umap',cols=clust.cp.separate)+NoAxes()#&
    
    DimPlot(data1, reduction = 'umap', group.by = 'orig.ident',cols=LibCP,
            order = (levels(data1@meta.data$orig.ident)))+NoAxes()+labs(title = 'tissue origin')
    
    data1 <- FindNeighbors(object = data1,reduction = "pca",k.param = 10,
   dims = d,  
   nn.method = 'annoy',
   annoy.metric = 'cosine' )
    
    data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0) # 1 is necessary to split BW cM and PM
    data1 <- BuildClusterTree(object = data1, reorder = TRUE, 
      reorder.numeric = TRUE, dims = d)
    
    DimPlot(data1, reduction = 'umap',cols=clust.cp.separate)+NoAxes()
    
    clusterNames<- read_excel("D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/SupplementaryData1.xlsx", 
      sheet = 'DataS9') 
    clusterNames=clusterNames[1:17,] #remove the BW gene marker
    goi = clusterNames$gene_short_name
    # assign cluster ID to the individual libraries
    data1<-ScaleData(data1,split.by = 'orig.ident',features=  goi)
    
    cl <-length(levels(data1@active.ident))
    C.suffix <-seq(1:cl)
    
    g=length(goi)
    clName = vector()
    m=matrix(0L,g,cl)
    for (j in 1:cl)
    {
      for (i in 1:g)
        m[i,j]=mean(data1@assays$RNA@scale.data[goi[i],WhichCells(data1,idents = C.suffix[j])])
      clName[j]=as.integer(which.max(m[,j]))
    }
    sort(clName)
    clusterNames$label[(clName)]
    clusterNames$label_merge[(clName)]
    clName=c(4,5,5)
    # arrange the clusters with the differentiated muscle populations first
    data1@active.ident = factor(data1@active.ident,
        levels(data1@active.ident)[order(clName)])
    #update names
    levels(data1@active.ident) = clusterNames$label_merge[clName][order(clName)]
    
    DimPlot(object = data1,reduction = 'umap', pt.size = 1,
            cols =clust.cp.separate, label = T, label.size = 5)  +NoLegend()
    
    #add names to Muscle7: both fine and coarse
    Muscle7=SetIdent(Muscle7,value = 'IDs.fine')
    Muscle7=SetIdent(Muscle7, cells = coi, value = data1@active.ident[coi])
    levels(Muscle7@active.ident)
    Muscle7@active.ident = factor(Muscle7@active.ident,
          levels(Muscle7@active.ident)[c(3:5,1,2,6:15)])
    Muscle7$IDs.fine = Muscle7@active.ident
    
    Muscle7=SetIdent(Muscle7,value = 'IDs')
    Muscle7=SetIdent(Muscle7, cells = coi, value = data1@active.ident[coi])
    levels(Muscle7@active.ident)
    Muscle7@active.ident = factor(Muscle7@active.ident,
          levels(Muscle7@active.ident)[c(3,4,1,2,5:9)])
    Muscle7$IDs = Muscle7@active.ident
    
    
    
    DimPlot(Muscle7,group.by = 'IDs',cols=muscle.cp[c(1:7,8,13)],label=T)&NoLegend()&
      DimPlot(Muscle7,group.by = 'IDs.fine',cols=c('blue',muscle.cp),label=T)&NoLegend()
    
    goi=c('MYPH-like8','Nve-MyHC-st',
          'NvGATA','TCF24-like','Nve-MELC3','Nve-MELC5',
          'Nve-MRLC2','MTPN-like1','MP20-like3',
          'COA-like1','GELS1-like2','Nve-MELC4','NvVAX-EMX-like','NvFoxA'
    )
    Muscle7=Rmagic::magic(Muscle7,genes=goi) 
    if (save.files){
      save (Muscle7, file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/Robj/muscle.Robj') 
    }
  }

## generate nem64.dataset and Seurat Analysis ----
  if (Mutant)
  #generate Tissues object
  {  
    raw.nem64 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/NVE_mapping/Nem64_mutant/filtered_feature_bc_matrix")
    rownames(raw.nem64) <- genes$gene_short_name
    nem64  <- CreateSeuratObject(counts = raw.nem64, project = "nem64.mutant")
    nem64  <- RenameCells(nem64, add.cell.id = "nem64")
    
    nem64=nem64
    nem64@active.assay='RNA'
    nem64 <- NormalizeData(nem64,scale.factor = 5000)
    nem64 <- FindVariableFeatures(nem64,nfeatures = 2000, verbose = FALSE) 
    nem64 <- ScaleData(nem64)
    nem64 <- RunPCA(nem64)
    nem64 <- RunUMAP(nem64, n.neighbors = 10L, spread = 1, seed.use = 0,     
                     dims = 1:16, reduction.name = "umap", min.dist = 0.6, local.connectivity = 10)
    nem64=FindNeighbors(object = nem64, reduction = "pca", 
                        dims = c(1:12),  nn.method = "annoy", 
                        annoy.metric = "cosine", k.param = 40)
    nem64=FindClusters(nem64,resolution = 0.1)
    
    nem64=BuildClusterTree(nem64,reorder = T,reorder.numeric = T,dims = c(1:12))
    DimPlot(nem64, reduction = 'umap',cols=clust.cp.separate)+NoAxes()
    #use the same gene set as the Tissue analysis to annotate the clusters
    # assign cluster ID
    clusterNames<- read_excel("D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/SupplementaryData1.xlsx", 
                              sheet = 'DataS8') 
    goi = clusterNames$Marker
    DotPlot(nem64,features = goi)+RotatedAxis()
    cl <-length(levels(nem64@active.ident))
    C.suffix <-seq(1:cl)
    nem64<- ScaleData(nem64,split.by = 'orig.ident',model.use = 'linear', 
                      use.umi = F, features= goi)
    
    g=length(goi)
    clName = vector()
    m=matrix(0L,g,cl)
    for (j in 1:cl)
    {
      for (i in 1:g)
        m[i,j]=mean(nem64@assays$RNA@scale.data[goi[i],WhichCells(nem64,idents = C.suffix[j])])
      clName[j]=as.integer(which.max(m[,j]))
    }
    sort(clName)
    clusterNames$ID[clName]
    DimPlot(nem64,cols=clust.cp.separate)&NoAxes()
    nem64@active.ident = factor(nem64@active.ident,
                                levels(nem64@active.ident)[order(clName)])
    levels(nem64@active.ident) = clusterNames$ID[clName][order(clName)]
    
    #save the IDs in metadata:
    nem64@meta.data$IDs = nem64@active.ident
    
    DimPlot(nem64,label = T, cols = clust.cp.separate, repel = T,
            reduction='umap')+NoAxes() 
    
    if (save.files){
      save (nem64,file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/RObj/nem64.Robj')
    }  
  }
  ## generate nem64.endo subset ----
  {
    cells = WhichCells(nem64,idents =  c('gastrodermis'))# ,'intermuscular membrane'))
    nem64.endo = CreateSeuratObject(nem64@assays$RNA@counts[,cells],)
    nem64.endo$IDs <- 'nem64'
    #process this subset:
    nem64.endo=NormalizeData(nem64.endo,scale.factor = 5000)
    nem64.endo=FindVariableFeatures(nem64.endo)
  
  ### Transfer labels from wildtype dataset ----
    load (file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/Robj/muscle.Robj')
    Muscle7 = SetIdent(Muscle7,value='IDs')
    # re-set the variable features:
    Muscle7@assays$RNA@var.features=Muscle7@misc$variable.genes.all
    
    # transfer cell IDs and UMAP coordinates to the mutant:
    
    anchor.set=FindTransferAnchors(
      Muscle7,
      nem64.endo,
      normalization.method = "LogNormalize",
      #use mutant and WT variable features:
      features = intersect(nem64.endo@assays$RNA@var.features,Muscle7@assays$RNA@var.features), 
      reduction ='pcaproject',
      reference.reduction = 'pca',
      reference.assay	= 'RNA'
    )
    
    nem64.endo <- MapQuery(
      anchorset = anchor.set, #
      query = nem64.endo,
      reference = Muscle7,
      reduction.model = "umap",
      refdata='IDs'
    )
    
    nem64.endo<- SetIdent(nem64.endo,value = 'predicted.id')
    nem64.endo$IDs = nem64.endo@active.ident
    
    nem64.endo$IDs <- paste(nem64.endo@active.ident, "nem64m", sep=" | ")
    nem64.endo$IDs <- as.factor(nem64.endo$IDs)
    nem64.endo$IDs = factor(nem64.endo$IDs,
                              levels(nem64.endo$IDs)[c(7,8,1,6,4,3,5,2)])
    nem64.endo<- SetIdent(nem64.endo,value = 'IDs')
    
    DimPlot(nem64.endo,reduction = 'ref.umap',cols=muscle.cp.coarse[2:9])        
    
    if (save.files){
      save(nem64.endo,file='D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/RObj/nem64.gastrodermis.RObj')
    }
  }
  ### make merged Muscle dataset ---- 
  {
    Muscle7 = SetIdent(Muscle7,value = 'IDs')
    nem64.endo = SetIdent(nem64.endo,value = 'IDs')
    data1 <- merge(nem64.endo, Muscle7)
    data1@assays$RNA@var.features=Muscle7@misc$variable.genes.all
    lib.order=c('nem64',levels(as.factor(Muscle7$orig.ident)))
    lib.ind=match(lib.order,levels(as.factor(data1$orig.ident)))
    lib.ind
    data1$orig.ident=as.factor(data1$orig.ident)
    data1$orig.ident = factor(data1$orig.ident,
                              levels(data1$orig.ident)[lib.ind])

    orderM=c(levels(Muscle7),levels(nem64.endo))
    indM=match(orderM,levels(as.factor(data1$IDs)))
    indM

    data1@active.ident = factor(data1@active.ident,
                       levels(data1@active.ident)[indM])
    
    
    # data1<-SetIdent(data1,value = 'IDs')
    #Import UMAP from Seurat label transfer...    
    data1@reductions$refumap = merge(nem64.endo@reductions$ref.umap,Muscle7@reductions$umap)
    DimPlot(data1,reduction = 'refumap',cols=c(muscle.cp.coarse,muscle.cp.coarse[2:9]))
    
  
  nem64.Muscle7 = data1
  if (save.files){
    save(nem64.Muscle7,file = 'D:/Alison/manuscripts/MusclePaper/NatureCom.2023/accepted/RObj/nem64.muscle7.Robj')
  }
  }
  
  # clean up workspace:
  rm(data1,m,raw.nem64,Tissues7.Raw,vargenelist,anchor.set,clusterNames,lib.ind,
     lib.order,list,order,x,j,i,g,goi,ind,cl,cells,C.suffix,clName,coi,muscle.cells)
  