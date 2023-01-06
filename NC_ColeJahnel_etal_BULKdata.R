
#' Calculate Transcripts Per X over a matrix
#'
#' General verison of the standard transcripts per million (TPM) values for each gene over a matrix. The columns
#' always add up to 1e+06.
#'
#'
#' @param countsMat Matrix of counts
#' @param lengths vector of gene lengths corresponding to the rows in R
#' @param x normalization factor (e.g. 1e6 for TPM)
#' @return transformed matrix with the TPX values
#' @export
calculateTPX <- function(countsMat, lengths, x) {
  kblengths <- lengths/1e3
  R <- countsMat/kblengths
  Z <- colSums(R)/x
  t(t(R)/Z)
}

# Here is the lower-order TPM function for simplicity. My guess is that we will have to adjust the X in TPX to fit both distributions.
calculateTPM <- function(countsMat, lengths) calculateTPX(countsMat, lengths, 1e6)

# First let's derive speaking names so we can use them if we want to point them out in the table
gene.info <- read.table('D:/Dropbox/Nematostella/analysis/R_AGC/nve_master_ortholog_table.20180724.txt',
                        header=T, quote = '', sep='\t', comment.char = '', stringsAsFactors = F, strip.white=T
                        )[,c('GeneId','deM_TF','Curated.ID','Trinotate.Name','Emapper.Name')]

gene.info$SpeakingName <- with(
  gene.info,
  ifelse(nchar(Curated.ID), Curated.ID,
         ifelse(nchar(Emapper.Name), Emapper.Name,
                ifelse(nchar(Trinotate.Name), Trinotate.Name, GeneId))))

# Next let's get the lengths of the genes so we can use them for the TPX calculations
nves <- read.table('D:/Dropbox/Nematostella/analysis/R_AGC/nveGenes.good.130208.longCDS.bed', quote = '', sep='\t',
                   col.names=c('chrom','start','end','name','score','strand',
                              'cdsStart','cdsEnd','itemRGB',
                              'blockCount','blockStarts','blockSizes'))[,c(2,3,4)]
nves$length <- nves$end-nves$start

# Now merge both the length and gene speaking names into a single table
gene.info.length <- merge(gene.info[,c('GeneId','deM_TF','SpeakingName')], nves[,c('name','length')],
                   by.x='GeneId',by.y='name')

# Add the DE info (see reproduce_stefan.R for how this table was generated)
gene.info.length.de <- merge(gene.info.length, read.table('D:/Dropbox/Nematostella/analysis/R_AGC/unfiltered_results.txt'), by.x='GeneId', all.x=T, by.y=0)


# Here is where we begin handling the counts themselves
counts <- read.table('D:/Dropbox/Nematostella/analysis/R_AGC/htseq_combined_mod.counts.matrix', skip=1, row.names=1, com='',
                     col.names=c('','muscle1','non.muscle1','muscle2','non.muscle2'))[,c(1,3,2,4)]
ann.counts <- merge(gene.info.length.de, counts, by.x='GeneId', by.y=0)

# calculate a matrix of tpm values
tpm <- calculateTPX(ann.counts[,grep('muscle',colnames(ann.counts))],
                    ann.counts$length, 1e6)

# Get the mean TPM values for each of the cell types
mean.tpm <- sapply(c('muscle','non.muscle'),
                   function(n) rowMeans(tpm[,grep(paste0('^',n), colnames(tpm))]))
# Take the logarithm
lmean.tpm <- log2(1+mean.tpm)

# Combine the gene info and the final gene info (including the most relevant columns)
ann.mean.tpm <- as.data.frame(cbind(gene.info.length.de[,c('GeneId','SpeakingName','deM_TF','logFC','FDR')], lmean.tpm))


# here is an example of how I might plot it
library(ggplot2)
library(reshape2)

# before melting, we would subset the table here -- this is an example of just plotting transcription factors
genes.to.plot <- which(nchar(ann.mean.tpm$deM_TF) > 0 & ann.mean.tpm$logFC < -2)
ann.mean.tpm.m <- melt(ann.mean.tpm[genes.to.plot,],
                       measure.vars=c('muscle','non.muscle'), variable.name='CellType')

# this is arnau's color scheme and I like it
cols <- c("white","orange","red","purple","black")

ggplot(ann.mean.tpm.m,aes(x=CellType, y=SpeakingName)) +
  geom_raster(aes(fill=value)) +
  theme_bw() +
  scale_x_discrete(expand=c(0,0)) +
  theme(axis.title = element_blank(),# axis.text.y=element_blank(),
        panel.border=element_rect(color="black",fill=NA,size=1),
        legend.position = 'top', legend.margin=margin(0,0,0,0),
        legend.justification='center')+
  scale_fill_gradientn(colours=cols,
                       values=c(0,seq(.Machine$double.eps,1,length.out=length(cols)-1))) +
  guides(fill=guide_colorbar(title='Normalized log expression',
                             title.position='top',
                             title.theme = element_text(size=10, hjust=0.5),
                             ticks=F, barheight = 0.5, barwidth=8,
                             frame.colour='black',frame.linewidth = 1,
                             label.theme = element_text(size=9)))

