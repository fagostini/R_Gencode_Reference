## ----setupEnv, include=FALSE---------------------------------------------
working.folder = "/Users/agostif/Desktop/GencodeReference"

knitr::opts_knit$set(root.dir = working.folder)
knitr::opts_chunk$set(echo = FALSE)

setwd(working.folder)

chromosomes = paste0("chr", c(1:22, "X", "Y", "M"))

species = "Homo sapiens" # or 'Mus musculus' for mouse
genome = "hg38" # or mm10 for mouse
version = 27 # or M16 for mouse

## ----sessionInfo, echo=FALSE---------------------------------------------
cat(sessionInfo()$R.version$version.string, fill=TRUE)
cat(paste("Platform", sessionInfo()$platform), fill=TRUE)
cat(paste("Running under", sessionInfo()$running), fill=TRUE)
cat(paste("Last knitted on", date()), fill=TRUE)
cat(paste("Working directory set to", getwd()), fill=TRUE)
cat(paste("Species: ", species), fill=TRUE)
cat(paste("Genome assembly:", genome), fill=TRUE)
cat(paste("Gencode version:", version), fill=TRUE)

## ----loadPackages, message=FALSE, results="hide"-------------------------
require("AnnotationHub")
require("biomaRt")
require("GenomicFeatures")
require("rtracklayer")
require("data.table")
require("knitr")
require("ggplot2")

## ----downloadGencode, results="hide", message=FALSE, warning=FALSE-------
if( species %in% "Homo sapiens" ){
   if( !file.exists(paste0("gencode.v", version, ".annotation.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
                    paste0("gencode.v", version, ".annotation.gff3.gz"), quiet = TRUE)
   
   if( !file.exists(paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", version, "/gencode.v", version, ".long_noncoding_RNAs.gff3.gz"),
                    paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz"), quiet = TRUE)
   
   if( !file.exists(paste0("gencode.v", version, ".tRNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", version, "/gencode.v", version, ".tRNAs.gff3.gz"),
                    paste0("gencode.v", version, ".tRNAs.gff3.gz"), quiet = TRUE)
   
   if( !file.exists(paste0(genome, ".chrom.sizes")) )
      download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/bigZips/", genome, ".chrom.sizes"),
                    paste0(genome, ".chrom.sizes"), quiet = TRUE)
}else{
   if( !file.exists(paste0("gencode.v", version, ".annotation.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
                    paste0("gencode.v", version, ".annotation.gff3.gz"), quiet = TRUE)
   
   if( !file.exists(paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".long_noncoding_RNAs.gff3.gz"),
                    paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz"), quiet = TRUE)
   
   if( !file.exists(paste0("gencode.v", version, ".tRNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".tRNAs.gff3.gz"),
                    paste0("gencode.v", version, ".tRNAs.gff3.gz"), quiet = TRUE)
   
   if( !file.exists(paste0(genome, ".chrom.sizes")) )
      download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/bigZips/", genome, ".chrom.sizes"),
                    paste0(genome, ".chrom.sizes"), quiet = TRUE)
}

## ----readAnnotation, message=FALSE, warning=FALSE------------------------
chrominfo = sortSeqlevels(with(fread(paste0(genome, ".chrom.sizes"))[V1%in%chromosomes,],
                               Seqinfo(V1, seqlengths = V2, genome = genome)))

TxDb = makeTxDbFromGFF(file = paste0("gencode.v", version, ".annotation.gff3.gz"),
                       format = "gff3",
                       dataSource = paste("Gencode version", version),
                       organism = species,
                       chrominfo = chrominfo)

gff = import.gff3(paste0("gencode.v", version, ".annotation.gff3.gz"))
seqinfo(gff) = chrominfo

# Remove the duplicated genes on the chromosome Y
x.genes = unique(gff[seqnames(gff)%in%"chrX"]$gene_id)
y.genes = unique(gff[seqnames(gff)%in%"chrY"]$gene_id)
y.genes = y.genes[y.genes%in%x.genes]
gff = gff[!(seqnames(gff)%in%"chrY" & gff$gene_id%in%y.genes)]

# # Remove the ribosomal RNA genes
# ribo.gff = gff[gff$gene_type%in%"rRNA",]
# gff = gff[!gff$gene_type%in%"rRNA",]

lnc.gff = import.gff3(paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz"))
seqinfo(lnc.gff) = chrominfo

# Remove the duplicated genes on the chromosome Y
x.genes = unique(lnc.gff[seqnames(lnc.gff)%in%"chrX"]$gene_id)
y.genes = unique(lnc.gff[seqnames(lnc.gff)%in%"chrY"]$gene_id)
y.genes = y.genes[y.genes%in%x.genes]
lnc.gff = lnc.gff[!(seqnames(lnc.gff)%in%"chrY" & lnc.gff$gene_id%in%y.genes)]

# if(interactive()) {
#    saveDb(TxDb, file=paste0(genome, "_Gencode", version, ".sqlite"))
# }

## ------------------------------------------------------------------------
tab = data.table(as.data.frame(gff))[type%in%"gene", .N, by=c("gene_type", "level")]

tab = dcast.data.table(tab, gene_type~level, value.var="N", fill=0)[order(-(`1`+`2`+`3`))]

setnames(tab, "gene_type", "Biotype")

kable(head(data.frame(tab), 20), col.names=colnames(tab), format='markdown', digits=1,
      caption="Table 1: Gene biotypes (top 20) per annotation level.")

## ----extractMetadata-----------------------------------------------------
# Genes, transcripts and exons
genes = genes(TxDb)
genes = genes[!grepl("PAR", names(genes)),]
txs = transcriptsBy(TxDb, by="gene")
txs = txs[!grepl("PAR", names(txs)),]
exons = exonsBy(TxDb, by="tx", use.names=TRUE)
exons = exons[!grepl("PAR", names(exons)),]

# Gene biotypes
bio.names = names(sort(table(gff[gff$type%in%"gene"]$gene_type), decreasing=TRUE))
geneBT = lapply(bio.names,
                function(x)
                   unique(gff$gene_id[gff$gene_type%in%x]) )
names(geneBT) = bio.names

gene.metadata = data.table(as.data.frame(gff))[type%in%"gene", 
                                               list(gene_id, gene_type, gene_name, level, seqnames, start, end, width, strand)]

# Transcripts biotypes
bio.names = names(sort(table(gff[gff$type%in%"transcript"]$transcript_type), decreasing=TRUE))
txsBT =  lapply(bio.names,
                function(x) unique(gff$transcript_id[gff$transcript_type%in%x]) )
names(txsBT) = bio.names

txs.metadata = data.table(as.data.frame(gff))[type%in%"transcript",
                                              list(transcript_id, transcript_type, gene_id, gene_type, gene_name, level, seqnames, start, end, width, strand)]

## ----longestGeneTxs------------------------------------------------------
# GRanges of protein-coding and ncRNA genes
protein.coding = names(geneBT)%in%c("protein_coding",
                                    paste(rep(c("IG","TR"), each=4), c("C","D","J","LV","V"), "gene", sep="_"))
protein.coding = as.vector(unlist(geneBT[protein.coding]))
genes.pc.granges = genes[names(genes) %in% protein.coding,]

ncRNA = unique(c(grep("RNA", names(geneBT), value=TRUE),
                 unique(lnc.gff[lnc.gff$type%in%"gene"]$gene_type)))
ncRNA = as.vector(unlist(geneBT[ncRNA]))
genes.nc.granges = genes[names(genes) %in% ncRNA,]

# Get the transcript lengths
txs.length = transcriptLengths(TxDb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
protein.coding = names(txsBT) %in% c("protein_coding",
                                           paste(rep(c("IG","TR"), each=4), c("C","D","J","V"), "gene", sep="_"))
protein.coding = as.vector(unlist(txsBT[protein.coding]))
txs.length.pc = txs.length[txs.length$tx_name %in% protein.coding,]

ncRNA = unique(c(grep("RNA", names(txsBT), value=TRUE),
                 unique(lnc.gff[lnc.gff$type%in%"gene"]$gene_type)))
ncRNA = as.vector(unlist(txsBT[ncRNA]))
txs.length.nc = txs.length[txs.length$tx_name %in% ncRNA,]

# Order by size
txs.length.pc = txs.length.pc[order(txs.length.pc$tx_len, decreasing=TRUE),]
txs.length.nc = txs.length.nc[order(txs.length.nc$tx_len, decreasing=TRUE),]

# Keep only the longest transcript per gene
txs.longest.pc = txs.length.pc[!duplicated(txs.length.pc$gene_id),]
txs.longest.nc = txs.length.nc[!duplicated(txs.length.nc$gene_id),]

# GRanges of protein-coding and ncRNA transcripts
txs.longest.pc.pre.granges = unlist(txs)
txs.longest.pc.pre.granges = txs.longest.pc.pre.granges[txs.longest.pc.pre.granges$tx_name %in% txs.longest.pc$tx_name,]
mcols(txs.longest.pc.pre.granges) = txs.longest.pc[match(txs.longest.pc.pre.granges$tx_name, txs.longest.pc$tx_name),]
names(txs.longest.pc.pre.granges) = NULL
txs.longest.pc.pre.granges = split(txs.longest.pc.pre.granges, txs.longest.pc.pre.granges$gene_id)

txs.longest.nc.pre.granges = unlist(txs)
txs.longest.nc.pre.granges = txs.longest.nc.pre.granges[txs.longest.nc.pre.granges$tx_name %in% txs.longest.nc$tx_name,]
mcols(txs.longest.nc.pre.granges) = txs.longest.nc[match(txs.longest.nc.pre.granges$tx_name, txs.longest.nc$tx_name),]
names(txs.longest.nc.pre.granges) = NULL
txs.longest.nc.pre.granges = split(txs.longest.nc.pre.granges, txs.longest.nc.pre.granges$gene_id)

# GRanges of protein-coding and ncRNA transcripts
# Including the intronic regions (pre-RNA)
txs.longest.pc.rna.granges = unlist(exons)
txs.longest.pc.rna.granges = txs.longest.pc.rna.granges[names(txs.longest.pc.rna.granges) %in% txs.longest.pc$tx_name,]
mcols(txs.longest.pc.rna.granges) = cbind(mcols(txs.longest.pc.rna.granges), txs.longest.pc[match(names(txs.longest.pc.rna.granges), txs.longest.pc$tx_name),])
names(txs.longest.pc.rna.granges) = NULL
txs.longest.pc.rna.granges = split(txs.longest.pc.rna.granges, txs.longest.pc.rna.granges$gene_id)

# Excluding the intronic regions (RNA)
txs.longest.nc.rna.granges = unlist(exons)
txs.longest.nc.rna.granges = txs.longest.nc.rna.granges[names(txs.longest.nc.rna.granges) %in% txs.longest.nc$tx_name,]
mcols(txs.longest.nc.rna.granges) = cbind(mcols(txs.longest.nc.rna.granges), txs.longest.nc[match(names(txs.longest.nc.rna.granges), txs.longest.nc$tx_name),])
names(txs.longest.nc.rna.granges) = NULL
txs.longest.nc.rna.granges = split(txs.longest.nc.rna.granges, txs.longest.nc.rna.granges$gene_id)

exonsByGene = exonsBy(TxDb, by="gene")
# Remove duplicated exon entries within individual genes
exonsByGene = exonsByGene[!duplicated(exonsByGene)]
exonsByTxs = exonsBy(TxDb, by="tx", use.names=TRUE)

# Merge the coding exons for each protein-coding gene
temp = unlist(exonsByTxs[names(exonsByTxs) %in% protein.coding])
temp$gene_id = gff[match(names(temp), gff$transcript_id),]$gene_id
names(temp) = NULL
exons.pc.granges = reduce(split(temp, temp$gene_id))

# Merge the non-coding exons for each non-coding gene
temp = unlist(exonsByTxs[names(exonsByTxs) %in% ncRNA])
temp$gene_id = gff[match(names(temp), gff$transcript_id),]$gene_id
names(temp) = NULL
exons.nc.granges = reduce(split(temp, temp$gene_id))

save(genes, txs,
     gene.metadata, txs.metadata,
     exonsByGene, exonsByTxs,
     genes.pc.granges, genes.nc.granges, 
     exons.pc.granges, exons.nc.granges,
     txs.longest.pc.pre.granges, txs.longest.nc.pre.granges, 
     txs.longest.pc.rna.granges, txs.longest.nc.rna.granges, 
     file=paste0(genome, "_Gencode", version, "_annotations.RData"))

## ------------------------------------------------------------------------
tab = gene.metadata

bio_order = tab[, list(median=median(width)), by="gene_type"][order(-median), gene_type]

tab[, gene_type := factor(gene_type, bio_order)]

gg = ggplot(tab, aes(x=gene_type, y=width)) +
    ggtitle("Gene length distribution") +
    geom_boxplot() +
    scale_y_log10("Gene width") +
    scale_x_discrete("") +
    theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave(gg, filename="img/gene_length.png", width=12, height=6)

tab = txs.metadata

bio_order = tab[, list(median=median(width)), by="gene_type"][order(-median), gene_type]

tab[, gene_type := factor(gene_type, bio_order)]

gg = ggplot(tab, aes(x=gene_type, y=width)) +
    ggtitle("Trascript length distribution (Longest transcript per gene)") +
    geom_boxplot() +
    scale_y_log10("Transcript width") +
    scale_x_discrete("") +
    theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave(gg, filename="img/txs_length.png", width=12, height=6)

## ----codingGeneRegions---------------------------------------------------
cds = cdsBy(TxDb, by=c("tx"), use.names=TRUE)
cds = cds[!grepl("PAR", names(cds)),]
utr5 = fiveUTRsByTranscript(TxDb, use.names=TRUE)
utr5 = utr5[!grepl("PAR", names(utr5)),]
utr3 = threeUTRsByTranscript(TxDb, use.names=TRUE)
utr3 = utr3[!grepl("PAR", names(utr3)),]

txs.sel = txs.longest.pc$tx_name

# function to subset and unlist
convert <- function(gr, region){
    # select transcripts
    gr = gr[names(gr)%in%txs.sel]
    # add tx.id
    gr = unlist(gr)
    mcols(gr) = NULL
    mcols(gr)$tx.id = names(gr)
    names(gr) = NULL
    # add region type
    mcols(gr)$region = region
    return(gr)
}

txs.longest.pc.regions = c(convert(cds, "CDS"), convert(utr5, "UTR5"), convert(utr3, "UTR3"))

txs.longest.pc.regions = sort.GenomicRanges(txs.longest.pc.regions, ignore.strand=TRUE)

txs.longest.pc.regions = split(txs.longest.pc.regions, txs.longest.pc.regions$tx.id)

save(txs.longest.pc.regions, file=paste0(genome, "_Gencode", version, "_annotations.pc.transcript.regions.RData"))

## ------------------------------------------------------------------------
tab = data.table(as.data.frame(unlist(txs.longest.pc.regions), row.names=NULL))

setnames(tab, "region", "Region")

tab = tab[, sum(width), by=c("tx.id", "Region")][, c(as.list(summary(V1)), Coverage=sum(V1)), by="Region"]

kable(data.frame(tab), col.names=colnames(tab), format='markdown', digits=1,
      caption="Table 2: Summary of the non-overlapping protein-coding regions.")

## ----genomicRegions------------------------------------------------------
# get transcript regions (all genes)
# hierarchy: ncRNA > cds > utr3 > utr5 > intron > other > intergenic

# identify protein-coding genes, ncRNAs and other
protein.coding = which(names(txsBT) %in% c("protein_coding",
                                           paste(rep(c("IG","TR"),each=4), c("C","D","J","V"), "gene", sep="_")))
ncRNA =  unique(c(grep("RNA", names(geneBT), value=TRUE),
                  unique(lnc.gff[lnc.gff$type%in%"gene"]$gene_type)))
other = setdiff( 1:length(txsBT), c(ncRNA, protein.coding) )
protein.coding = as.vector(unlist(txsBT[protein.coding]))
ncRNA = as.vector(unlist(txsBT[ncRNA]))
other = as.vector(unlist(txsBT[other]))

# ncRNAs: get transcripts
ncRNA.tx = unlist(exonsBy(TxDb, by="tx", use.names=TRUE))
ncRNA.tx = subset(ncRNA.tx, names(ncRNA.tx) %in% ncRNA)
ncRNA.tx = reduce(ncRNA.tx)
mcols(ncRNA.tx)$region = "ncRNA"

# # Protein_coding genes: get UTR, CDS and intron coordinates
# tx.pc = ensemblTable$ensembl_transcript_id[ensemblTable$ensembl_gene_id %in% protein.coding]

# CDS
cds = unlist(cdsBy(TxDb, by="tx", use.names=TRUE))
cds = subset(cds, names(cds) %in% protein.coding)
cds = reduce(cds)
cds = setdiff(cds, ncRNA.tx)
mcols(cds)$region = "CDS"

# 3'UTR
utr3 = threeUTRsByTranscript(TxDb, use.names=TRUE)
utr3 = subset(utr3, names(utr3) %in% protein.coding)
utr3 = reduce( unlist(utr3))
utr3 = setdiff(utr3, union(cds, ncRNA.tx))
mcols(utr3)$region = "UTR3"

# 5'UTR
utr5 = fiveUTRsByTranscript(TxDb, use.names=TRUE)
utr5 = subset(utr5, names(utr5) %in% protein.coding)
utr5 = reduce( unlist(utr5))
utr5 = setdiff(utr5, union(cds, ncRNA.tx))
utr5 = setdiff(utr5, utr3)
mcols(utr5)$region = "UTR5"

# introns
introns = intronsByTranscript(TxDb, use.names=TRUE)
introns = subset(introns, names(introns) %in% protein.coding)
introns = reduce(unlist(introns))
introns = setdiff(introns, union(cds, ncRNA.tx))
introns = setdiff(introns, union(utr3, utr5))
mcols(introns)$region = "intron"

# other
other.tx = unlist(transcriptsBy(TxDb, by="gene"))
other.tx = subset(other.tx, mcols(other.tx)$tx_name %in% other)
other.tx = reduce(unlist(other.tx))
other.tx = setdiff(other.tx, union(cds, ncRNA.tx))
other.tx = setdiff(other.tx, union(utr3, utr5))
other.tx = setdiff(other.tx, introns)
mcols(other.tx)$region = "other"

# combine all regions (in order of hierarchy)
tx.regions = c(ncRNA.tx, cds, utr3, utr5, introns, other.tx)
mcols(tx.regions)$region = factor(mcols(tx.regions)$region, levels=c("ncRNA", "CDS", "UTR3", "UTR5", "intron", "other"))

# # sanity check
# length(tx.regions) == length(reduce(tx.regions, min.gapwidth = 0L))

save(tx.regions, file=paste0(genome, "_Gencode", version, "_annotations.all.genes.transcript.regions.RData"))

## ------------------------------------------------------------------------
tab = data.table(as.data.frame(unlist(tx.regions), row.names=NULL))

setnames(tab, "region", "Region")

tab = tab[, c(as.list(summary(width)), Coverage=sum(width)), by="Region"]

kable(data.frame(tab), col.names=colnames(tab), format='markdown', digits=1,
      caption="Table 3: Summary of all non-overlapping genomic regions.")

