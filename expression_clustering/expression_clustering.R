source("https://bioconductor.org/biocLite.R")
# biocLite("sva")
# biocLite("ballgown")
library(ballgown)
library(sva)
library(pheatmap)
library(mixtools)

# Extract expression estimates
data_directory = Sys.getenv('DATA_DIR')
setwd(data_directory)
bg = ballgown(dataDir=data_directory, samplePattern='^SCLC', meas='all')
# save(bg, file='bg.rda')
# bg = load(file='bg.rda')
gene_fpkm = gexpr(bg)
genes = rownames(gene_fpkm)
samples = colnames(gene_fpkm)

# Annotate by library prep
relapse = samples[!grepl('tumor', samples)]
relapse.numbers = sub('FPKM.SCLC([[:digit:]]+)_.*','\\1', relapse, perl=TRUE)
relapse.first = relapse[as.integer(relapse.numbers) < 19]
relapse.second = relapse[as.integer(relapse.numbers) > 18]

# Filter by expressed status
relapse.first.expressed = apply(gene_fpkm[,relapse.first] > 1, 1, (function (x) table(x)['TRUE']))
relapse.first.expressed[is.na(relapse.first.expressed)] = 0
relapse.second.expressed = apply(gene_fpkm[,relapse.second] > 1, 1, (function (x) table(x)['TRUE']))
relapse.second.expressed[is.na(relapse.second.expressed)] = 0
expressed_genes = relapse.first.expressed >= 1 & relapse.second.expressed >= 1
expressed_gene_fpkm = log2(gene_fpkm[expressed_genes,] + 1)

# Normalize across batches
m = matrix(nrow=18, ncol=1)
batch = data.frame(m, row.names = relapse)
batch[relapse.first, 1] = 1
batch[relapse.second, 1] = 2
batch.colnames = c('batch')
modcombat = model.matrix(~1, data=batch)
combat_edata = ComBat(dat=expressed_gene_fpkm, batch=batch$m, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
save(combat_edata, file = 'corrected_expr.rda')

# ENSG IDs for DOI: 10.1016/j.ccell.2016.12.005
neuro.marker.genes = data.frame(
  ensembl = c(
    'ENSG00000139352',
    'ENSG00000162992',
    'ENSG00000102003',
    'ENSG00000173404',
    'ENSG00000100604',
    'ENSG00000089199',
    'ENSG00000171951',
    'ENSG00000134443',
    'ENSG00000149294',
    'ENSG00000154277',
    'ENSG00000110680',
    'ENSG00000175868'
  ),
  HGNC = c(
    'ASCL1',
    'NEUROD1',
    'SYP',
    'INSM1',
    'CHGA',
    'CHGB',
    'SCG2',
    'GRP',
    'NCAM1',
    'UCHL1',
    'CALCA',
    'CALCB'
  )
)

neuro.marker.genes.filtered = neuro.marker.genes[neuro.marker.genes$ensembl %in% rownames(combat_edata),]
neuro.marker.fpkm = combat_edata[as.character(neuro.marker.genes.filtered$ensembl),]

colnames(neuro.marker.fpkm) = gsub('^FPKM.', '', colnames(neuro.marker.fpkm))
rownames(neuro.marker.fpkm) = neuro.marker.genes.filtered$HGNC

pheatmap(neuro.marker.fpkm, cluster_rows = FALSE, cutree_cols = 3)

### The below code may be used to make annotation tracks that annotate genes with high/low expression status
# intersect <- function(m1, s1, m2, s2, prop1, prop2){    
#   B <- (m1/s1^2 - m2/s2^2)
#   A <- 0.5*(1/s2^2 - 1/s1^2)
#   C <- 0.5*(m2^2/s2^2 - m1^2/s1^2) - log((s1/s2)*(prop2/prop1))
#   (-B + c(1,-1)*sqrt(B^2 - 4*A*C))/(2*A)
# }
# 
# 
# high_low <- function(edata){
#   m = normalmixEM(edata, k=2)
#   i = intersect(m$mu[1], m$sigma[1], m$mu[2], m$sigma[2], m$lambda[1], m$lambda[2])
#   low = min(m$mu[1], m$mu[2])
#   high = max(m$mu[1], m$mu[2])
#   if (i[1] > low && i[1] < high){
#     b = i[1]
#   }else if(i[2] > low && i[2] < high){
#     b = i[2]
#   }else{
#     stop('Bad intersects.')
#   }
#   return(ifelse(m$x > b, 'High', 'Low'))
# }

# x.apc = high_low(combat_edata['ENSG00000134982', ])
# x.myc = high_low(combat_edata['ENSG00000136997', ])
# x.ctnnb1 = high_low(combat_edata['ENSG00000168036', ])
# x.axin2 = high_low(combat_edata['ENSG00000168646',])
# x.mycl = high_low(combat_edata['ENSG00000116990', ])
# x.mycn = high_low(combat_edata['ENSG00000134323', ])
# x.dll3 = high_low(combat_edata['ENSG00000090932', ])
# x.slfn11 = high_low(combat_edata['ENSG00000172716', ])
# x.ezh2 = high_low(combat_edata['ENSG00000106462', ])
# x.cd44 = high_low(combat_edata['ENSG00000026508', ])
# x.hes1 = high_low(combat_edata['ENSG00000114315', ])
# x.rest = high_low(combat_edata['ENSG00000084093', ])
# x.notch1 = high_low(combat_edata['ENSG00000148400', ])
# x.notch2 = high_low(combat_edata['ENSG00000134250', ])
# x.notch3 = high_low(combat_edata['ENSG00000074181', ])

# annotation_col = data.frame(
#   CD44 = x.cd44,
#   HES1 = x.hes1,
#   REST = x.rest,
#   NOTCH1 = x.notch1,
#   NOTCH2 = x.notch2,
#   NOTCH3 = x.notch3
#   Myc = x.myc,
#   Apc = x.apc,
#   Ctnnb1 = x.ctnnb1,
#   Axin2 = x.axin2,
#   Mycl = x.mycl,
#   Mycn = x.mycn,
#   Dll3 = x.dll3
# )

# rownames(annotation_col) = colnames(neuro.marker.fpkm)
# 
# ann_colors = list(
#   CD44 = c(High = 'red', Low = 'blue'),
#   HES1 = c(High = 'red', Low = 'blue'),
#   REST = c(High = 'red', Low = 'blue'),
#   NOTCH1 = c(High = 'red', Low = 'blue'),
#   NOTCH2 = c(High = 'red', Low = 'blue'),
#   NOTCH3 = c(High = 'red', Low = 'blue')
#   Myc = c(High = 'red', Low = 'blue'),
#   Apc = c(High = 'red', Low = 'blue'),
#   Ctnnb1 = c(High = 'red', Low = 'blue'),
#   Axin2 = c(High = 'red', Low = 'blue'),
#   Mycl = c(High = 'red', Low = 'blue'),
#   Mycn = c(High = 'red', Low = 'blue'),
#   Dll3 = c(High = 'red', Low = 'blue')
# )
