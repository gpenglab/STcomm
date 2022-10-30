
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STcomm

<!-- badges: start -->
<!-- badges: end -->

Welcome to STcomm, an R package to illustrate the spatially resolved
cell interactions by combined the spatial cellular colocalization with
their enriched ligand-receptor co-expression patterns inferred from both
spatial and single-cell transcriptomic data.

## Installation

STcomm R package can be easily installed from Github using devtools:

``` r
devtools::install_github("Vanjia-lee/STcomm")
```

## Quick Guide to Getting Started with stComm

Firstly, `rctd4weights` function will help you tidy the confident
weights for every pixel after RCTD analysis.

``` r
library(STcomm)
rctd.multi <- rctd4weight(rctd_obj, rctd_mode = 'multi', conf = TRUE)
```

Secondly, `cellColocation` function will help you quantify the colocalization 
of cell type pairs within spots by calculating Pearson correlation coefficient 
or Jaccard similarity coefficient based on cell type composition predicted by 
RCTD.

```r
# By calculating Pearson correlation coefficient based on the RCTD object  
colocal_ctps1 <- cellColocation(rctd_obj, rctd_mode = 'multi', 
                                method = 'pcc', pcc = 0.06, 
                                pval = 0.05, padj = 0.05)

# By calculating Pearson correlation coefficient based on the data.frame of RCTD 
confident weights (with conf = TRUE when RCTD) for each pixel
colocal_ctps2 <- cellColocation(rctd_df, method = 'pcc', 
                               pcc = 0.06, pval = 0.05, padj = 0.05)

# By calculating Jaccard similarity coefficient based on the RCTD object
colocal_ctps3 <- cellColocation(rctd_obj, rctd_mode = 'multi', 
                                method = 'jac', jac = 0.05)

# By calculating Jaccard similarity coefficient based on the data.frame of RCTD 
confident weights (with conf = TRUE when RCTD) for each pixel
colocal_ctps4 <- cellColocation(rctd_df, method = 'jac', jac = 0.05)
```

Then, you can identify significant co-occurrence cell type groups belonging to 
the same spot from the cell type colocalization network. 

Thirdly, based on the spatial data, you can obtain significantly co-expressed 
Ligand-Receptor (LR) pairs for spatially co-localized cell types by performing 
Fisher's exact test on binarized co-localized cell type pairs and co-expressed 
LR pairs at spot level. Next, you calculate significant communication between 
LR pairs in co-localized cell type pairs based on the refrence single cell 
transcriptomic data. Finally, to get high confidence and spatial aware cell-cell 
communication, you can keep only spatially relevant communication information 
based on the above Fisher exact significancy. `st_comm` function can help you 
characterize confident spatially resolved cell-cell interaction with the tissue 
organization.

```r
# load the CellChat object for the refrence single cell transcriptomic data
cellchat <- readRDS(cellchat)
st_net <- st_comm(st_obj, weights.df = rctd_multi, ctpairs = colocal_ctps1, cellchat = cellchat)

# or you would like to prepare a data.frame tidyed frome the CellChat object
net.df <- subsetCommunication(cellchat)
net.df$ct_pairs <- paste0(net.df$source, "_", net.df$target)
st_net <- st_comm(st_obj, weights.df = rctd_multi, ctpairs = colocal_ctps1, cellchat = net.df)
```
