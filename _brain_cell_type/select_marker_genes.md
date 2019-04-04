# Selecting top 100 genes comment

If we selected marker genes per cell type on one variable or criteria, we can simply order all genes by the variable and select. However, we'd like to select 100 marker genes on two criteria, specifically 

```{r}
logFC > logFC_thres and FDR < FDR_thres
```

As a remark, the current list of genes considered already satisfy pre-defined logFC thresholds and FDR thresholds. From among this list, we'd like to further select optimal genes on both criteria, we can do a grid search with one quantile threshold to select the intersection of top genes. In R code, we'd like to obtain a quantile threshold "qq" such that

```{r}
sum( logFC > quantile(logFC,qq) & FDR < quantile(FDR,1-qq) )
```

is about 100 per cell type.
