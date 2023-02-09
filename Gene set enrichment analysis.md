In a typical Rna-seq experiment thousands of genes are found to be differentially expressed (DE). Each gene may be involved in one or more biological processes, may be part of one or more pathways, or may share and/or part of more than one relevent biological phenomena. If we just go about documenting the effect caused by each gene individually to make the whole story, it will certainly be a cumborsome task if not impossible, and that too with a risk of concluding a contrasting/confusing results. To save us of this dilemma some smart people from [GO consortium](http://geneontology.org/)  came with a brilliant idea. The idea was to create a set of genes (a collection of genes) with a shared biological property, a comman pathways, and unsurprisingly a comman GO term etc. So instead of looking at the Diffrential Expression of individual genes between two groups (e.g. case and control) of samples, we look for the what GO Terms are enriched in our experiments. This gives us a broader view of the results at hand. Looking at the enriched GO terms we can make some conclusion about what pathway/process/components have been affected and draw a mechanistic perpective of the experiment. So GSEA essentially is a kind of statitical test (more aptly known as enrichment analyses) of functional terms that appear associated to the given set of differentially expressed genes more often than expected by chance. 

So lets begin by selecting genes that are significantly DE between CTRL and Case samples. Let’s extract genes that have an adjusted p-value below 0.1 and that show a 2-fold change (either negative or positive) in the case compared to CTRL.

```r
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
# eliminating all the genes for which p.values is NA
DE <- DEresults[!is.na(DEresults$padj),]

# Extracting all the genes having p.values < 0.1
DE <- DE[DE$padj < 0.1,]

# Extracting all the genes having a log fold change > 2
DE <- DE[abs(DE$log2FoldChange) > 1,]


```


