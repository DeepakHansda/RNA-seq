In a typical Rna-seq experiment thousands of genes are found to be differentially expressed (DE). Each gene may be involved in one or more biological processes, may be part of one or more pathways, or may share and/or part of more than one relevent biological phenomena. If we just go about documenting the effect caused by each gene individually to make the whole story, it will certainly be a cumbersome task if not impossible, and that too with a risk of concluding a contrasting/confusing results. To save us of this dilemma some smart people from [GO consortium](http://geneontology.org/)  came with a brilliant idea. The idea was to create a set of genes (a collection of genes) with a shared biological property, a comman pathways, and unsurprisingly a comman GO term etc. So instead of looking at the Diffrential Expression of individual genes between two groups (e.g. case and control) of samples, we look for the what GO Terms are enriched in our experiments. This gives us a broader view of the results at hand. Looking at the enriched GO terms we can make some conclusion about what pathway/process/components have been affected and draw a mechanistic perpective of the experiment. So GSEA essentially is a kind of statitical test (more aptly known as enrichment analyses) of functional terms that appear associated to the given set of differentially expressed genes more often than expected by chance. 

So lets begin by selecting genes that are significantly DE between CTRL and Case samples. Let’s extract genes that have an adjusted p-value below 0.1 and that show a 2-fold change (either negative or positive) in the case compared to CTRL.

```r
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
# eliminating all the genes for which p.values is NA
DE <- DEresults[!is.na(DEresults$padj),]

# Extracting all the genes having p.values < 0.1
DE <- DE[DE$padj < 0.1,]

# Extracting all the genes having a log fold change > 2
DE <- DE[abs(DE$log2FoldChange) > 1,]

# getting the set of genes 
genesOfInterest <- rownames(DE)

# we have 4367 DE genes
> length(genesOfInterest)
[1] 4367

library(DESeq2)
library(gprofiler2)
library(ggfortify)

# calculate enriched GO terms

go_terms <- gost(query = genesOfInterest,
                        organism = "hsapiens",
                        ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE,
                        exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, 
                        evcodes = TRUE, 
                        user_threshold = 0.05, 
                        correction_method = "g_SCS", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = 'GO', as_short_link = FALSE
                        )



publish_gosttable(go_terms, 
                  highlight_terms = head(go_terms$result[order(go_terms$result$p_value),],10),
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)

```
![image12](https://user-images.githubusercontent.com/85447250/217903209-00037714-95cd-4ff6-bd6b-7dcf5d54803e.png)

Fig. Table showing top 10 Go terms based on the p.values.

We also can help ourselves with some Manhattan style plot for the different kinds of annotations we can obtain with **gprofiler2** package. We can take full advantage of this package to find annotations from all other avilable sources such as KEGG, REAC, CORUM etc.

```r
# first we obtain annotations from all sources by setting sources=NULL in gost()
terms_all <- gost(query = genesOfInterest,
                 organism = "hsapiens",
                 ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE,
                 exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, 
                 evcodes = FALSE, 
                 user_threshold = 0.05, 
                 correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL, as_short_link = FALSE
)

> plot_object_static <- gostplot(terms_all, capped = TRUE, interactive = F)


> plot_obj_table <- publish_gostplot(plot_object_static, highlight_terms = c("GO:0015850", "REAC:R-HSA-1474244", "KEGG:04080"),
                                                                 width = NA, height = NA, filename = NA)

> plot_obj_table

```

![image14](https://user-images.githubusercontent.com/85447250/217913585-4f1094e9-ecc2-469f-8ed7-a859bfbc6463.png)

Fig. Annotations from different sources. Particularly, three from GO, KEGG, and REAC has been shown.


A similar aspect of GSEA is to do an enrichment analysis of a well defined gene sets (for example it may be subset of genes from count matrix with highest variations; and a gene set created with random sampling)  in CASE and CTRL and check if the well defined gene set is clearly enriched in CASE vs CTRL study as compared to the gene set created with random sampling. What essentially we are going to do is take two gene sets (as defined above) and do a group comparison between the case
samples with respect to the control samples (some thing like asking a question; is the well defined set of genes is clearly upregulated or downregulated when we look for the expression values in CASE vs CTRL samples? is the random set of genes shows any difference in expression values in CASE vs CTRL samples). We do it with the help of **gage** package.

We begin by getting those two sets of genes; 1) known genes set: we extract it from top GO terms from previous analysis (`go_terms`). 2) a gene set from random sampling of genes

```r
library(gage)

#significant GO terms found in the GO analysis. order go results by p.value
sig_go_terms <- go_terms$result[order(go_terms$result$p_value),]

#restrict the terms that have at most 100 genes overlapping with the query
go_top <- sig_go_terms[sig_go_terms$intersection_size < 100,]

go_top <- go_top[order(go_top$intersection),]

#use the top term from this table to create a gene set
geneSet1 <- unlist(strsplit(go_top[1,]$intersection, ','))

#define another gene set by just randomly selecting 25 genes from the counts tables
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)
geneSet2 <- sample(rownames(normalizedCounts), 25)

# a list of two gene sets
geneSets <- list('top_GO_term' = geneSet1,
                 'random_set' = geneSet2)


gseaResults <- gage(exprs = log2(normalizedCounts+1),
                    ref = match(rownames(colData[colData$group == 'CTRL',]),
                                colnames(normalizedCounts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalizedCounts)),
                    gsets = geneSets, compare = 'as.group')
                    
> gseaResults$greater
               p.geomean stat.mean        p.val        q.val set.size         exp1
top_GO_term 1.835290e-07 5.8490273 1.835290e-07 3.670581e-07       31 1.835290e-07
random_set  2.494493e-01 0.6821332 2.494493e-01 2.494493e-01       25 2.494493e-01

> gseaResults$less
            p.geomean stat.mean     p.val     q.val set.size      exp1
random_set  0.7505507 0.6821332 0.7505507 0.9999998       25 0.7505507
top_GO_term 0.9999998 5.8490273 0.9999998 0.9999998       31 0.9999998
```
We can see that the random gene set shows no signiﬁcant up- or down-regulation, while the gene set we deﬁned using the top GO term shows a signiﬁcant up-regulation (adjusted p-value < 1.835290e-07). At this point it is worthwhile to visualize these systematic changes in a heatmap.

```r
# get the expression data for the gene set of interest
M <- normalizedCounts[rownames(normalizedCounts) %in% geneSet1, ]

pheatmap(log2(M+1),
         annotation_col = colData,
         show_rownames = TRUE,
         fontsize_row = 8,
         scale = 'row',
         cutree_cols = 2,
         cutree_rows = 2)
```
![image15](https://user-images.githubusercontent.com/85447250/218023331-a51c3040-20b3-43c8-846e-b5c876070612.png)

Fig. Heatmap of expression value from the genes with the top GO term. 

We can see that almost all genes from this gene set display an increased level of expression in the case samples compared to the controls.




