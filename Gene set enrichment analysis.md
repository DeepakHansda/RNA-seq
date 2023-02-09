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
                        evcodes = FALSE, 
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

Fig. Annotations from different sources. Particularly, three from GO, KEGG, and REAC has been dipicted.




