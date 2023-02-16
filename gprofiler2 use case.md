
### Functional enrichment analysis

 We will use `gost` from **gprofiler2** package to do gene set enrichment analysis. 

```r
term_results <- gost(query = genesOfInterest,
                        organism = "hsapiens",
                        ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE,
                        exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, 
                        evcodes = FALSE, 
                        user_threshold = 0.05, 
                        correction_method = "g_SCS", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
                        
 # gost() produces a named list; the list includes a data.frame (result) and a metadata object (meta)
 
 > names(term_results)
[1] "result" "meta" 

> dim(term_results$result)
[1] 1688   14

> head(term_results$result,3)
    query significant       p_value term_size query_size intersection_size
1 query_1        TRUE  1.471339e-03         6        526                 6
2 query_1        TRUE  4.991540e-02         6        526                 5
3 query_1        TRUE 9.576567e-111      7463       3790              1957
    precision    recall    term_id source
1 0.011406844 1.0000000 CORUM:7268  CORUM
2 0.009505703 0.8333333 CORUM:7265  CORUM
3 0.516358839 0.2622270 GO:0032501  GO:BP
                                              term_name effective_domain_size
1    HEPACAM-MLC1-Na,K-ATPase-Kir4.1-AQP4-TRPV4 complex                  3385
2 MLC1-Na,K-ATPase-Kir4.1-AQP4-TRPV4-syntrophin complex                  3385
3                      multicellular organismal process                 21092
  source_order       parents
1         2733 CORUM:0000000
2         2730 CORUM:0000000
3         8794    GO:0008150

```

Once we are done with the gene set enrichment analysis, we can visualize what GO, REAC, KEGG are enriched in our gene sets.

```r
plot_object <- gostplot(term_results, capped = TRUE, interactive = TRUE)

plot_obj_table <- publish_gostplot(plot_object, highlight_terms = c("GO:0015850", 
                                                        "REAC:R-HSA-1474244",
                                                        "KEGG:04080"), 
                       width = NA, height = NA, filename = NULL )
plot_obj_table

```
![image9](https://user-images.githubusercontent.com/85447250/217650364-f34224e0-aa58-405e-a2ac-e0ed73cc9d31.png)

Fig. GSEA analysis showing some enriched GO, KEGG, and REAC pathways

We can also get the enriched terms in our data sets as a table. Table below shows top ten terms based on p.values

```r
publish_gosttable(term_results, highlight_terms = head(term_results$result[order(term_results$result$p_value),],10),
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)
                  
```

![image11](https://user-images.githubusercontent.com/85447250/217881068-180c7be1-37f0-472f-9cf4-b58f4149e751.png)

Fig. Table of enriched terms. For representation purpose only top ten are being shown here.
