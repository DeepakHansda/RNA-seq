First we need to get the data. For this purpose we can extract data from SRA SRP029880 (Kim et al 2016). The objective of this study is to identify a prognostic signature in colorectal cancer (CRC) patients with diverse progression and heterogeneity of CRCs. Authors of this study have tried to identify upregulated genes associated with colorectal cancer (CRC) liver metastasis (CLM). They have conducted 54 Rna-seq samples from 18 CRC patients. For the the purpose of this workflow we will take only 10 samples (5 metastasized colorectal cancer samples and 5 normal colon samples). count data `counts` will be in the form of a matrix with rows representing genes and columns representing samples (5 CRC and 5 normal colon, with one column `width` indicating the genomic length of each gene) 


```r

library(pheatmap)
library(stats)
library(ggplot2)
library(ggfortify)

counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))

> dim(counts)
[1] 19719    11

> counts[1:3,]
       CASE_1 CASE_2 CASE_3 CASE_4 CASE_5 CTRL_1 CTRL_2 CTRL_3 CTRL_4 CTRL_5 width
TSPAN6 776426 371725 612244 456147 513335 559544 489653 332084 238516 634115 12883
TNMD     1483    806   2995    297   1095   4631   1884   4484   1961   3976 15084
DPM1   364919 274342 248740 371045 325628 211173 123204 113606  67338 198331 23689

# Computing CPM

cpm <- apply(subset(counts, select = c(-width)), 2,
             function(x) x/sum(as.numeric(x)) * 10^6)

# as per the definition of cpm, column should add up to a million
> colSums(cpm)
CASE_1 CASE_2 CASE_3 CASE_4 CASE_5 CTRL_1 CTRL_2 CTRL_3 CTRL_4 CTRL_5 
 1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06 

# Computing  RPKM            

# First we need to creat a vector for gene lengts
geneLengths <- as.vector(subset(counts, select = c(width)))

rpkm <- apply(subset(counts, select = c(-width)), 2,
              function(x) {10^9 * x / geneLengths / sum(as.numeric(x))
              })

# Computing TPM

# First we need to find the normalized values of read counts normalized by the gene lengths (geneLengths/1000; i.e., # of 
# reads per base kilo base pair)

rpk <- apply( subset(counts, select = c(-width)), 2,
              function(x) x/(geneLengths/1000))

# Now we can normalize rpk by the sample size and figure out how many transcripts are there out a million reads

tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

# Now again, as per the definition of tpm, the sample size should add up to a million which indeed is the case.

> colSums(tpm)
CASE_1 CASE_2 CASE_3 CASE_4 CASE_5 CTRL_1 CTRL_2 CTRL_3 CTRL_4 CTRL_5 
 1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06  1e+06 
```

### Explorartory analysis of read count table
+ clustering
  
  It is generally not a good idea to do clustering using all the available genes. It can take lot of time and resources, so instead we do clustering with the help of 100 most variable genes in `tpm` data.
  
  ```r
  # Computing the variance of genes across samples in tpm
  V <- apply(tpm, 1, var)
  
  # sort the V by decreasing order of variance and select top 100 genes.
  
  selectedGenes <- names(V[order(V, decreasing = T)][1:100])
  
  library(pheatmap)
  
  pheatmap(tpm[selectedGenes,], scale = "row", show_rownames = FALSE)
  ```
  
 
 ![image1](https://user-images.githubusercontent.com/85447250/216706661-1fe21cfc-173e-4c6b-bbb3-1e46766a6154.png)
  
  Fig. Clustering and visualization of the topmost variable genes as a heatmap. Columns shows different samples while rows shows 100 topmost variable genes.
  
  We can put some annotations on the heatmap above, primary purpose being to check if the replicates of a sample cluster together or not. Ideally, replicates should cluster together.
  
  ```r
  coldata <- read.table(coldata_file, header = T, sep = '\t',
                      stringsAsFactors = TRUE)
                      
      pheatmap(tpm[selectedGenes,], scale = "row", show_rownames = FALSE, 
         annotation_col = coldata)
 ```
 
 
 
 ![image2](https://user-images.githubusercontent.com/85447250/216751591-c7ee930e-abef-4032-9cad-eeacdb7fa187.png)
 
 Fig. Clustering and visualization of the topmost variable genes as a heatmap with annotation.
 
 
 + Dimensionality reduction: PCA
  
  Apart from heatmap we can also do PCA to confirm if the replicates of a type cluster together. Here too, we should see a clear separation between `CASE` and `CTRL` samples.   
 
 ```r
 # we need to transpose the count matrix, because that is the way a generic PCA calculation function works
 M <- t(tpm[selectedGenes,])
 
 # Making sure the the log function does not have to evaluate a $0$
 
M <- log2(M+1)

pcaResults <- prcomp(M)

pca.plot <- autoplot(pcaResults,
                          data = coldata,
                          colour = 'group')

pca.plot
 ```
![image3](https://user-images.githubusercontent.com/85447250/216831642-37265991-a6d0-433c-9076-52dd1311770c.png)

Fig. PCA plot of CASE and CTRL data using two largest PCs. We can see a clear separation between them.

+ Correlation plot

Another similar approach to look for the similarity between replicates is to compute the correlation between samples.

```r
library(stats)

correlationMatrix <- cor(tpm)
library(pheatmap)
pheatmap(correlationMatrix, annotation_col = coldata)
```

![image4](https://user-images.githubusercontent.com/85447250/216833751-7e383b90-ce01-4fac-9215-4b2701d031d6.png)

Fig. Pairwise correlation of samples displayed as heatmap


### Differential gene expression analysis 

Differential gene analysis compares a gene to thousands of other genes with the null hypothesis that the expression level of gene is same in two different samples. 

```r

library(DEseq2)
# to do a differential expression analysis using 'DESeq2' we need $3$ input objects, the count matrix, metadata of count 
# matrix, and a design formula

countData <- as.matrix(subset(counts, select = c(-width)))

colData <- read.table(coldata_file, header = TRUE, sep = "\t",
                      stringsAsFactors = TRUE) 
designformula <- "~ group"

# Creating a DEseqDataSet object 
dds <- DESeqDataSetFromMatrix(countData, colData = colData, 
                              design = as.formula(designformula) )

# we are taking out all such genes which does not have even have total sum of 1 across all samples (10 samples)

dds <- dds[rowSums(DESeq2::counts(dds)) > 1,]

# Now we can use DEseq() function from DESeq2 which will normalize the counts, estimate the dispersion values, and 
# compute the generalized linear model with the provided design formula. It does return an updated dds object.

dds <- DESeq(dds)
```

Now we can contrast and compare the sample based upon some variable. For instance, we have here "group" variable which describe
the `case` or `control` label for the samples. 

```r
DEresults = results(dds, contrast = c("group", "CASE", "CTRL"))

# arranging the result with decreasing order of p.values

DEresults <- DEresults[order(DEresults$pvalue),]
```

Above we have obtained a table containing the differential expression status of `case` samples compared to the `control` samples. It is important to note that the sequence of the elements provided in the contrast argument determines which group of samples are to be used as the `control`. This impacts the way the results are interpreted, for instance, if a gene is found up-
regulated (has a positive log2 fold change), the up-regulation status is only relative to the factor that is provided as `control`. In this case, we used samples from the “CTRL” group as `control` and contrasted the samples from the “CASE” group with respect to the “CTRL” samples. Thus genes with a positive log2 fold change are called up-regulated in the case samples with respect to the `control`, while genes with a negative log2 fold change are down-regulated in the `case` samples. Whether the deregulation is signiﬁcant or not, warrants assessment of the adjusted p-values. We can have a look at DEresult below.

```r
> DEresults
log2 fold change (MLE): group CASE vs CTRL 
Wald test p-value: group CASE vs CTRL 
DataFrame with 19097 rows and 6 columns
            baseMean log2FoldChange     lfcSE       stat       pvalue         padj
           <numeric>      <numeric> <numeric>  <numeric>    <numeric>    <numeric>
CYP2E1       4829889        9.36024  0.215223    43.4909  0.00000e+00  0.00000e+00
FCGBP       10349993       -7.57579  0.186433   -40.6355  0.00000e+00  0.00000e+00
ASGR2         426422        8.01830  0.216207    37.0863 4.67898e-301 2.87741e-297
GCKR          100183        7.82841  0.233376    33.5442 1.09479e-246 5.04945e-243
APOA5         438054       10.20248  0.312503    32.6477 8.64906e-234 3.19133e-230
...              ...            ...       ...        ...          ...          ...
CCDC195      20.4981      -0.215607   2.89255 -0.0745386           NA           NA
SPEM3        23.6370     -22.154765   3.02785 -7.3170030           NA           NA
AC022167.5   21.8451      -2.056240   2.89545 -0.7101618           NA           NA
BX276092.9   29.9636       0.407326   2.89048  0.1409199           NA           NA
ETDC         22.5675      -1.795274   2.89421 -0.6202983           NA           NA
> 
```

The ﬁrst three lines in this output show the contrast and the statistical test that were used to compute these results. Below these lines is the actual table with 6 columns:

+ **baseMean** represents the average normalized expression of the gene across all considered samples.
+ **log2FoldChange** represents the base-2 logarithm of the fold change of the normalized expression of the gene in the given contrast.
+ **lfcSE** represents the standard error of log2 fold change estimate.
+ **stat** is the statistic calculated in the contrast.
+ **pvalue** represent the pvalue.
+ **padj** represents the pvalue adjusted for multiple testing.


### Some more diagnostic plots

At this point, before going to further downstream analysis, it is important to look for the quality of our data in hand and to check if that has improved or not.

+ MA plot:

An MA plot is a useful tool to check if the data normalization worked well. The MA plot is a scatter plot where the x-axis denotes the average of normalized counts across samples and the y-axis denotes the log fold change in the given
contrast (here case vs control). Most points are expected to be on the horizontal 0 line (most genes are not expected to be differentially expressed).

```r
DESeq2::plotMA(object = dds, ylim = c(-5, 5))
```

![image5](https://user-images.githubusercontent.com/85447250/216890585-78fc1103-8854-46c5-b060-1575bf820b96.png)

Fig. MA plot of differential expression results.

+ p-value distribution

It is also important to observe the distribution of raw p-values (Figure 8.7). We expect to see a peak around low p-values and a uniform distribution at P-values above 0.1. Otherwise, adjustment for multiple testing would not work and the results
are not meaningful.

```r
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) +
  geom_histogram(bins = 100)
```

![image6](https://user-images.githubusercontent.com/85447250/216891706-db24f9b4-a9d8-4509-bcf7-77af4d14c1dc.png)

Fig. P-value distribution genes before adjusting for multiple testing











  
  
  

  
  
  
  





