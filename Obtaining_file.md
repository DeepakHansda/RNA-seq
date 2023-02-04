First we need to get the data. For this purpose we can extract data from SRA SRP029880 (Kim et al 2016). The objective of this study is to identify a prognostic signature in colorectal cancer (CRC) patients with diverse progression and heterogeneity of CRCs. Authors of this study have tried to identify upregulated genes associated with colorectal cancer (CRC) liver metastasis (CLM). They have conducted 54 Rna-seq samples from 18 CRC patients. For the the purpose of this workflow we will take only 10 samples (5 metastasized colorectal cancer samples and 5 normal colon samples). count data `counts` will be in the form of a matrix with rows representing genes and columns representing samples (5 CRC and 5 normal colon, with one column `width` indicating the genomic length of each gene) 


```r
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

  
  
  

  
  
  
  





