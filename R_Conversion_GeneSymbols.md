# Conversion of Entrez IDs to Gene Symbols.

So we have columns of data with Entrez IDs but no corresponding Gene Symbols.

> How on Earth are we going to solve it ?

The solution is simple , anyone who is familiar with R and a bit of bioinformatics can do it using a few steps. 

## STEP - 1 :  Install the prerequisites.

We are going to use BiocManager and a few packages that I'll mention in the code snippet below. 

```r
# This package helps you to get Gene Symbols from ENSEMBL IDs. 
BiocManager::install("EnsDb.Hsapiens.v79")

# This package helps you get Gene Symbols from Entrez IDs.
BiocManager::install('org.Hs.eg.db')

#Helps with Annotation
BiocManager::install("annotate")
```

## STEP - 2 : Get started

Now we need the genes that we need to convert and import the libraries we need. 

```r
library(org.Hs.eg.db)
library(annotate)

IDs <- c('3815', '3816', '2341')
```

## STEP - 3 : Convert

Now just rip the bandage off and write a line of code.

```r
getSYMBOL(IDs, data='org.Hs.eg')
```

This script can be used with ENSEMBL as well.

## STEP - 4 : (Optional) Reverse Conversion

This will give out ENTREZ IDs from Gene Symbols.

```r
library('org.Hs.eg.db')

# you will have your own list here
symbols <- c('AHNAK', 'BOD1L1', 'HSPB1', 'SMARCA4', 'TRIM28')

# use mapIds method to obtain Entrez IDs
mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
```
