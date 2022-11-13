# Finding Common Gene Symbols from two or more different sources. 

## Step 1 

First we must have the genes we need : 

```r
GeneList1 <- c("G1" , "G2" , "G3" , "G4" , "G5")
GeneList2 <- c("G3" , "G2" , "G6" ,"G8" , "G5")
```

> These are fictional genes for the sake of the tutorial , do not use them in your workplace.

## Step 2

We use mathematics to retrieve the common genes.

This is called intersection in Set Theory , represented with an inverted U.

Lets say :

Set A = (1,2,3,4,5)
Set B = (3,4,5,6,7)

Set C is Set A intersection with Set B

Set C will be common elements in the set A and set B

Hence , Set C will be (3,4,5)

```r
CommonGenes <- intersect(GeneList1,GeneList2)
```

## Step 3 (Optional)

Lets try it out with 3 Lists. 

```r
GeneList1 <- c("G1", "G2" , "G3" , "G4")
GeneList2 <- c("G2", "G3" , "G4" , "G5")
GeneList3 <- c("G3", "G4" , "G5" , "G6")
```
So the intersection among these three lists would be :

```r
intersection(intersection(GeneList1, GeneList2), GeneList3)
```