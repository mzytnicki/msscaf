# msscaf
Multiple-source scaffolder

This tool scaffold contigs using several sources.
To date, msscaf supports:
 - long reads (mapped to the contigs, in `.paf` format, for instance mapped with `minimap2`);
 - 10X (mapped with LongRanger);
 - Hi-C (in `.hic` format).

Download package:
```
devtools::install_github("mzytnicki/msscaf")
```

Load package:
```
library("msscaf")
```

Read data:
```
resolution  <- 10000
hicFileName <- "/path/to/hic/file.hic"
hicData     <- parseHicFile(hicFileName, resolution)
ontFileName <- "/path/to/long/reads.paf"
ontData     <- parsePafFile(ontFileName, resolution, 500, 10, 10)
bamFileName <- "/path/to/long/10X.bam"
bamData     <- parseBamFile(bamFileName, resolution)
allObject   <- tenxchecker(resolution)
allObject   <- addExp(allObject, ontData, "ONT")
allObject   <- addExp(allObject, hicData, "HiC")
allObject   <- addExp(allObject, bamData, "10X")
```

Estimate distributions:
```
allObject <- estimateDistributions(allObject)
allObject <- cleanData(allObject)
```

Find breaks:
```
pvalue    <- 0.01
allObject <- findBreaks(allObject, pvalue)
allObject <- splitChromosomes(allObject)
```

Find merges:
```
allObject <- findJoins(allObject, pvalue)
allObject <- scaffold(allObject)
```

Write scaffolds to file:
```
outputFileName <- "/path/to/output/scaffolds.fa"
Biostrings::writeXStringSet(Biostrings::DNAStringSet(allObject@sequences), outputFileName)
```
