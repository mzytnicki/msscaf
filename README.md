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
ontData     <- parsePafFile(ontFileName, resolution)
bamFileName <- "/path/to/long/10X.bam"
bamData     <- parseBamFile(bamFileName, resolution)
```

Load data:
```
contigFileName <- "/path/to/contigs.fa"
scaffolder     <- msscaf(contigFileName, resolution)
scaffolder     <- addExp(scaffolder, ontData, "ONT")
scaffolder     <- addExp(scaffolder, hicData, "HiC")
scaffolder     <- addExp(scaffolder, bamData, "10X")
```

Estimate distributions:
```
scaffolder <- estimateDistributions(scaffolder)
scaffolder <- cleanData(scaffolder)
```

Find breaks:
```
pvalue     <- 0.01
scaffolder <- findBreaks(scaffolder, pvalue)
scaffolder <- splitChromosomes(scaffolder)
```

Find merges:
```
scaffolder <- findJoins(scaffolder, pvalue)
scaffolder <- scaffold(scaffolder)
```

Write scaffolds to file:
```
outputFileName <- "/path/to/output/scaffolds.fa"
Biostrings::writeXStringSet(Biostrings::DNAStringSet(scaffolder@sequences), outputFileName)
```
