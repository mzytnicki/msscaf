# msscaf
Multiple-source scaffolder

This tool scaffold contigs using several sources.
To date, msscaf supports:
 - long reads (mapped to the contigs, in `.paf` format, for instance mapped with `minimap2`);
 - 10X (mapped with LongRanger);
 - Hi-C (in `.hic` format).

msscaf is written in R.
No previous knowledge of R is required.
You can adapt the following pipe-line to suit your needs.

Simply get R installed, and start it with the command
```
 R
```

## Pipe-line

Download package:
```
if (!("devtools" %in% rownames(installed.packages()))) install.packages("devtools")
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
# Alternatively, the "parseCoolFile" is available for the .cool/.mcool formats
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


## Accessing data

The breaks can be accessed with:
```
scaffolder@breaks
```
The joins can be accessed with:
```
scaffolder@joins
```


## Visualization

The whole matrices can be displayed with:
```
plot.msscaf(scaffolder)
```

One contig/scaffold (here, `ctg1`) can be displayed with:
```
plot.msscaf(scaffolder, ref1 = ctg1)
```

The interaction between two contig/scaffolds (here, `ctg1` and `ctg2`) can be displayed with:
```
plot.msscaf(scaffolder, ref1 = ctg1, ref2 = ctg2)
```

One break (located at position `ctg1:100`) can be displayed with:
```
plot.msscafBreak(scaffolder, ref = ctg1, bin = 100)
```

One join (joining the right ends of `ctg1` and `ctg2`) can be displayed with:
```
plot.msscafJoin(scaffolder, ref1 = ctg1, ref2 = ctg2, after1 = TRUE, after2 = TRUE)
```
