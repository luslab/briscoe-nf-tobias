# briscoe-nf-tobias
Analysis pipeline for Joaquinas work with the tobias analysis tool

## Pipeline outline

The pipeline runs the different steps of the [TOBIAS footprinting tool](https://github.com/loosolab/TOBIAS). For detailed descriptions of each step visit their [wiki](https://github.com/loosolab/TOBIAS/wiki).

Given a design file, a motif file and mapped files (bam format) it will run TOBIAS' ATACorrect, ScoreBigwig, BINDetect of all against all conditions and it will generate PlotAggregate metaplots of corrected insertions for each motif in the file.


## Usage

The pipeline can be run in the following example:

```
nextflow run luslab/briscoe-nf-tobias \
  -r master \
  -profile crick \
  --genome genome.fa \
  --regions regions.bed \
  --peaks regions.bed \
  --blacklist blacklist.bed \
  --motifs motifs \
  --design design.csv

  ```

### --genome
Path to genome in fasta format used in ATACorrect and BINDetect.

### --regions
Path to BED file with the open chromatin regions for these experiment. Used for ATACorrect.

### --peaks 
Can be the same BED file as `--regions` or a different BED file of regions of interest. Used for BINDetect.

### --blacklist
Path to Blacklisted regions in BED format. Used for ATACorrect

### --motifs
Path to Motif file in JASPAR or MEME format.

### --design
Csv file with paths and condition grouping with the following header:
`sample_id,data1,data2` 

 - `sample_id` correspond to the sample name
 - `data1` corresponds to the path to the .bam
 - `data2` corresponds to the path to the .bam.bai, if available

