# briscoe-nf-tobias
Analysis pipeline for Joaquinas work with the tobias analysis tool

## Pipeline outline

The pipeline runs the different steps of the [TOBIAS footprinting tool](https://github.com/loosolab/TOBIAS). For detailed descriptions of each step visit their [wiki](https://github.com/loosolab/TOBIAS/wiki).

Given a design file, a motif file and mapped files (bam format) it will run TOBIAS' ATACorrect, ScoreBigwig, BINDetect of all against all conditions and it will generate PlotAggregate metaplots of corrected insertions for each motif in the file.


## Usage

The release is set to the latest version, to run the development version of the pipeline, change to `-r dev`

The pipeline can be run as in the following example:

```

#!/bin/sh

export NXF_WORK="/path/to/work"
export NXF_SINGULARITY_CACHEDIR="/path/to/sing"

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.04.0
ml Singularity/3.4.2
ml Graphviz

## UPDATE PIPLINE
nextflow pull luslab/briscoe-nf-tobias

## RUN PIPLINE
nextflow run luslab/briscoe-nf-tobias \
  -r master \
  -profile crick_tobias \
  --genome genome.fa \
  --regions regions.bed \
  --peaks regions.bed \
  --blacklist blacklist.bed \
  --motifs motifs \
  --design design.csv \
  --skip_bam_index true
  
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

### --skip_bam_index
This should be `true` if you are providing paths to your .bam.bai

### --output_peaks <BED> (optional)
Provide this parameter as a bed file. From the manual:
  ```
    --output-peaks <bed>             Gives the possibility to set the output peak set
                                   differently than the input --peaks. This will limit all
                                   analysis to the regions in --output-peaks. NOTE:
                                   --peaks must still be set to the full peak set!
  ```
