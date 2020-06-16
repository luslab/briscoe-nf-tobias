#!/bin/bash
TOBIAS ATACorrect --bam /camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_FullFiltered/results/bwa/mergedReplicate/D3_0_NMP.mRp.clN.sorted.bam --genome /camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_Full/results/reference_genome/genome.fa --peaks /camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_Full/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed --blacklist /camp/home/delasj/.nextflow/assets/nf-core/atacseq/assets/blacklists/mm10-blacklist.bed --outdir ATACorrect_mergedReps --cores 16 
