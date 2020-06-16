#!/bin/bash
GENOME=/camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_Full/results/reference_genome/genome.fa
BED=/camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_Full/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed
BLACKLIST=/camp/home/delasj/.nextflow/assets/nf-core/atacseq/assets/blacklists/mm10-blacklist.bed

echo "#!/bin/bash" > go.sh

for i in /camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_FullFiltered/results/bwa/mergedReplicate/*.bam; do

cleanName=${i%%.mRp.clN.sorted.bam}
cleanName=${cleanName##*/}

echo $cleanName

echo "#!/bin/bash" > RunATACorrect_$cleanName.sh;
echo "TOBIAS ATACorrect --bam $i --genome $GENOME --peaks $BED --blacklist $BLACKLIST --outdir ATACorrect_mergedReps --cores 16 " >> RunATACorrect_$cleanName.sh;

echo "sbatch --time=10:00:0 --cpus-per-task=16 --job-name=ACorrect_${cleanName} --mail-type=END,FAIL --mail-user=joaquina.delas@crick.ac.uk \
    --output=ACorr_${cleanName}.o --error=ACorr_${cleanName}.e RunATACorrect_$cleanName.sh" >> go.sh;

done

echo "Done: run 'go.sh' to submit jobs"
