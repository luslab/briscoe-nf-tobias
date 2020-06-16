#!/bin/bash
BED=/camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_Full/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed

echo "#!/bin/bash" > go.sh

for i in /camp/lab/briscoej/working/Joaquina/hpc_camp/TOBIAS_ATAC_Filtered2/ATACorrect_mergedReps/*_corrected.bw; do

cleanName=${i%%.mRp.clN.sorted_corrected.bw}
cleanName=${cleanName##*/}

echo $cleanName

echo "#!/bin/bash" > RunFootprint_$cleanName.sh;
echo "TOBIAS FootprintScores --signal $i \
	--regions $BED  \
	--output ATACFootprint_mergedReps/${cleanName}_footprints.bw --cores 16" >> RunFootprint_$cleanName.sh;

echo "sbatch --time=10:00:0 --cpus-per-task=16 --job-name=AFootprint_${cleanName} --mail-type=END,FAIL --mail-user=joaquina.delas@crick.ac.uk \
    --output=AFootprint_${cleanName}.o --error=AFootprint_${cleanName}.e RunFootprint_$cleanName.sh" >> go.sh;

done

echo "Done: run 'go.sh' to submit jobs"
