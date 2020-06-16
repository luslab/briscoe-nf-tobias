
#!/bin/bash
GENOME=/camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_Full/results/reference_genome/genome.fa
MOTIF=/camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_2conditions_TOBIAS/database_motifs/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt
BED=/camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_Full/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed

echo "#!/bin/bash" > go.sh

for i in /camp/lab/briscoej/working/Joaquina/hpc_camp/TOBIAS_ATAC_Filtered2/ATACFootprint_mergedReps/*.bw; do

cleanName=${i%%_footprints.bw}
cleanName=${cleanName##*/}

echo $cleanName

echo "#!/bin/bash" > RunBINDetect_$cleanName.sh;
echo "
QT_QPA_PLATFORM=offscreen TOBIAS BINDetect --motifs $MOTIF \
	--signals $i \
	--genome $GENOME \
	--peaks $BED \
    --outdir BINDetect_output_${cleanName} \
	--cond_names $cleanName \
	--cores 8 " >> RunBINDetect_$cleanName.sh;

echo "sbatch --time=10:00:0 --cpus-per-task=16 --job-name=ABDetect_${cleanName} --mail-type=END,FAIL --mail-user=joaquina.delas@crick.ac.uk \
    --output=ABDetect_${cleanName}.o --error=ABDetect_${cleanName}.e RunBINDetect_$cleanName.sh" >> go.sh;

done

echo "Done: run 'go.sh' to submit jobs"

