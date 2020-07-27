#!/usr/bin/env nextflow
/*
========================================================================================
                         TOBIAS wrapper module
========================================================================================
 #### Homepage / Documentation

----------------------------------------------------------------------------------------
Module decription
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
*/

// Specify DSL2
nextflow.preview.dsl = 2

// dedup reusable component
process tobias_atacorrect {
    publishDir "${params.outdir}/tobias_atacorrect",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-tobias:latest'
    
	input:
		tuple val(sample_id), path(bam), path(genome), path(bed), path(blacklist)

	output:
		path "*.log", emit: report
		tuple val(sample_id), path("*_corrected.bw"), emit: corrected

    script:
    """
    TOBIAS ATACorrect --bam $bam --genome $genome --peaks $bed --blacklist $blacklist --outdir . --cores ${task.cpus} > ${sample_id}_atacorrect.log
    """
}

process tobias_footprint {
    publishDir "${params.outdir}/tobias_footprint",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-tobias:latest'

	input: tuple val(sample_id), path(corrected), path(bed)

	output:
		path "*.log", emit: report 
		tuple val(sample_id), path("*_footprints.bw"), emit: footprints
		val sample_id, emit: sampleIds
		path "*_footprints.bw", emit: footprintsDetect

	script:
	"""
	TOBIAS FootprintScores --signal $corrected --regions $bed --output ${sample_id}_footprints.bw --cores ${task.cpus} > ${sample_id}_footprint.log
	"""	
}

process tobias_bindetect {
	publishDir "${params.outdir}/tobias_bindetect",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-tobias:latest'

	input: tuple val(sample_ids), path(signals), path(motifs), path(genome), path(peaks)

	output:
		path "*.log", emit: report
		path "bindetect_*"
		path "*/*_overview.txt"
		//path "*/beds/**"
		path "*/plots/**"
		path "*/beds/*_all.bed", emit: motifBeds
		path "*/beds/*_bound.bed", emit: boundBeds
		path "*/beds/*_unbound.bed", emit: unboundBeds
		path "motiflist.txt", emit: motifList
		//tuple val(sample_ids), path(signals), path(motifs), path(genome), path(peaks), emit: bindetect
//${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${PEAK_TYPE}","")}
	script:
	"""
	awk '{if(\$1 ~ />/) {a=\$2"_"\$1; gsub(/>/,"",a); gsub("::","",a); gsub(/[()]/,"",a); print a}}' $motifs | sort > motiflist.txt 
	TOBIAS BINDetect --motifs $motifs --signal $signals --peaks $peaks --genome $genome --outdir . --cond_names ${sample_ids.join(' ')} --cores ${task.cpus} > bindectect.log
	"""	
}
process tobias_plotaggregate {
	publishDir "${params.outdir}/tobias_bindetect",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-tobias:latest'

	input: tuple val(motifname), path(motifBeds), path(correctedInsertions)

	output:
		path "$motifname/*.pdf", emit: plotaggregated

	script:
	"""
	TOBIAS PlotAggregate --TFBS $motifBeds --signals ${correctedInsertions.join(' ')} --output $motifname/${motifname}_plotaggregate.pdf --share_y both --plot_boundaries --signal-on-x
	"""	
}

//TOBIAS PlotAggregate --TFBS $motifBeds  --signals $correctedInsertions --output  --share_y both --plot_boundaries --signal-on-x

// TOBIAS PlotAggregate --TFBS BINDetect_output_D6_0_1/OLIG2_MA0678.1/beds/OLIG2_MA0678.1_all.bed \
// 	--signals ATACorrect_mergedReps/D6_0_1.mRp.clN.sorted_corrected.bw \
// 		ATACorrect_mergedReps/D6_10_2.mRp.clN.sorted_corrected.bw \
// 		ATACorrect_mergedReps/D6_100_M.mRp.clN.sorted_corrected.bw \
// 		ATACorrect_mergedReps/D6_500_3.mRp.clN.sorted_corrected.bw \
// 	--output OLIG2_footprint_comparison_D6.pdf \
// 	--share_y both --plot_boundaries --signal-on-x

	/*if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		cores = 4

		readString = ""

		// Options we add are
		bowtie2_options = bowtie2_args
		bowtie2_options +=  " --no-unal "  // We don't need unaligned reads in the BAM file
		
		// single-end / paired-end distinction. Might also be handled via params.single_end 
		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
		}
		else {
			readString = "-U " + reads
		}*/
		