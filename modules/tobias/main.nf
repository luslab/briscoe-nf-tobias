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

// Local default params
params.internal_outdir = params.outdir
params.internal_process_name = 'tobias'

// dedup reusable component
process tobiasatacorrect {
    publishDir "${params.internal_outdir}/${params.internal_process_name}",
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

process tobiasfootprint {
    publishDir "${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-tobias:latest'

	input: tuple val(sample_id), path(corrected), path(bed)

	output:
		path "*.log", emit: report // can these have the same name as above and also *.log only?
		tuple val(sample_id), path("*_footprints.bw"), emit: footprints

	script:
	"""
	TOBIAS FootprintScores --signal $corrected --regions $bed --output ${sample_id}_footprints.bw --cores ${task.cpus} > ${sample_id}_footprint.log
	"""	
}


//
// TOBIAS FootprintScores --signal $i \
//	--regions $BED  \
//	--output ATACFootprint_mergedReps/${cleanName}_footprints.bw --cores 16
//TOBIAS ATACorrect --bam $i --genome $GENOME --peaks $BED --blacklist $BLACKLIST --outdir ATACorrect_mergedReps --cores 16 

// D3_0_NMP.mRp.clN.sorted_AtacBias.pickle
// -rw-r--r-- 1 delasj domain_users  34K Jun 11 23:42 D3_0_NMP.mRp.clN.sorted_atacorrect.pdf
// -rw-r--r-- 1 delasj domain_users 1.4G Jun 11 23:40 D3_0_NMP.mRp.clN.sorted_bias.bw
// -rw-r--r-- 1 delasj domain_users 806M Jun 11 23:42 D3_0_NMP.mRp.clN.sorted_corrected.bw
// -rw-r--r-- 1 delasj domain_users 771M Jun 11 23:41 D3_0_NMP.mRp.clN.sorted_expected.bw
// -rw-r--r-- 1 delasj domain_users 167M Jun 11 23:38 D3_0_NMP.mRp.clN.sorted_uncorrected.bw

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
		