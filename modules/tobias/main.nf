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
process tobias {
    publishDir "${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-tobias:latest'
    //input:
    //  tuple val(sample_id), path(bai), path(bam)
       
   // output:
   //   tuple val(sample_id), path(bam), emit: dedupBam

    shell:
    """
    TOBIAS --version
    """
}

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
		