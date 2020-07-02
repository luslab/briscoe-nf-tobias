#!/usr/bin/env nextflow

/*
========================================================================================
                         Briscoe NF Tobias Wrapper
========================================================================================
 #### Homepage / Documentation
 https://github.com/luslab/briscoe-nf-tobias
----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.preview.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Params
-------------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { tobiasatacorrect; tobiasfootprint } from './modules/tobias/main.nf'
include { build_debug_param_summary; luslab_header; check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'
include simple_metadata from './luslab-nf-modules/tools/metadata/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/

// Main data parameters
//params.input = ''
//params.bowtie_index = ''

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/

// Show banner
log.info luslab_header()

if(params.verbose){
    log.info build_debug_param_summary()
}

// Show work summary
//def summary = [:]
/*summary['Classpath'] = params.classpath
summary['Output directory'] = params.outdir
summary['Trace directory'] = params.tracedir
summary['Max CPUs'] = params.max_cpus
summary['Max memory'] = params.max_memory
summary['Max time'] = params.max_time
summary['Bowtie index path'] = params.bowtie_index
summary['Star index path'] = params.star_index
summary['Genome path'] = params.genome
summary['Genome index path'] = params.genome_fai
summary['Segmentation path'] = params.segmentation
summary['Regions path'] = params.peka_regions
summary['Metadata path'] = params.input
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"*/

check_params(['genome','regions','blacklist','design'])

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

ch_genome = Channel.fromPath(params.genome, checkIfExists: true)
ch_regions = Channel.fromPath(params.regions, checkIfExists: true)
ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true)

// Run workflow
workflow {
    simple_metadata(params.design)

    ch_atacorrectinput = simple_metadata.out.metadata
        .map{row -> [row[0], file(row[1][0])]}
        .combine(ch_genome)
        .combine(ch_regions)
        .combine(ch_blacklist)

    tobiasatacorrect( ch_atacorrectinput )

    tobiasfootprint( tobiasatacorrect.out.corrected.combine(ch_regions) )
}

// workflow.onComplete {
//     log.info "\nPipeline finished executing\n"
// }

/*------------------------------------------------------------------------------------*/