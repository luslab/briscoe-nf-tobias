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

include tobias from './modules/tobias/main.nf'

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
//log.info luslabHeader()

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

// Check params
/*if (!params.input) {
    exit 1, "No metadata input provided"
}
if (!params.bowtie_index) {
    exit 1, "No bowtie index provided"
}
if (!params.star_index) {
    exit 1, "No star index provided"
}
if (!params.genome) {
    exit 1, "No genome provided"
}
if (!params.genome_fai) {
    exit 1, "No genome index provided"
}
if (!params.segmentation) {
    exit 1, "No segmentation provided"
}
if (!params.peka_regions) {
    exit 1, "No regions provided"
}*/

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    tobias()
}

workflow.onComplete {
    log.info "\nPipeline finished executing\n"
}

/*------------------------------------------------------------------------------------*/