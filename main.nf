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
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Params
-------------------------------------------------------------------------------------------------------------------------------*/

params.verbose = true

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { tobias_atacorrect; tobias_footprint; tobias_bindetect; tobias_plotaggregate } from './modules/tobias/main.nf'
include { build_debug_param_summary; luslab_header; check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'
include { fastq_metadata } from './luslab-nf-modules/tools/metadata/main.nf'
include { awk } from './luslab-nf-modules/tools/luslab_linux_tools/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/

// Main data parameters
//params.input = ''
//params.bowtie_index = ''

// 'awk' {
//             args             = ""
//             outfile_name     = ""
//             publish_dir      = "awk"
//             publish_results  = "all"
//         }

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

//check_params(['genome','regions','blacklist','design','motifs','peaks'])

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

ch_genome = Channel.fromPath(params.genome, checkIfExists: true)
ch_regions = Channel.fromPath(params.regions, checkIfExists: true)
ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true)
ch_peaks = Channel.fromPath(params.peaks, checkIfExists: false)

motifs_format = [
      [[:], params.motifs]
]

Channel
    .from(motifs_format)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_motifs}

workflow {
    fastq_metadata(params.design)

    awk(params.modules['motifsplit_awk'], ch_motifs)

    tobias_atacorrect( fastq_metadata.out.metadata, ch_genome, ch_regions, ch_blacklist )

    tobias_footprint( tobias_atacorrect.out.corrected, ch_regions )

    // tobias_footprint.out.footprints.flatten().branch {
    //     meta: it instanceof LinkedHashMap
    //     footprints: it instanceof String
    //     }
    //     .set{ch_footprint_split}

    // ch_footprint_split.meta
    //     .map { row -> row.sample_id }
    //     .set { ch_sample_ids }

    // tobias_bindetect( 
    //     ch_sample_ids.collect(),
    //     ch_footprint_split.footprints.collect(),
    //     awk.out.file_no_meta,
    //     ch_genome,
    //     ch_peaks
    //     )

    // tobias_footprint.out.footprints | view

    // ch_bind = tobias_footprint.out.footprints.collect{ it[0] }
    //     .merge(tobias_footprint.out.footprints.collect{ it[1] }) {a,b -> tuple(a,b)}
    //     .subscribe {log.info("$it")}
    

    // ch_bindtest = tobias_footprint.out.footprints.collect{ it[0] }
    //     .merge(tobias_footprint.out.footprints.collect{ it[1] }) {a,b -> tuple(a,b)}
    //     .combine(ch_motifs)
    //     .combine(ch_genome)
    //     .combine(ch_peaks)

    // tobias_bindetect( ch_bindtest )
    
    // ch_plotagreggate = tobias_bindetect.out.motifList
    //     .splitCsv(header:false)
    //     .merge(tobias_bindetect.out.motifBeds.toSortedList().flatten())
    //     .combine(tobias_atacorrect.out.corrected.collect{ it[1] }.map{row -> [row]}) 
   

    // tobias_plotaggregate(ch_plotagreggate)


   
}


//     fastq_metadata.out.metadata | view


// ch_bindtest.view { "$it" }
 // .subscribe {log.info("$it")}

    // tobias_bindetect.out.motifBeds
    //     .flatten()
    //     .subscribe {log.info("$it")}

  // .subscribe {log.info("$it")}

    // tobias_bindetect.out.motifBeds
    //     .flatten()
    //     .subscribe {log.info("$it")}


 // ch_plotagreggate = tobias_bindetect.out.motifBeds
    //     .combine(tobias_atacorrect.out.corrected.collect{ it[1] }) 
    // // ch_plotagreggate = tobias_atacorrect.out.corrected.collect{ it[1] }
    // //     .join(tobias_bindetect.out.motifBeds)

    
    // tobias_plotaggregate( ch_plotagreggate )
    // //tobias_bindetect.out.motifBeds | view
    // //tobias_atacorrect.out.corrected | view

// tuple val(sample_ids), path(signals), path(motifs), path(genome), path(peaks)

// workflow.onComplete {
//     log.info "\nPipeline finished executing\n"
// }

/*------------------------------------------------------------------------------------*/