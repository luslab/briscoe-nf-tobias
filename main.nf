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
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { tobias_atacorrect; tobias_footprint; tobias_bindetect; tobias_plotaggregate } from './modules/tobias/main.nf'
include { build_debug_param_summary; luslab_header; check_params } from './luslab-nf-modules/tools/luslab_util/main.nf'
include { bam_metadata } from './luslab-nf-modules/tools/metadata/main.nf'
include { sort_index_bam } from './luslab-nf-modules/workflows/bam_flows/main.nf'
include { samtools_faidx } from './luslab-nf-modules/tools/samtools/main.nf'
include { awk } from './luslab-nf-modules/tools/luslab_linux_tools/main.nf'

/*-----------------------------------------------------------------------------------------------------------------------------
Initialisation
-------------------------------------------------------------------------------------------------------------------------------*/

// Show banner
log.info luslab_header()

if(params.verbose){
    log.info build_debug_param_summary()
}

// Show work summary
def summary = [:]
summary['Show debug'] = params.verbose
summary['Output directory'] = params.outdir
summary['Trace directory'] = params.tracedir
summary['Max CPUs'] = params.max_cpus
summary['Max memory'] = params.max_memory
summary['Max time'] = params.max_time
summary['Design file'] = params.design
summary['Genome file'] = params.genome
summary['Genome index file'] = params.genome_index
summary['Regions file'] = params.regions
summary['Blacklist file'] = params.blacklist
summary['Peaks file'] = params.peaks
summary['Motifs file'] = params.motifs
summary['Skip bam indexing'] = params.skip_bam_index
summary['Skip genome indexing'] = params.skip_genome_index
summary['Motif bundle size'] = params.motif_bundle_count
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"

// Check parameters
check_params(['genome','regions','blacklist','design','motifs','peaks'])
if(params.skip_genome_index) check_params(['genome_index'])

/*-----------------------------------------------------------------------------------------------------------------------------
Channel Initialisation
-------------------------------------------------------------------------------------------------------------------------------*/

motifs_format = [[[:], params.motifs]]

ch_regions = Channel.value(file(params.regions, checkIfExists: true))
ch_blacklist = Channel.value(file(params.blacklist, checkIfExists: true))
ch_peaks = Channel.value(file(params.peaks, checkIfExists: false))

Channel.from(motifs_format)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_motifs}

if(!params.skip_genome_index) {
    ch_genome = Channel.value(file(params.genome, checkIfExists: true))
}
else {
    genome_files = [params.genome, params.genome_index]
    Channel.from(genome_files)
        .map { row -> [ file(row[0], checkIfExists: true), file(row[1], checkIfExists: true)]}
        .set {ch_genome}
}

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

workflow {
    bam_metadata(params.design)
    //bam_metadata.out.metadata | view

    // Check if we need to index the bam
    ch_bam_bai = bam_metadata.out.metadata
    if(!params.skip_bam_index) {
        sort_index_bam(bam_metadata.out.metadata)
        ch_bam_bai = sort_index_bam.out.bam_bai
        //ch_bam_bai | view
    }

    // Check if we need to index the genome
    if(!params.skip_genome_index) {
        samtools_faidx( params.modules['cust_samtools_faidx'], ch_genome )
        ch_genome = samtools_faidx.out.indexedFasta
        //ch_genome | view
    }

    awk(params.modules['motifsplit_awk'], ch_motifs)
    //awk.out.file_no_meta | view

    tobias_atacorrect( params.modules['tobias_atacorrect'], ch_bam_bai, ch_genome, ch_regions, ch_blacklist )
    //tobias_atacorrect.out.corrected | view

    tobias_footprint( params.modules['tobias_footprint'], tobias_atacorrect.out.corrected, ch_regions )
    //tobias_footprint.out.footprints | view

    tobias_footprint.out.footprints.flatten().branch {
            meta: it.getClass().toString() == "class java.util.LinkedHashMap"
            footprints: it.getClass().toString() == "class sun.nio.fs.UnixPath"
         }
         .set{ch_footprint_split}
    //ch_footprint_split.meta | view
    //ch_footprint_split.footprints | view

    ch_footprint_split.meta
        .map { row -> row.sample_id }
        .set { ch_sample_ids }
    //ch_sample_ids | view

    // Group motifs into bundles for processing
    ch_motif_bundles = awk.out.file_no_meta.flatten().collate(params.motif_bundle_count, false)
    //ch_motif_bundles | view

    tobias_bindetect(
        params.modules['tobias_bindetect'],
        ch_sample_ids.collect(),
        ch_footprint_split.footprints.collect(),
        ch_motif_bundles,
        ch_genome,
        ch_peaks
    )
    // tobias_bindetect.out.report | view
    // tobias_bindetect.out.bindetect_files | view
    // tobias_bindetect.out.overview | view
    // tobias_bindetect.out.images | view
    // tobias_bindetect.out.plots | view
    // tobias_bindetect.out.motif_beds | view
    // tobias_bindetect.out.bound_beds | view
    // tobias_bindetect.out.unbound_beds | view
    // tobias_bindetect.out.motif_list | view

    ch_sorted_motif_names = tobias_bindetect.out.motif_list
        .splitCsv(header:false)
        .flatten()
    //ch_sorted_motif_names.collect() | view

    ch_sorted_beds = tobias_bindetect.out.motif_beds
        .flatten()
    //ch_sorted_beds.collect() | view

    tobias_plotaggregate(
        params.modules['tobias_plotaggregate'],
        ch_sorted_motif_names,
        ch_sorted_beds,
        tobias_atacorrect.out.corrected_no_meta.collect())
    //tobias_plotaggregate.out.pdf | view
}

// Syntax help

// fastq_metadata.out.metadata | view

// ch_bindtest.view { "$it" }

// .subscribe {log.info("$it")}

/*------------------------------------------------------------------------------------*/