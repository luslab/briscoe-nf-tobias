#!/usr/bin/env nextflow
/*
========================================================================================
                         TOBIAS wrapper module
========================================================================================
 #### Homepage / Documentation

 https://github.com/loosolab/TOBIAS

----------------------------------------------------------------------------------------
Module decription
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
*/

// Specify DSL2
nextflow.enable.dsl = 2

process tobias_atacorrect {
    tag "${meta.sample_id}"
    label 'h_cpu'
    time = 4.h

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-tobias:latest'
    
    input:
        val opts
        tuple val(meta), path(bam), path(bai)
        tuple path(genome), path(genome_index)
        path bed
        path blacklist

    output:
        path "*.log", emit: report
        tuple val(meta), path("*_corrected.bw"), emit: corrected
        path "*_corrected.bw", emit: corrected_no_meta
        tuple val(meta), path("*_expected.bw"), emit: expected
        tuple val(meta), path("*_uncorrected.bw"), emit: uncorrected

    script:
        command = "TOBIAS ATACorrect --bam $bam --genome $genome --peaks $bed --blacklist $blacklist --outdir . --cores ${task.cpus} > ${meta.sample_id}_atacorrect.log"
        if (params.verbose){
          println ("[MODULE] tobias/atacorrect command: " + command)
        }
    """
    ${command}
    """
}

process tobias_footprint {
    tag "${meta.sample_id}"
    label 'h_cpu'
    time = 4.h

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-tobias:latest'

    input:
        val opts
        tuple val(meta), path(corrected)
        path bed

    output:
        val opts
        path "*.log", emit: report 
        tuple val(meta), path("*_footprints.bw"), emit: footprints
        path "*_footprints.bw", emit: footprintsDetect

    script:
        command = "TOBIAS FootprintScores --signal $corrected --regions $bed --output ${meta.sample_id}_footprints.bw --cores ${task.cpus} > ${meta.sample_id}_footprint.log"
        if (params.verbose){
          println ("[MODULE] tobias/footprint command: " + command)
        }
    """
    ${command}
    """
}

process tobias_bindetect {
    tag "${run_hash}"
    label 'mx_cpu'
    time = 24.h

    publishDir "${params.outdir}/${opts.publish_dir}/${run_hash}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-tobias:latest'

    input:
        val opts
        val sample_ids
        path signals
        path motifs
        tuple path(genome), path(genome_index)
        path peaks

    output:
        path "*.log", emit: report
        path "bindetect_*", emit: bindetect_files
        path "*/*_overview.txt", emit: overview
        path "*/*.png", emit: images
        path "*/plots/**", emit: plots
        path "*/beds/*_all.bed", emit: motif_beds
        path "*/beds/*_bound.bed", emit: bound_beds
        path "*/beds/*_unbound.bed", emit: unbound_beds
        path "motiflist.txt", emit: motif_list

    script:
        run_hash = "${motifs}".md5()
        command = "TOBIAS BINDetect --motifs $motifs --signal $signals --peaks $peaks --genome $genome --outdir . --cond_names ${sample_ids.join(' ')} --cores ${task.cpus} > bindectect.log"
        if (params.verbose){
          println ("[MODULE] tobias/bindetect command: " + command)
        }

        dev_command = "fil-profile run /home/TOBIAS/tobias/tools/bindetect.py --motifs $motifs --signal $signals --peaks $peaks --genome $genome --outdir . --cores ${task.cpus} > bindectect.log"
    """
    awk '{if(\$1 ~ />/) {a=\$2"_"\$1; gsub(/>/,"",a); gsub("::","",a); gsub(/[()]/,"",a); print a}}' $motifs | sort > motiflist.txt
    source activate nfcore-module-tobias
    ${dev_command}
    """
}
process tobias_plotaggregate {
    tag "${motif_name}"
    label 'mn_cpu'
    time = 4.h

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-tobias:latest'

    input:
        val opts
        val motif_name
        path motif_beds
        path corrected_insertions

    output:
        path "$motif_name/*.pdf", emit: pdf

    script:
    """
    mkdir $motif_name
    TOBIAS PlotAggregate --TFBS $motif_beds --signals ${corrected_insertions.join(' ')} --output $motif_name/${motif_name}_plotaggregate.pdf --share_y both --plot_boundaries --signal-on-x
    """
}