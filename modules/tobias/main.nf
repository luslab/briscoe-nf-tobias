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

    container 'quay.io/biocontainers/tobias:0.12.10--py37h77a2a36_1'
    
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

/*
----------------------------------------------------------------------------------------
*/

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

    container 'quay.io/biocontainers/tobias:0.12.10--py37h77a2a36_1'

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

/*
----------------------------------------------------------------------------------------
*/

// Usage:
// TOBIAS BINDetect --signals <bigwig1> (<bigwig2> (...)) --motifs <motifs.txt> --genome
// <genome.fasta> --peaks <peaks.bed>
// Output files:
// - <outdir>/<prefix>_figures.pdf
// - <outdir>/<prefix>_results.{txt,xlsx}
// - <outdir>/<prefix>_distances.txt
// - <outdir>/<TF>/<TF>_overview.{txt,xlsx} (per motif)
// - <outdir>/<TF>/beds/<TF>_all.bed (per motif)
// - <outdir>/<TF>/beds/<TF>_<condition>_bound.bed (per motif-condition pair)
// - <outdir>/<TF>/beds/<TF>_<condition>_unbound.bed (per motif-condition pair)
// ------------------------------------------------------------------------------------------
// Required arguments:
//   --signals [<bigwig> [<bigwig> ...]]
//                                    Signal per condition (.bigwig format)
//   --peaks <bed>                    Peaks.bed containing open chromatin regions across all
//                                    conditions
//   --motifs [<motifs> [<motifs> ...]]
//                                    Motif file(s) in pfm/jaspar/meme format
//   --genome <fasta>                 Genome .fasta file
// Optional arguments:
//   --cond-names [<name> [<name> ...]]
//                                    Names of conditions fitting to --signals (default:
//                                    prefix of --signals)
//   --peak-header <file>             File containing the header of --peaks separated by
//                                    whitespace or newlines (default: peak columns are named
//                                    "_additional_<count>")
//   --naming <string>                Naming convention for TF output files ('id', 'name',
//                                    'name_id', 'id_name') (default: 'name_id')
//   --motif-pvalue <float>           Set p-value threshold for motif scanning (default:
//                                    1e-4)
//   --bound-pvalue <float>           Set p-value threshold for bound/unbound split (default:
//                                    0.001)
//   --pseudo <float>                 Pseudocount for calculating log2fcs (default: estimated
//                                    from data)
//   --time-series                    Will only compare signals1<->signals2<->signals3 (...)
//                                    in order of input, and skip all-against-all comparison.
//   --skip-excel                     Skip creation of excel files - for large datasets, this
//                                    will speed up BINDetect considerably
//   --output-peaks <bed>             Gives the possibility to set the output peak set
//                                    differently than the input --peaks. This will limit all
//                                    analysis to the regions in --output-peaks. NOTE:
//                                    --peaks must still be set to the full peak set!
//   --prefix <prefix>                Prefix for overview files in --outdir folder (default:
//                                    bindetect)
// Run arguments:
//   --outdir <directory>             Output directory to place TFBS/plots in (default:
//                                    bindetect_output)
//   --cores <int>                    Number of cores to use for computation (default: 1)
//   --split <int>                    Split of multiprocessing jobs (default: 100)
//   --verbosity <int>                Level of output logging (0: silent, 1: errors/warnings,
//                                    2: info, 3: stats, 4: debug, 5: spam) (default: 3)

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

    container 'quay.io/biocontainers/tobias:0.12.10--py37h77a2a36_1'

    input:
        val opts
        val list_command
        path signals
        path motifs
        tuple path(genome), path(genome_index)
        path peaks
        path output_peaks

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
    def command = ""
    if(params.output_peaks) {
        command = "TOBIAS BINDetect --motifs $motifs --signal $signals --peaks $peaks --genome $genome --output_peaks $output_peaks --outdir . --cores ${task.cpus} > bindetect.log"
    } else {
        command = "TOBIAS BINDetect --motifs $motifs --signal $signals --peaks $peaks --genome $genome --outdir . --cores ${task.cpus} > bindetect.log" // --cond_names ${sample_ids.join(' ')}
    }
    
    if (params.verbose){
        println ("[MODULE] tobias/bindetect command: " + command)
    }
    """
    $list_command $motifs | sort > motiflist.txt
    ${command}
    """
}

/*
----------------------------------------------------------------------------------------
*/

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

    container 'quay.io/biocontainers/tobias:0.12.10--py37h77a2a36_1'

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