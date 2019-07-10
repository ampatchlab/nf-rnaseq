#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-rnaseq: RNA-Seq nextflow pipeline
 *
 * Copyright (C) 2019 QIMR Berghofer Medical Research Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


import nextflow.splitter.CsvSplitter
import nextflow.config.ConfigParser

nextflow_config = file("${baseDir}/nextflow.config").text
parsed_config = new ConfigParser().setIgnoreIncludes(true).parse(nextflow_config)
defaults = parsed_config.params

check_params()


// log required arguments
log.info("Input csv: ${params.csv}")
log.info("Reference genome: ${params.genome ?: null}")

// log optional arguments
log.info("Paired-end readgroups: ${String.valueOf(params.paired_end)}")
log.info("ReadGroup ID separator: ${params.rgid_sep}")
log.info("Adapters to trim: ${params.adapters}")

// log reference genome params
log.info("Reference FASTA file: ${params.fasta}")
log.info("Reference GTF file: ${params.gtf}")

// log STAR genome generate params
log.info("STAR size of the bins for genome storage: ${params.star_genome_chr_bin_n_bits}")
log.info("STAR length of the SA pre-indexing string: ${params.star_genome_sa_index_n_bases}")
log.info("STAR SJDB overhang: ${params.star_sjdb_overhang}")

// log Cutadapt params
log.info("Cutadapt R1 adapter: ${params.cutadapt_r1_adapter ?: null}")
log.info("Cutadapt R2 adapter: ${params.cutadapt_r2_adapter ?: null}")
log.info("Cutadapt base quality cutoff: ${params.cutadapt_base_qual_cutoff}")
log.info("Cutadapt min read length: ${params.cutadapt_min_read_length}")

// log RSEM mparams
log.info("RSEM fragment length min: ${params.rsem_fragment_length_min}")
log.info("RSEM fragment length max: ${params.rsem_fragment_length_max}")
log.info("RSEM fragment length mean: ${params.rsem_fragment_length_mean}")
log.info("RSEM fragment length sd: ${params.rsem_fragment_length_sd}")
log.info("RSEM estimate RSPD: ${params.rsem_estimate_rspd}")
log.info("RSEM number of RSPD bins: ${params.rsem_num_rspd_bins}")
log.info("RSEM seed length: ${params.rsem_seed_length}")

// log RSEM Gibbs params
log.info("RSEM Gibbs burn-in: ${params.gibbs_burnin}")
log.info("RSEM Gibbs number of samples: ${params.gibbs_num_samples}")
log.info("RSEM Gibbs sampling gap: ${params.gibbs_sampling_gap}")

// log RSEM CalcCI params
log.info("RSEM CalcCI credibility level: ${params.ci_credibility_level}")
log.info("RSEM CalcCI number of samples per CV: ${params.ci_num_samples_per_count_vector}")

// log MultiQC params
log.info("MultiQC config: ${params.multiqc_config}")

// log output params
log.info("Reference directory: ${params.refdir}")
log.info("Output directory: ${params.outdir}")

// log reports params
log.info("Execution report: ${config.report.file}")
log.info("Trace report: ${config.trace.file}")
log.info("Timeline report: ${config.timeline.file}")
log.info("Flowchart: ${config.dag.file}")

// log AWS Batch params
log.info("AWS Batch JobQueue: ${params.aws_queue}")
log.info("AWS Region: ${params.aws_region}")

validate_input_csv()


// input files
csv_file = file(params.csv)
fasta_file = file(params.fasta)
gtf_file = file(params.gtf)

// asset files
multiqc_cfg = file(params.multiqc_config)

// reference directories
STAR = params.genome ?: 'STAR'
RSEM = params.genome ?: 'RSEM'


/*
 * PREPROCESSING - Gunzip reference FASTA file
 */
process gunzip_fasta {
    storeDir params.refdir

    input:
    file fasta from fasta_file

    output:
    file "${fasta.getBaseName()}" into gunzipped_fasta

    when:
    fasta.getExtension() == "gz"

    """
    gzip -dc "${fasta}" > "${fasta.getBaseName()}"
    """
}


/*
 * PREPROCESSING - Gunzip reference GTF file
 */
process gunzip_gtf {
    storeDir params.refdir

    input:
    file gtf from gtf_file

    output:
    file "${gtf.getBaseName()}" into gunzipped_gtf

    when:
    gtf.getExtension() == "gz"

    """
    gzip -dc "${gtf}" > "${gtf.getBaseName()}"
    """
}


/*
 * PREPROCESSING - Create a value channel for the reference FASTA file
 */
gunzipped_fasta
    .ifEmpty { fasta_file }
    .set { ref_fasta }


/*
 * PREPROCESSING - Create a value channel for the reference GTF file
 */
gunzipped_gtf
    .ifEmpty { gtf_file }
    .set { ref_gtf }


/*
 * PREPROCESSING - STAR index
 */
process star_index {
    storeDir "${params.refdir}/STAR"

    label 'star'

    input:
    file ref_fasta
    file ref_gtf

    output:
    file "${STAR}" into star_index

    """
    mkdir "${STAR}"
    STAR \\
        --runThreadN "${task.cpus}" \\
        --runMode genomeGenerate \\
        --genomeDir "${STAR}" \\
        --genomeFastaFiles "${ref_fasta}" \\
        --genomeChrBinNbits "${params.star_genome_chr_bin_n_bits}" \\
        --genomeSAindexNbases "${params.star_genome_sa_index_n_bases}" \\
        --sjdbGTFfile "${ref_gtf}" \\
        --sjdbOverhang "${params.star_sjdb_overhang}"
    """
}


/*
 * PREPROCESSING - SAMtools faidx
 */
process samtools_faidx {
    storeDir params.refdir

    label 'samtools'

    input:
    file ref_fasta

    output:
    file "${ref_fasta}.fai" into ref_faidx

    """
    samtools faidx "${ref_fasta}"
    """
}


/*
 * PREPROCESSING - For RNA-SeQC, gene features must have a 'transcript_id' attribute
 */
process create_rnaseqc_gtf {
    storeDir "${params.refdir}/RNA-SeQC"

    label 'gffutils'

    input:
    file ref_gtf

    output:
    file "${ref_gtf.getBaseName()}.rnaseqc.gtf" into rnaseqc_ref_gtf

    """
    #!/usr/bin/env python3
    import gffutils

    with open("${ref_gtf.getBaseName()}.rnaseqc.gtf", 'w') as output_gtf:
        di = gffutils.DataIterator("${ref_gtf}")
        for directive in di.directives:
            print('##' + directive, file=output_gtf)
        for feature in di:
            if feature.featuretype == 'gene' and 'transcript_id' not in feature.attributes:
                print(feature, 'transcript_id "";', file=output_gtf)
            else:
                print(feature, file=output_gtf)
    """
}


/*
 * PREPROCESSING - For RSeQC, create a BED12 input file
 */
process gtftogenepred {
    storeDir "${params.refdir}/UCSC"

    label 'ucsc_gtftogenepred'

    input:
    file ref_gtf

    output:
    file "${ref_gtf.getBaseName()}.pred" into rseqc_ref_pred

    """
    gtfToGenePred "${ref_gtf}" "${ref_gtf.getBaseName()}.pred"
    """
}


process genepredtobed {
    storeDir "${params.refdir}/UCSC"

    label 'ucsc_genepredtobed'

    input:
    file rseqc_ref_pred

    output:
    file "${rseqc_ref_pred.getBaseName()}.bed12" into rseqc_ref_bed12

    """
    genePredToBed "${rseqc_ref_pred}" "${rseqc_ref_pred.getBaseName()}.bed12"
    """
}


/*
 * PREPROCESSING - RSEM extract reference transcripts
 */
process rsem_extract_reference_transcripts {
    storeDir "${params.refdir}/RSEM"

    label 'rsem'

    input:
    file ref_gtf
    file ref_fasta

    output:
    file "${RSEM}.{grp,ti,chrlist}" into rsem_index_files
    file "${RSEM}.transcripts.fa" into rsem_transcripts_fasta

    """
    rsem-extract-reference-transcripts \\
        "${RSEM}" \\
        0 \\
        "${ref_gtf}" \\
        None \\
        0 \\
        "${ref_fasta}"
    """
}


/*
 * PREPROCESSING - RSEM prepare reference
 */
process rsem_preref {
    storeDir "${params.refdir}/RSEM"

    label 'rsem'

    input:
    file rsem_transcripts_fasta

    output:
    file "${RSEM}.seq" into rsem_seq_index
    file "${RSEM}.idx.fa"
    file "${RSEM}.n2g.idx.fa"

    """
    rsem-preref \\
        "${rsem_transcripts_fasta}" \\
        1 \\
        "${RSEM}"
    """
}


/*
 * PREPROCESSING - Extract gene and transcript annotations from the input GTF
 */
process extract_annotations {
    storeDir "${params.refdir}/annotation"

    label 'gffutils'

    input:
    file ref_gtf

    output:
    file 'genes.tsv' into gene_annotation
    file 'transcripts.tsv' into transcript_annotation

    """
    #!/usr/bin/env python3
    import csv
    import gffutils

    def iter_features(gtf, featuretype, fieldnames):
        for feature in gffutils.DataIterator(gtf):
            if feature.featuretype == featuretype:
                yield { fieldname: ','.join(feature[fieldname]) for fieldname in fieldnames }

    gene_fieldnames = ['gene_id', 'gene_name', 'gene_type']
    transcript_fieldnames = ['transcript_id', 'transcript_name', 'transcript_type']

    kwargs = {'delimiter': '\\t', 'lineterminator': '\\n'}

    with open('genes.tsv', 'w') as genes_tsv:
        writer = csv.DictWriter(genes_tsv, fieldnames=gene_fieldnames, **kwargs)
        writer.writeheader()
        writer.writerows(iter_features("${ref_gtf}", 'gene', gene_fieldnames))

    with open('transcripts.tsv', 'w') as transcripts_tsv:
        writer = csv.DictWriter(transcripts_tsv, fieldnames=transcript_fieldnames, **kwargs)
        writer.writeheader()
        writer.writerows(iter_features("${ref_gtf}", 'transcript', transcript_fieldnames))
    """
}


/*
 * PREPROCESSING - Tag the genome and transcript FASTA files
 */
ref_fasta
    .map { tuple('genome', it) }
    .set { genome_sizes_inputs }

rsem_transcripts_fasta
    .map { tuple('transcript', it) }
    .set { transcript_sizes_inputs }


/*
 * PREPROCESSING - Create genome and transcript chrom.sizes files
 */
process fatotwobit {
    tag { tag }

    storeDir "${params.refdir}/UCSC"

    label 'ucsc_fatotwobit'

    input:
    set tag, file(fasta) from genome_sizes_inputs.mix(transcript_sizes_inputs)

    output:
    set tag, file("${fasta.getBaseName()}.2bit") into tagged_twobits

    """
    faToTwoBit "${fasta}" "${fasta.getBaseName()}.2bit"
    """
}


process twobitinfo {
    tag { tag }

    storeDir "${params.refdir}/UCSC"

    label 'ucsc_twobitinfo'

    input:
    set tag, file(twobit) from tagged_twobits

    output:
    set tag, file("${twobit.getBaseName()}.chrom.sizes") into tagged_chrom_sizes

    """
    twoBitInfo "${twobit}" stdout | sort -k2rn > "${twobit.getBaseName()}.chrom.sizes"
    """
}


/*
 * PREPROCESSING - Parse in the readgroups and create the inputs for FastQC and Cutadapt
 */
(rgids, r1_inputs, r2_inputs) = Channel
    .fromPath( params.csv )
    .splitCsv( header:true )
    .ifEmpty { exit 1, "No readgroups found post-validation. Exiting." }
    .separate(3) { row ->
        def rgid = row.readgroup ? [row.sample, row.readgroup].join(params.rgid_sep) : row.sample

        // files to stage
        def fq1 = params.paired_end ? file(row.fastq1) : file(row.fastq)
        def fq2 = params.paired_end ? file(row.fastq2) : 'NO_FILE'

        // filenames to stage with
        def fq1_fn = params.paired_end ? "${rgid}.1.${get_fastq_extn(fq1)}" : "${rgid}.${get_fastq_extn(fq1)}"
        def fq2_fn = params.paired_end ? "${rgid}.2.${get_fastq_extn(fq2)}" : 'NO_FILE'

        tuple(tuple(row.sample, rgid), tuple(fq1_fn, fq1), tuple(fq2_fn, fq2))
    }

rgids.into { samples; fastqc_raw_rgids; cutadapt_rgids }
r1_inputs.into { fastqc_raw_r1_inputs; cutadapt_r1_inputs }
r2_inputs.into { fastqc_raw_r2_inputs; cutadapt_r2_inputs }


/*
 * STEP 1 - Run FastQC on the raw reads
 */
process fastqc_raw {
    tag { rgid }

    label 'fastqc'

    publishDir "${params.outdir}/FastQC/${rgid}/raw"

    input:
    set sample, rgid from fastqc_raw_rgids
    set fq1, file("${fq1}") from fastqc_raw_r1_inputs
    set fq2, file("${fq2}") from fastqc_raw_r2_inputs

    output:
    file "*_fastqc.{zip,html}" into fastqc_raw_results

    script:
    if (params.paired_end) {

        """
        fastqc -q "${fq1}" "${fq2}"
        """

    } else {

        """
        fastqc -q "${fq1}"
        """
    }
}


/*
 * STEP 2 - Trim adapters using Cutadapt
 */
process cutadapt {
    tag { rgid }

    label 'cutadapt'

    publishDir "${params.outdir}/Cutadapt/${rgid}"

    input:
    set sample, rgid from cutadapt_rgids
    set fq1, file("${rgid}/${fq1}") from cutadapt_r1_inputs
    set fq2, file("${rgid}/${fq2}") from cutadapt_r2_inputs

    output:
    set rgid, file("*.fastq.gz") into fastqc_trimmed_inputs, trimmed_readgroups
    file "*.log" into cutadapt_logs

    script:
    if (params.paired_end) {

        """
        cutadapt \\
            -a "${params.cutadapt_r1_adapter}" \\
            -A "${params.cutadapt_r2_adapter}" \\
            -q "${params.cutadapt_base_qual_cutoff}" \\
            -m "${params.cutadapt_min_read_length}" \\
            --trim-n \\
            -o "${rgid}.1.fastq.gz" \\
            -p "${rgid}.2.fastq.gz" \\
            "${rgid}/${fq1}" \\
            "${rgid}/${fq2}" \\
            > "${rgid}.log"
        """

    } else {

        """
        cutadapt \\
            -a "${params.cutadapt_r1_adapter}" \\
            -q "${params.cutadapt_base_qual_cutoff}" \\
            -m "${params.cutadapt_min_read_length}" \\
            --trim-n \\
            -o "${rgid}.fastq.gz" \\
            "${rgid}/${fq1}" \\
            > "${rgid}.log"
        """
    }
}


/*
 * STEP 3 - Run FastQC on the trimmed reads
 */
process fastqc_trimmed {
    tag { rgid }

    label 'fastqc'

    publishDir "${params.outdir}/FastQC/${rgid}/trimmed"

    input:
    set rgid, file(fastqs) from fastqc_trimmed_inputs

    output:
    file "*_fastqc.{zip,html}" into fastqc_trimmed_results

    script:
    if (params.paired_end) {

        """
        fastqc -q "${rgid}.1.fastq.gz" "${rgid}.2.fastq.gz"
        """

    } else {

        """
        fastqc -q "${rgid}.fastq.gz"
        """
    }
}


/*
 * STEP 4 - Group the trimmed readgroups by sample
 */
samples
    .groupTuple()
    .map { sample, rgids ->
        tuple( groupKey(sample, rgids.size()), rgids )
    }
    .transpose()
    .map { sample, readgroup ->
        tuple(readgroup, sample)
    }
    .join( trimmed_readgroups )
    .groupTuple( by:1 )
    .map { rgids, sample, fastqs ->
        tuple(sample.toString(), rgids, fastqs.flatten())
    }
    .set { star_inputs }


/*
 * STEP 5 - Align the trimmed reads using STAR
 */
process star {
    tag { sample }

    label 'star'

    publishDir "${params.outdir}/STAR/${sample}"

    input:
    set sample, rgids, file(fastqs) from star_inputs
    file index from star_index.collect()

    output:
    set sample, file("*.Aligned.out.bam") into star_alignments
    set sample, file("*.Aligned.sortedByCoord.out.bam") into star_csorted_bam_files, star_sorted_alignments
    set sample, file("*.Aligned.toTranscriptome.out.bam") into star_transcriptome_alignments
    set sample, file("*.ReadsPerGene.out.tab") into star_gene_counts
    set sample, file("*.out.bg") into star_signal_output
    set sample, file("*.SJ.out.tab") into star_splice_junctions
    file "*.out" into star_logs

    script:
    def rgs = rgids.collect { /"ID:${it}" "SM:${sample}"/ }.join(" , ")

    if (params.paired_end) {

        def fq1 = rgids.collect { "${it}.1.fastq.gz" }.join(',')
        def fq2 = rgids.collect { "${it}.2.fastq.gz" }.join(',')

        """
        STAR \\
            --genomeDir "${index}" \\
            --quantMode TranscriptomeSAM GeneCounts \\
            --readFilesCommand zcat \\
            --readFilesIn "${fq1}" "${fq2}" \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outSAMattributes All \\
            --outSAMattrRGline ${rgs} \\
            --outSAMtype BAM Unsorted SortedByCoordinate \\
            --outSAMunmapped Within KeepPairs \\
            --outWigType bedGraph \\
            --outFileNamePrefix "${sample}."
        """

    } else {

        def fq = rgids.collect { "${it}.fastq.gz" }.join(',')

        """
        STAR \\
            --genomeDir "${index}" \\
            --quantMode TranscriptomeSAM GeneCounts \\
            --readFilesCommand zcat \\
            --readFilesIn "${fq}" \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outSAMattributes All \\
            --outSAMattrRGline ${rgs} \\
            --outSAMtype BAM Unsorted SortedByCoordinate \\
            --outSAMunmapped Within KeepPairs \\
            --outWigType bedGraph \\
            --outFileNamePrefix "${sample}."
        """
    }
}


/*
 * STEP 5.1 - Index the coordinate-sorted alignments
 */
process samtools_index {
    tag { sample }

    label 'samtools'

    publishDir "${params.outdir}/STAR/${sample}"

    input:
    set sample, file(bam) from star_csorted_bam_files

    output:
    file "*.bai"

    """
    samtools index "${bam}"
    """
}


/*
 * STEP 6 - Mark duplicates using Picard
 */
process mark_duplicates {
    tag { sample }

    label 'picard'

    publishDir "${params.outdir}/MarkDuplicates/${sample}", saveAs: { fn ->
        fn.endsWith(".bai") ? "${sample}.bam.bai" : "${fn}"
    }

    input:
    set sample, file(bam) from star_sorted_alignments

    output:
    set sample, file("*.bam"), file("*.bai") into rseqc_inputs, rnaseqc_inputs
    file "*.metrics.txt" into mark_duplicates_metrics

    """
    picard \\
        -Xmx${task.memory.toGiga()}g \\
        -XX:+UseSerialGC \\
    MarkDuplicates \\
        INPUT="${bam}" \\
        OUTPUT="${sample}.bam" \\
        METRICS_FILE="${sample}.metrics.txt" \\
        ASSUME_SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT
    """
}


/*
 * STEP 7 - RSeQC
 */
rseqc_inputs
    .into {
        rseqc_bam_stat_inputs
        rseqc_gene_body_coverage_inputs
        rseqc_infer_experiment_inputs
        rseqc_inner_distance_inputs
        rseqc_junction_annotation_inputs
        rseqc_junction_saturation_inputs
        rseqc_read_distribution_inputs
        rseqc_read_duplication_inputs
        rseqc_read_gc_inputs
    }


/*
 * STEP 7.1 - bam_stat.py
 */
process bam_stat {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/bam_stat"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_bam_stat_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_bam_stat_results

    """
    bam_stat.py -i "${bam}" > "${sample}.txt"
    """
}


/*
 * STEP 7.2 - geneBody_coverage.py
 */
process gene_body_coverage {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/gene_body_coverage"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_gene_body_coverage_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_gene_body_coverage_results

    """
    geneBody_coverage.py -i "${bam}" -r "${bed}" -o "${sample}"
    """
}


/*
 * STEP 7.3 - infer_experiment.py
 */
process infer_experiment {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/infer_experiment"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_infer_experiment_inputs
    file bed from rseqc_ref_bed12

    output:
    set sample, file("${sample}.txt") into strandedness_inputs
    file "${sample}.*" into rseqc_infer_experiment_results

    """
    infer_experiment.py -i "${bam}" -r "${bed}" > "${sample}.txt"
    """
}


/*
 * STEP 7.4 - inner_distance.py
 */
process inner_distance {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/inner_distance"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_inner_distance_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_inner_distance_results

    """
    inner_distance.py -i "${bam}" -r "${bed}" -o "${sample}"
    """
}


/*
 * STEP 7.5 - junction_annotation.py
 */
process junction_annotation {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/junction_annotation"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_junction_annotation_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_junction_annotation_results

    """
    junction_annotation.py -i "${bam}" -r "${bed}" -o "${sample}" 2> "${sample}.txt"
    """
}


/*
 * STEP 7.6 - junction_saturation.py
 */
process junction_saturation {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/junction_saturation"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_junction_saturation_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_junction_saturation_results

    """
    junction_saturation.py -i "${bam}" -r "${bed}" -o "${sample}"
    """
}


/*
 * STEP 7.7 - read_distribution.py
 */
process read_distribution {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/read_distribution"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_read_distribution_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_read_distribution_results

    """
    read_distribution.py -i "${bam}" -r "${bed}" > "${sample}.txt"
    """
}


/*
 * STEP 7.8 - read_duplication.py
 */
process read_duplication {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/read_duplication"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_read_duplication_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_read_duplication_results

    """
    read_duplication.py -i "${bam}" -o "${sample}"
    """
}


/*
 * STEP 7.9 - read_GC.py
 */
process read_gc {
    tag { sample }

    label 'rseqc'

    publishDir "${params.outdir}/RSeQC/read_gc"

    input:
    set sample, file(bam), file("${bam}.bai") from rseqc_read_gc_inputs
    file bed from rseqc_ref_bed12

    output:
    file "${sample}.*" into rseqc_read_gc_results

    """
    read_GC.py -i "${bam}" -o "${sample}"
    """
}


/*
 * STEP 8 - RNA-SeQC
 */
process rna_seqc {
    tag { sample }

    label 'rna_seqc'

    publishDir "${params.outdir}/RNA-SeQC/samples"

    input:
    set sample, file(bam), file("${bam}.bai") from rnaseqc_inputs
    file ref_gtf from rnaseqc_ref_gtf
    file ref_fasta
    file ref_faidx

    output:
    file sample into rnaseqc_results, rnaseqc_summary_inputs

    script:
    if (params.paired_end) {

        """
        rna-seqc \\
            -Xmx${task.memory.toGiga()}g \\
            -XX:+UseSerialGC \\
            -s "${sample}\\|${bam}\\|-" \\
            -r "${ref_fasta}" \\
            -t "${ref_gtf}" \\
            -o "${sample}"
        """

    } else {

        """
        rna-seqc \\
            -Xmx${task.memory.toGiga()}g \\
            -XX:+UseSerialGC \\
            -singleEnd \\
            -s "${sample}\\|${bam}\\|-" \\
            -r "${ref_fasta}" \\
            -t "${ref_gtf}" \\
            -o "${sample}"
        """
    }
}


/*
 * STEP 8.1 - Summarize the RNA-SeQC metrics
 */
process rna_seqc_summary {

    label 'python'

    publishDir "${params.outdir}/RNA-SeQC"

    input:
    file 'samples/*' from rnaseqc_summary_inputs.collect()

    output:
    file 'metrics.tsv'

    """
    #!/usr/bin/env python3
    import csv
    import os

    fieldnames = []
    data = dict()

    for sample in os.listdir('samples'):
        tsv = "samples/{}/metrics.tsv".format(sample)
        with open(tsv, 'r') as input_tsv:
            reader = csv.DictReader(input_tsv, delimiter='\\t')
            if not fieldnames:
                fieldnames = sorted(reader.fieldnames)
                fieldnames.insert(0, fieldnames.pop(fieldnames.index('Sample')))
                fieldnames.insert(1, fieldnames.pop(fieldnames.index('Note')))
            row = next(reader)
            data [ row['Sample'] ] = row

    with open('metrics.tsv', 'w') as output_tsv:
        writer = csv.DictWriter(output_tsv, fieldnames=fieldnames, delimiter='\\t', lineterminator='\\n')
        writer.writeheader()
        writer.writerows(data[sample] for sample in sorted(data))
    """
}


/*
 * STEP 9 - Infer a value for library strandedness
 */
process infer_strandedness {
    tag { sample }

    label 'python'

    publishDir "${params.outdir}/Strandedness/samples"

    input:
    set sample, file(inferred_experiment) from strandedness_inputs

    output:
    set sample, file("${sample}.strandedness") into strandedness_results
    file("${sample}.strandedness") into strandedness_summary

    """
    #!/usr/bin/env python3
    import re

    def infer_strandedness(sense, antisense):
        if 0.4 <= sense <= 0.6 and 0.4 <= antisense <= 0.6:
            return 'none'
        if sense >= 0.7 and antisense <= 0.3:
            return 'forward'
        if sense <= 0.3 and antisense >= 0.7:
            return 'reverse'
        return None

    regexes = {
        'pe_sense': r"\\"1\\+\\+,1--,2\\+-,2-\\+\\": (\\d\\.\\d+)",
        'pe_antisense': r"\\"1\\+-,1-\\+,2\\+\\+,2--\\": (\\d\\.\\d+)",
        'se_sense': r"\\"\\+\\+,--\\": (\\d\\.\\d+)",
        'se_antisense': r"\\"\\+-,-\\+\\": (\\d\\.\\d+)",
    }

    with open('${inferred_experiment}', 'r') as infile:
        data = infile.read()

    d = dict()
    for k, r in regexes.items():
        m = re.search(r, data, re.MULTILINE)
        if m:
            d[k] = float(m.group(1))

    pe_strandedness = None
    if ('pe_sense' in d and 'pe_antisense' in d):
        pe_strandedness = infer_strandedness(d['pe_sense'], d['pe_antisense'])
    se_strandedness = None
    if ('se_sense' in d and 'se_antisense' in d):
        se_strandedness = infer_strandedness(d['se_sense'], d['se_antisense'])

    strandedness = None
    if (pe_strandedness and not se_strandedness):
        strandedness = pe_strandedness
    if (se_strandedness and not pe_strandedness):
        strandedness = se_strandedness

    with open('${sample}.strandedness', 'w') as outfile:
        outfile.write(strandedness or 'undetermined')
    """
}


/*
 * STEP 9.1 - Read the inferred strandedness values
 */
strandedness_results
    .map { sample, strandedness_file ->
        tuple(sample, strandedness_file.text)
    }
    .set {
        rsem_strandedness
    }


/*
 * STEP 9.2 - Summarize sample strandedness values
 */
process summarize_strandedness {

    label 'python'

    publishDir "${params.outdir}/Strandedness"

    input:
    file 'samples/*' from strandedness_summary.collect()

    output:
    file 'summary.tsv'

    """
    #!/usr/bin/env python3
    import csv
    import os

    samples = dict()
    for f in os.listdir('samples'):
        sample = os.path.splitext(f)[0]
        with open("samples/{}".format(f), 'r') as infile:
            samples[sample] = infile.readlines().pop()

    with open('summary.tsv', 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\\t', lineterminator='\\n')
        writer.writerow(['sample', 'strandedness'])
        writer.writerows([s, samples[s]] for s in sorted(samples))
    """
}


/*
 * STEP 10 - RSEM
 */
star_transcriptome_alignments
    .into {
        rsem_transcriptome_alignments
        rsem_input_alignments
    }


/*
 * STEP 10.1 - Create a parameters file
 */
process rsem_mparams {
    tag { sample }

    publishDir "${params.outdir}/RSEM/samples/${sample}/stats"

    input:
    set sample, strandedness from rsem_strandedness

    output:
    set sample, file("*.mparams") into rsem_mparams_files

    script:
    def rsem_mate_fragment_length_min = params.rsem_fragment_length_min
    def rsem_mate_fragment_length_max = params.rsem_fragment_length_max

    def forward_probability = '0.5'
    if (strandedness == 'forward') {
        forward_probability = '1.0'
    }
    if (strandedness == 'reverse') {
        forward_probability = '0.0'
    }

    """
    cat << EOF > "${sample}.mparams"
    ${params.rsem_fragment_length_min} ${params.rsem_fragment_length_max}
    ${forward_probability}
    ${params.rsem_estimate_rspd ? 1 : 0}
    ${params.rsem_num_rspd_bins}
    ${rsem_mate_fragment_length_min} ${rsem_mate_fragment_length_max}
    ${params.rsem_fragment_length_mean} ${params.rsem_fragment_length_sd}
    ${params.rsem_seed_length}
    EOF
    """
}


/*
 * STEP 10.2 - Parse alignments
 */
process rsem_parse_alignments {
    tag { sample }

    label 'rsem'

    publishDir "${params.outdir}/RSEM/samples/${sample}/reads", pattern: '*.fq'
    publishDir "${params.outdir}/RSEM/samples/${sample}/stats", pattern: '*.{cnt,data,omit}'

    input:
    set sample, file(bam) from rsem_transcriptome_alignments
    file "${RSEM}/*" from rsem_index_files.mix(rsem_seq_index).collect()

    output:
    set sample, file("*_alignable{,_1,_2}.fq") into rsem_alignable_reads
    set sample, file("*.{cnt,dat,fq}") into rsem_parsed_alignments
    set sample, file("*.cnt") into rsem_alignment_stats
    set sample, file("*.omit") into rsem_omit_info

    script:
    def read_type = params.paired_end ? '3' : '1'

    """
    rsem-parse-alignments \\
        "${RSEM}/${RSEM}" \\
        "${sample}" \\
        "${sample}" \\
        "${bam}" \\
        "${read_type}"
    """
}


/*
 * STEP 10.3 - Build read indexes
 */
process rsem_build_read_index {
    tag { sample }

    label 'rsem'

    publishDir "${params.outdir}/RSEM/samples/${sample}/reads"

    input:
    set sample, file(fastqs) from rsem_alignable_reads

    output:
    set sample, file("*.ridx") into rsem_read_indexes

    script:
    if (params.paired_end) {

        def (fq1, fq2) = fastqs

        """
        rsem-build-read-index \\
            32 \\
            1 \\
            0 \\
            "${fq1}" \\
            "${fq2}"
        """

    } else {

        """
        rsem-build-read-index \\
            32 \\
            1 \\
            0 \\
            "${fastqs}"
        """
    }
}


/*
 * STEP 10.4 - Create the inputs for EM
 */
rsem_input_alignments
    .join( rsem_parsed_alignments )
    .join( rsem_read_indexes )
    .join( rsem_mparams_files )
    .map { sample, bam, reads, indexes, mparams ->
        tuple(sample, bam, [reads, indexes, mparams].flatten())
    }
    .set { rsem_em_inputs }


/*
 * STEP 10.5 - Run EM
 */
process rsem_run_em {
    tag { sample }

    label 'rsem'

    publishDir "${params.outdir}/RSEM/samples/${sample}/stats", pattern: '*.{model,theta,ofg}'
    publishDir "${params.outdir}/RSEM/samples/${sample}/bam", pattern: '*.bam'

    input:
    set sample, file(bam), file("*") from rsem_em_inputs
    file "${RSEM}/*" from rsem_index_files.mix(rsem_seq_index).collect()

    output:
    set sample, file("*.{model,theta,ofg}") into rsem_em_outputs, rsem_model
    set sample, file("*.transcript.bam") into transcript_bams
    set sample, file("*.gene_res") into rsem_em_gene_res
    set sample, file("*.iso_res") into rsem_em_iso_res

    script:
    def read_type = params.paired_end ? '3' : '1'

    """
    rsem-run-em \\
        "${RSEM}/${RSEM}" \\
        "${read_type}" \\
        "${sample}" \\
        "${sample}" \\
        "${sample}" \\
        -p "${task.cpus}" \\
        -b "${bam}" 0 \\
        --gibbs-out
    """
}


/*
 * STEP 10.5.1 - Run Gibbs sampling
 */
process rsem_run_gibbs {
    tag { sample }

    label 'rsem'

    input:
    set sample, file("*"), file("*") from rsem_em_outputs.join(rsem_omit_info)
    file "${RSEM}/*" from rsem_index_files.mix(rsem_seq_index).collect()

    output:
    set sample, file("${sample}.countvectors*") into rsem_count_vectors
    set sample, file("*.gene_res") into rsem_gibbs_gene_res
    set sample, file("*.iso_res") into rsem_gibbs_iso_res

    """
    rsem-run-gibbs \\
        "${RSEM}/${RSEM}" \\
        "${sample}" \\
        "${sample}" \\
        "${params.gibbs_burnin}" \\
        "${params.gibbs_num_samples}" \\
        "${params.gibbs_sampling_gap}" \\
        -p "${task.cpus}"
    """
}


/*
 * STEP 10.5.2 - Calculate credibility intervals
 */
process rsem_calc_ci {
    tag { sample }

    label 'rsem'

    input:
    set sample, file("*"), file("*") from rsem_count_vectors.join(rsem_model)
    file "${RSEM}/*" from rsem_index_files.mix(rsem_seq_index).collect()

    output:
    set sample, file("*.gene_res") into rsem_ci_gene_res
    set sample, file("*.iso_res") into rsem_ci_iso_res

    """
    rsem-calculate-credibility-intervals \\
        "${RSEM}/${RSEM}" \\
        "${sample}" \\
        "${sample}" \\
        "${params.ci_credibility_level}" \\
        "${params.gibbs_num_samples}" \\
        "${params.ci_num_samples_per_count_vector}" \\
        "${task.memory.toMega() - (task.attempt * 1024)}" \\
        -p "${task.cpus}"
    """
}


/*
 * STEP 10.5.3 - Collect results
 */
process collect_gene_results {
    tag { sample }

    label 'python'

    input:
    set sample, file('em.gene_res'), file('gibbs.gene_res'), file('ci.gene_res') \
        from rsem_em_gene_res.join(rsem_gibbs_gene_res).join(rsem_ci_gene_res)

    output:
    set sample, file("${sample}.gene.results") into rsem_gene_results

    """
    #!/usr/bin/env python3
    import csv
    import fileinput
    import sys

    csv.field_size_limit(sys.maxsize)

    em_fieldnames = [
        "gene_id",
        "transcript_id(s)",
        "length",
        "effective_length",
        "expected_count",
        "TPM",
        "FPKM",
    ]

    gibbs_fieldnames = [
        "posterior_mean_count",
        "posterior_standard_deviation_of_count",
        "pme_TPM",
        "pme_FPKM",
    ]

    ci_fieldnames = [
        "TPM_ci_lower_bound",
        "TPM_ci_upper_bound",
        "TPM_coefficient_of_quartile_variation",
        "FPKM_ci_lower_bound",
        "FPKM_ci_upper_bound",
        "FPKM_coefficient_of_quartile_variation",
    ]

    files = ["em.gene_res", "gibbs.gene_res", "ci.gene_res"]
    fieldnames = em_fieldnames + gibbs_fieldnames + ci_fieldnames

    with fileinput.input(files=files) as infile, open('${sample}.gene.results', 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\\t', lineterminator='\\n')
        writer.writerow(fieldnames)
        writer.writerows(zip(*csv.reader(infile, delimiter='\\t')))
    """
}


process collect_isoform_results {
    tag { sample }

    label 'python'

    input:
    set sample, file('em.iso_res'), file('gibbs.iso_res'), file('ci.iso_res') \
        from rsem_em_iso_res.join(rsem_gibbs_iso_res).join(rsem_ci_iso_res)

    output:
    set sample, file("${sample}.isoforms.results") into rsem_isoform_results

    """
    #!/usr/bin/env python3
    import csv
    import fileinput
    import sys

    csv.field_size_limit(sys.maxsize)

    em_fieldnames = [
        "transcript_id",
        "gene_id",
        "length",
        "effective_length",
        "expected_count",
        "TPM",
        "FPKM",
        "IsoPct",
    ]

    gibbs_fieldnames = [
        "posterior_mean_count",
        "posterior_standard_deviation_of_count",
        "pme_TPM",
        "pme_FPKM",
        "IsoPct_from_pme_TPM",
    ]

    ci_fieldnames = [
        "TPM_ci_lower_bound",
        "TPM_ci_upper_bound",
        "TPM_coefficient_of_quartile_variation",
        "FPKM_ci_lower_bound",
        "FPKM_ci_upper_bound",
        "FPKM_coefficient_of_quartile_variation",
    ]

    files = ["em.iso_res", "gibbs.iso_res", "ci.iso_res"]
    fieldnames = em_fieldnames + gibbs_fieldnames + ci_fieldnames

    with fileinput.input(files=files) as infile, open('${sample}.isoforms.results', 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\\t', lineterminator='\\n')
        writer.writerow(fieldnames)
        writer.writerows(zip(*csv.reader(infile, delimiter='\\t')))
    """
}


/*
 * STEP 10.5.4 - Annotate results
 */
process annotate_gene_results {
    tag { sample }

    label 'python'

    publishDir "${params.outdir}/RSEM/samples/${sample}/counts"

    input:
    set sample, file(results) from rsem_gene_results
    file genes from gene_annotation
    file transcripts from transcript_annotation

    output:
    file "${sample}.gene.results.tsv" into annotated_gene_results

    """
    #!/usr/bin/env python3
    import csv

    with open('${genes}', 'r') as genes_file:
        gene_reader = csv.DictReader(genes_file, delimiter='\\t')
        gene_fieldnames = gene_reader.fieldnames
        genes = { row['gene_id']: row for row in gene_reader }

    with open('${transcripts}', 'r') as transcripts_file:
        transcript_reader = csv.DictReader(transcripts_file, delimiter='\\t')
        transcript_fieldnames = transcript_reader.fieldnames
        transcripts = { row['transcript_id']: row for row in transcript_reader }

    with open("${results}", 'r') as gene_results, open('${sample}.gene.results.tsv', 'w') as outfile:
        results_reader = csv.DictReader(gene_results,  delimiter='\\t')
        fieldnames = list(results_reader.fieldnames)

        gene_id_idx = fieldnames.index('gene_id')
        fieldnames[gene_id_idx:gene_id_idx+1] = gene_fieldnames

        transcripts_id_idx = fieldnames.index('transcript_id(s)')
        fieldnames[transcripts_id_idx:transcripts_id_idx+1] = ['{}(s)'.format(k) for k in transcript_fieldnames]

        results_writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\\t', lineterminator='\\n')
        results_writer.writeheader()

        for row in results_reader:

            gene = genes[row['gene_id']]
            row['gene_name'] = gene['gene_name']
            row['gene_type'] = gene['gene_type']

            ids = row['transcript_id(s)'].split(',')
            row['transcript_name(s)'] = ','.join([transcripts[x]['transcript_name'] for x in ids])
            row['transcript_type(s)'] = ','.join([transcripts[x]['transcript_type'] for x in ids])

            results_writer.writerow(row)
    """
}


process annotate_isoform_results {
    tag { sample }

    label 'python'

    publishDir "${params.outdir}/RSEM/samples/${sample}/counts"

    input:
    set sample, file(results) from rsem_isoform_results
    file genes from gene_annotation
    file transcripts from transcript_annotation

    output:
    file "${sample}.isoform.results.tsv" into annotated_isoform_results

    """
    #!/usr/bin/env python3
    import csv

    with open('${genes}', 'r') as genes_file:
        gene_reader = csv.DictReader(genes_file, delimiter='\\t')
        gene_fieldnames = gene_reader.fieldnames
        genes = { row['gene_id']: row for row in gene_reader }

    with open('${transcripts}', 'r') as transcripts_file:
        transcript_reader = csv.DictReader(transcripts_file, delimiter='\\t')
        transcript_fieldnames = transcript_reader.fieldnames
        transcripts = { row['transcript_id']: row for row in transcript_reader }

    with open("${results}", 'r') as isoform_results, open('${sample}.isoform.results.tsv', 'w') as outfile:
        results_reader = csv.DictReader(isoform_results,  delimiter='\\t')
        fieldnames = list(results_reader.fieldnames)

        transcripts_id_idx = fieldnames.index('transcript_id')
        fieldnames[transcripts_id_idx:transcripts_id_idx+1] = transcript_fieldnames

        gene_id_idx = fieldnames.index('gene_id')
        fieldnames[gene_id_idx:gene_id_idx+1] = gene_fieldnames

        results_writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\\t', lineterminator='\\n')
        results_writer.writeheader()

        for row in results_reader:

            transcript = transcripts[row['transcript_id']]
            row['transcript_name'] = transcript['transcript_name']
            row['transcript_type'] = transcript['transcript_type']

            gene = genes[row['gene_id']]
            row['gene_name'] = gene['gene_name']
            row['gene_type'] = gene['gene_type']

            results_writer.writerow(row)
    """
}


/*
 * STEP 10.6 - Summarize results
 */
process summarize_genes {

    label 'python'

    publishDir "${params.outdir}/RSEM"

    input:
    file 'samples/*' from annotated_gene_results.collect()

    output:
    file "*.tsv"

    """
    #!/usr/bin/env python3
    import csv
    import os

    from contextlib import ExitStack

    output_types = [
        'expected_count',
        'FPKM',
        'TPM',
    ]

    output_fieldnames = [
        'gene_id',
        'gene_name',
        'gene_type',
        'transcript_id(s)',
        'transcript_name(s)',
        'transcript_type(s)',
    ]

    input_filenames = sorted(os.listdir('samples'))
    samples = [fn[:-len('.gene.results.tsv')] for fn in input_filenames]

    with ExitStack() as stack:
        output_files = [stack.enter_context(open("rsem.genes.{}.tsv".format(t), 'w'))
                        for t in output_types ]
        writers = {t: csv.DictWriter(f, fieldnames=output_fieldnames + samples, delimiter='\\t')
                   for t, f in zip(output_types, output_files)}

        for writer in writers.values():
            writer.writeheader()

        input_files = [stack.enter_context(open("samples/{}".format(fn))) for fn in input_filenames]
        readers = [csv.DictReader(input_file, delimiter='\\t') for input_file in input_files]
        rows = [next(row, None) for row in readers]

        while not all(v is None for v in rows):
            output = dict()

            for fieldname in output_fieldnames:
                fields = set([row[fieldname] for row in rows])
                assert len(fields) == 1
                output[fieldname] = fields.pop()

            output_dicts = {output_type: output.copy() for output_type in output_types}

            for sample, row in zip(samples, rows):
                for output_type in output_types:
                    output_dicts[output_type][sample] = row[output_type]

            for output_type, writer in writers.items():
                writer.writerow(output_dicts[output_type])

            rows = [next(row, None) for row in readers]
    """
}


process summarize_isoforms {

    label 'python'

    publishDir "${params.outdir}/RSEM"

    input:
    file 'samples/*' from annotated_isoform_results.collect()

    output:
    file "*.tsv"

    """
    #!/usr/bin/env python3
    import csv
    import os

    from contextlib import ExitStack

    output_types = [
        'expected_count',
        'FPKM',
        'TPM',
    ]

    output_fieldnames = [
        'transcript_id',
        'transcript_name',
        'transcript_type',
        'gene_id',
        'gene_name',
        'gene_type',
    ]

    input_filenames = sorted(os.listdir('samples'))
    samples = [fn[:-len('.isoform.results.tsv')] for fn in input_filenames]

    with ExitStack() as stack:
        output_files = [stack.enter_context(open("rsem.isoforms.{}.tsv".format(t), 'w'))
                        for t in output_types ]
        writers = {t: csv.DictWriter(f, fieldnames=output_fieldnames + samples, delimiter='\\t')
                   for t, f in zip(output_types, output_files)}

        for writer in writers.values():
            writer.writeheader()

        input_files = [stack.enter_context(open("samples/{}".format(fn))) for fn in input_filenames]
        readers = [csv.DictReader(input_file, delimiter='\\t') for input_file in input_files]
        rows = [next(row, None) for row in readers]

        while not all(v is None for v in rows):
            output = dict()

            for fieldname in output_fieldnames:
                fields = set([row[fieldname] for row in rows])
                assert len(fields) == 1
                output[fieldname] = fields.pop()

            output_dicts = {output_type: output.copy() for output_type in output_types}

            for sample, row in zip(samples, rows):
                for output_type in output_types:
                    output_dicts[output_type][sample] = row[output_type]

            for output_type, writer in writers.items():
                writer.writerow(output_dicts[output_type])

            rows = [next(row, None) for row in readers]
    """
}


/*
 * STEP 10.7 - Coordinate-sorted BAM and Wiggle generation
 */
transcript_bams
    .into {
        rsem_tbam2gbam_inputs
        rsem_get_unique_inputs
        rsem_transcript_bams
    }


/*
 * STEP 10.7.1 - Convert the transcript BAM files to genome BAM files
 */
process rsem_tbam2gbam {
    tag { sample }

    label 'rsem'

    publishDir "${params.outdir}/RSEM/samples/${sample}/bam"

    input:
    set sample, file(transcript_bam) from rsem_tbam2gbam_inputs
    file "${RSEM}/*" from rsem_index_files.mix(rsem_seq_index).collect()

    output:
    set sample, val('genome'), file("${sample}.genome.bam") into rsem_genome_bams

    """
    rsem-tbam2gbam \\
        "${RSEM}/${RSEM}" \\
        "${transcript_bam}" \\
        "${sample}.genome.bam" \\
        -p "${task.cpus}"
    """
}


/*
 * STEP 10.7.2 - Extract a subset of unique transcripts
 */
process rsem_get_unique {
    tag { sample }

    label 'rsem'

    publishDir "${params.outdir}/RSEM/samples/${sample}/bam"

    input:
    set sample, file(transcript_bam) from rsem_get_unique_inputs

    output:
    set sample, val('uniq.transcript'), file("${sample}.uniq.transcript.bam") into rsem_uniq_transcript_bams

    """
    rsem-get-unique \\
        "${task.cpus}" \\
        "${transcript_bam}" \\
        "${sample}.uniq.transcript.bam"
    """
}


/*
 * STEP 10.7.3 - Create the tagged inputs for sorting
 */
rsem_transcript_bams
    .map { sample, bam ->
        tuple(sample, 'transcript', bam)
    }
    .mix( rsem_genome_bams )
    .mix( rsem_uniq_transcript_bams )
    .set { rsem_bams }


/*
 * STEP 10.7.4 - Coordinate-sort the transcript and genome BAM files
 */
process rsem_sort_bam {
    tag { "${sample}.${tag}" }

    label 'samtools'

    publishDir "${params.outdir}/RSEM/samples/${sample}/bam"

    input:
    set sample, tag, file(bam) from rsem_bams

    output:
    set sample, tag, file("${sample}.${tag}.sorted.bam") into sorted_bams, rsem_sorted_bams

    """
    samtools sort \\
        -o "${sample}.${tag}.sorted.bam" \\
        -T "${sample}.${tag}.sorted" \\
        -@ "${task.cpus - 1}" \\
        "${bam}"
    """
}


/*
 * STEP 10.7.5 - Index the coordinate-sorted transcript and genome BAM files
 */
process rsem_index_bam {
    tag { "${sample}.${tag}" }

    label 'samtools'

    publishDir "${params.outdir}/RSEM/samples/${sample}/bam"

    input:
    set sample, tag, file(bam) from sorted_bams

    output:
    set sample, tag, file("${bam}.bai") into rsem_bam_indexes

    """
    samtools index "${bam}"
    """
}


/*
 * STEP 10.7.6 - Create the inputs for Wiggle generation and read-depth calculations
 */
rsem_sorted_bams
    .join( rsem_bam_indexes, by: [0,1] )
    .map { sample, tag, bam, idx ->
        tuple(sample, tag, [bam, idx])
    }
    .tap { rsem_bam2wig_inputs }
    .filter { sample, tag, bam_files ->
        tag in ['transcript', 'uniq.transcript']
    }
    .set {
        rsem_bam2readdepth_inputs
    }


/*
 * STEP 10.7.7 - Create a Wiggle file for each coordinate-sorted BAM
 */
process rsem_bam2wig {
    tag { "${sample}.${tag}" }

    label 'rsem'

    publishDir "${params.outdir}/RSEM/samples/${sample}/wig"

    input:
    set sample, tag, file(bam_files) from rsem_bam2wig_inputs

    output:
    set sample, tag, file("${sample}.${tag}.wig") into rsem_wig_output

    script:
    def (bam, idx) = bam_files

    """
    rsem-bam2wig \\
        "${bam}" \\
        "${sample}.${tag}.wig" \\
        "${sample}.${tag}"
    """
}


/*
 * STEP 10.7.7.1 - Add an entry for "uniq.transcript" to the list of chrom sizes
 */
tagged_chrom_sizes
    .tap { tagged_genome_sizes }
    .filter { tag, sizes ->
        tag == 'transcript'
    }
    .map { tag, sizes ->
        tuple("uniq.${tag}".toString(), sizes)
    }
    .mix( tagged_genome_sizes )
    .set { rsem_chrom_sizes }


/*
 * STEP 10.7.7.2 - Combine the wig files with the chrom sizes files
 */
rsem_wig_output
    .map { sample, tag, wig ->
        tuple(tag, sample, wig)
    }
    .combine( rsem_chrom_sizes, by:0 )
    .set { wig2bigwig_inputs  }


/*
 * STEP 10.7.7.3 - Convert wiggle files to bigwig format
 */
process rsem_wig2bigwig {
    tag { "${sample}.${tag}" }

    label 'ucsc_wigtobigwig'

    publishDir "${params.outdir}/RSEM/samples/${sample}/wig"

    input:
    set tag, sample, file(wig), file(chrom_sizes) from wig2bigwig_inputs

    output:
    file "*.bw"

    """
    wigToBigWig \\
        "${wig}" \\
        "${chrom_sizes}" \\
        "${sample}.${tag}.bw"
    """
}


/*
 * STEP 10.7.8 - Generate readdepth output for each transcript BAM
 */
process rsem_bam2readdepth {
    tag { "${sample}.${tag}" }

    label 'rsem'

    publishDir "${params.outdir}/RSEM/samples/${sample}/readdepth"

    input:
    set sample, tag, file(bam_files) from rsem_bam2readdepth_inputs

    output:
    set sample, tag, file("${sample}.${tag}.readdepth") into rsem_readdepth_output

    script:
    def (bam, idx) = bam_files

    """
    rsem-bam2readdepth \\
        "${bam}" \\
        "${sample}.${tag}.readdepth"
    """
}


/*
 * STEP 11 - Run MultiQC
 */
process multiqc {

    label 'multiqc'

    publishDir "${params.outdir}/MultiQC"

    input:
    file config from multiqc_cfg
    file 'fastqc-raw/*' from fastqc_raw_results.collect()
    file 'cutadapt/*' from cutadapt_logs.collect()
    file 'fastqc-trimmed/*' from fastqc_trimmed_results.collect()
    file 'star/*' from star_logs.collect()
    file 'picard/markduplicates/*' from mark_duplicates_metrics.collect()
    file 'rseqc/bam_stat/*' from rseqc_bam_stat_results.collect()
    file 'rseqc/gene_body_coverage/*' from rseqc_gene_body_coverage_results.collect()
    file 'rseqc/infer_experiment/*' from rseqc_infer_experiment_results.collect()
    file 'rseqc/inner_distance/*' from rseqc_inner_distance_results.collect()
    file 'rseqc/junction_annotation/*' from rseqc_junction_annotation_results.collect()
    file 'rseqc/junction_saturation/*' from rseqc_junction_saturation_results.collect()
    file 'rseqc/read_distribution/*' from rseqc_read_distribution_results.collect()
    file 'rseqc/read_duplication/*' from rseqc_read_duplication_results.collect()
    file 'rseqc/read_gc/*' from rseqc_read_gc_results.collect()
    file 'rnaseqc/*' from rnaseqc_results.collect()
    file 'rsem/*' from rsem_alignment_stats.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    """
    multiqc \\
        --config "${config}" \\
        -m fastqc \\
        -m cutadapt \\
        -m star \\
        -m picard \\
        -m rseqc \\
        -m rna_seqc \\
        -m rsem \\
        .
    """
}


/* Workflow */

workflow.onComplete {

    log.info "Workflow completed at: ${workflow.complete}"
    log.info "Time taken: ${workflow.duration}"
    log.info "Execution status: ${workflow.success ? 'success' : 'failed'}"
    log.info "Output directory: ${params.outdir}"
}

workflow.onError {

    log.info "Execution halted: ${workflow.errorMessage}"
}


/* Functions */

def usage() {

    log.info"""
    Usage:
        nextflow run ampatchlab/nf-rnaseq [options]


    Required arguments:

        --csv FILE
            Comma-separated list of sample and readgroup inputs

        --genome STR
            Reference genome name [Either: ${defaults.genomes.keySet().join(", ")}; Default: ${defaults.genome ?: null}]


    Optional arguments:

        --paired_end
            Expect entries for 'fastq1' and 'fastq2' in the input CSV

        --rgid_sep STR
            The separator used to create unique input readgroup IDs [Default: ${defaults.rgid_sep}]

        --adapters STR
            The adapters to trim [Either: ${defaults.seq_adapters.keySet().join(", ")}; Default: ${defaults.adapters}]


    Reference genome options:

        --fasta FILE
            Override the reference genome FASTA with FILE [Default: ${defaults.fasta ?: null}]

        --gtf FILE
            Override the reference genome GTF with FILE [Default: ${defaults.gtf ?: null}]


    STAR genome generate options:

        --star_genome_chr_bin_n_bits INT
            Size of the bins for genome storage [Default: ${defaults.star_genome_chr_bin_n_bits}]

        --star_genome_sa_index_n_bases INT
            Length (bases) of the SA pre-indexing string [Default: ${defaults.star_genome_sa_index_n_bases}]

        --star_sjdb_overhang INT
            Length of the donor/acceptor sequence on each side of the junctions [Default: ${defaults.star_sjdb_overhang}]

    Cutadapt options:

        --cutadapt_r1_adapter STR
            Sequence of the R1 adapter [Default: ${defaults.seq_adapters[params.adapters].r1}]

        --cutadapt_r2_adapter STR
            Sequence of the R2 adapter [Default: ${defaults.seq_adapters[params.adapters].r2}]

        --cutadapt_base_qual_cutoff [INT,]INT
            Trim low-quality bases from each read [Default: ${defaults.cutadapt_base_qual_cutoff}]

        --cutadapt_min_read_length INT[:INT]
            Discard reads shorter than INT [Default: ${defaults.cutadapt_base_qual_cutoff}]


    RSEM mparams:

        --rsem_fragment_length_min INT
            Minimum read/insert length allowed [default: ${defaults.rsem_fragment_length_min}]

        --rsem_fragment_length_max INT
            Maximum read/insert length allowed [default: ${defaults.rsem_fragment_length_max}]

        --rsem_fragment_length_mean INT (single-end data only)
            The mean of the fragment length distribution, which is assumed to
            be Gaussian. A value of -1, disables the use of the fragment length
            distribution [Default: ${defaults.rsem_fragment_length_mean}]

        --rsem_fragment_length_sd INT (single-end data only)
            The standard deviation of the fragment length distribution, which
            is assumed to be Gaussian. A value of 0, assumes that all fragments
            are of the same length [Default: ${defaults.rsem_fragment_length_sd}]

        --rsem_estimate_rspd
            Set this option if you want to estimate the read start position
            distribution (RSPD) from data. Otherwise, RSEM will use a uniform
            RSPD

        --rsem_num_rspd_bins INT
            Number of bins in the RSPD. Only relevant when '--rsem_estimate_rspd'
            is specified. The default value is recommended. [Default: ${defaults.rsem_num_rspd_bins}]

        --rsem_seed_length INT
            Seed length [Default: ${defaults.rsem_seed_length}]


    RSEM Gibbs options:

        --gibbs_burnin INT
            The number of burn-in rounds for RSEM's Gibbs sampler. Each round
            passes over the entire data set once. If RSEM can use multiple
            threads, multiple Gibbs samplers will start at the same time and
            all samplers share the same burn-in number [Default: ${defaults.gibbs_burnin}]

        --gibbs_num_samples INT
            The total number of count vectors RSEM will collect from its Gibbs
            samplers [Default: ${defaults.gibbs_num_samples}]

        --gibbs_sampling_gap INT
            The number of rounds between two succinct count vectors RSEM collects.
            If the count vector after round N is collected, the count vector after
            round N + INT will also be collected [Default: ${defaults.gibbs_sampling_gap}]


    RSEM CalcCI options:

        --ci_credibility_level FLOAT
            The credibility level for credibility intervals [Default: ${defaults.ci_credibility_level}]

        --ci_num_samples_per_count_vector INT
            The number of read generating probability vectors sampled per sampled
            count vector. The crebility intervals are calculated by first sampling
            P(C | D) and then sampling P(Theta | C) for each sampled count vector.
            This option controls how many Theta vectors are sampled per sampled
            count vector [Default: ${defaults.ci_num_samples_per_count_vector}]


    MultiQC options:

        --multiqc_config FILE
            MultiQC YAML config file [Default: ${defaults.multiqc_config}]


    Output options:

        --refdir DIR
            Path where the reference index files will be saved [Default: ${defaults.refdir}]

        --outdir DIR
            Path where the results will be saved [Default: ${defaults.outdir}]


    Report options

        --execution_report STR
            Name of the Nextflow execution report to generate [Default: ${defaults.execution_report}]

        --trace_report STR
            Name of the Nextflow trace report to generate [Default: ${defaults.trace_report}]

        --timeline_report STR
            Name of the Nextflow timeline report to generate [Default: ${defaults.timeline_report}]

        --flowchart STR
            Name of the Nextflow flowchart to generate [Default: ${defaults.flowchart}]


    AWS Batch options

        --aws_queue STR
            AWS Batch JobQueue definition [Default: ${defaults.aws_queue}]

        --aws_region STR
            AWS Region definition [Default: ${defaults.aws_region}]


    Standard options:

        --help
            Show this message and exit

        --version
            Show the pipeline version and exit
    """.stripIndent()
}

def die() {
    usage()
    exit 1
}

def check_params() {

    // standard options

    if (params.help) {
        usage()
        exit 0
    }

    if (params.version) {
        log.info(workflow.manifest.version)
        exit 0
    }


    // required arguments

    if (!params.csv) {
        log.error("A list of samples and readgroups is required. Please use the `--csv` option.")
        die()
    }

    if (file(params.csv).getExtension() != "csv") {
        log.error("Readgroup input file `${params.csv}` must be a CSV file with the '.csv' extension.")
        die()
    }


    // reference genome options

    if (params.genome && !(params.genome in params.genomes)) {
        log.error("Unknown `--genome` entry: `${params.genome}`")
        die()
    }

    params.fasta = params.genome ? params.genomes[ params.genome ].fasta : null
    params.gtf = params.genome ? params.genomes[ params.genome ].gtf : null

    if (!params.genome && !params.fasta && !params.gtf) {
        log.error("A reference genome is required. Please use the `--genome` option.")
        die()
    }

    if (!params.fasta) {
        log.error("A reference FASTA file is required. Please use the `--fasta` option.")
        die()
    }

    if (!params.gtf) {
        log.error("A reference GTF file is required. Please use the `--gtf` option.")
        die()
    }


    // STAR genome generate options

    if (!(params.star_genome_chr_bin_n_bits.toString().isInteger())) {
        log.error("Unknown `--star_genome_chr_bin_n_bits` entry: `${params.star_genome_chr_bin_n_bits}`")
        die()
    }

    if (!(params.star_genome_sa_index_n_bases.toString().isInteger())) {
        log.error("Unknown `--star_genome_sa_index_n_bases` entry: `${params.star_genome_sa_index_n_bases}`")
        die()
    }

    if (!(params.star_sjdb_overhang.toString().isInteger())) {
        log.error("Unknown `--star_sjdb_overhang` entry: `${params.star_sjdb_overhang}`")
        die()
    }


    // Cutadapt options

    if (params.adapters && !(params.adapters in params.seq_adapters)) {
        log.error("Unknown `--adapters` entry: `${params.adapters}`")
        die()
    }

    params.cutadapt_r1_adapter = params.adapters ? params.seq_adapters[ params.adapters ].r1 : null
    params.cutadapt_r2_adapter = params.adapters ? params.seq_adapters[ params.adapters ].r2 : null

    if (params.cutadapt_r1_adapter && !(params.cutadapt_r1_adapter instanceof String)) {
        log.error("Unknown `--cutadapt_r1_adapter` entry: `${params.cutadapt_r1_adapter}`")
        die()
    }

    if (params.cutadapt_r2_adapter && !(params.cutadapt_r2_adapter instanceof String)) {
        log.error("Unknown `--cutadapt_r2_adapter` entry: `${params.cutadapt_r2_adapter}`")
        die()
    }

    if (!(params.cutadapt_base_qual_cutoff.toString().isInteger())) {

        if (params.cutadapt_base_qual_cutoff.toString().contains(',')) {
            def (five_prime_cutoff, three_prime_cutoff) = params.cutadapt_base_qual_cutoff.split(',', 2)

            if (!five_prime_cutoff.isInteger() || !three_prime_cutoff.isInteger()) {
                log.error("Unknown `--cutadapt_base_qual_cutoff` entry: `${params.cutadapt_base_qual_cutoff}`")
                die()
            }
        }
        else {
            log.error("Unknown `--cutadapt_base_qual_cutoff` entry: `${params.cutadapt_base_qual_cutoff}`")
            die()
        }
    }

    if (!(params.cutadapt_min_read_length.toString().isInteger())) {

        if (params.cutadapt_min_read_length.toString().contains(':')) {
            def (r1_min_length, r2_min_length) = params.cutadapt_min_read_length.split(':', 2)

            if (!r1_min_length.isInteger() || !r2_min_length.isInteger()) {
                log.error("Unknown `--cutadapt_min_read_length` entry: `${params.cutadapt_min_read_length}`")
                die()
            }
        }
        else {
            log.error("Unknown `--cutadapt_min_read_length` entry: `${params.cutadapt_min_read_length}`")
            die()
        }
    }


    // RSEM mparams options

    if (!(params.rsem_fragment_length_min.toString().isInteger())) {
        log.error("Unknown `--rsem_fragment_length_min` entry: `${params.rsem_fragment_length_min}`")
        die()
    }

    if (params.rsem_fragment_length_min < 1) {
        log.error("The value specified using `--rsem_fragment_length_min` must be >= 1.")
        die()
    }

    if (!(params.rsem_fragment_length_max.toString().isInteger())) {
        log.error("Unknown `--rsem_fragment_length_max` entry: `${params.rsem_fragment_length_max}`")
        die()
    }

    if (params.rsem_fragment_length_min > params.rsem_fragment_length_max) {
        log.error("The minimum fragment length should be less than or equal to the maximum fragment length.")
        die()
    }

    if (!(params.rsem_fragment_length_mean.toString().isInteger())) {
        log.error("Unknown `--rsem_fragment_length_mean` entry: `${params.rsem_fragment_length_mean}`")
        die()
    }

    if (!(params.rsem_fragment_length_sd.toString().isInteger())) {
        log.error("Unknown `--rsem_fragment_length_sd` entry: `${params.rsem_fragment_length_sd}`")
        die()
    }

    if (!(params.rsem_num_rspd_bins.toString().isInteger())) {
        log.error("Unknown `--rsem_num_rspd_bins` entry: `${params.rsem_num_rspd_bins}`")
        die()
    }

    if (params.rsem_num_rspd_bins < 1) {
        log.error("The value specified using `--rsem_num_rspd_bins` must be >= 1.")
        die()
    }

    if (!(params.rsem_seed_length.toString().isInteger())) {
        log.error("Unknown `--rsem_seed_length` entry: `${params.rsem_seed_length}`")
        die()
    }

    if (params.rsem_seed_length < 5) {
        log.error("The value specified using `--rsem_seed_length` must be >= 5.")
        die()
    }


    // RSEM Gibbs options

    if (!(params.gibbs_burnin.toString().isInteger())) {
        log.error("Unknown `--gibbs_burnin` entry: `${params.gibbs_burnin}`")
        die()
    }

    if (params.gibbs_burnin < 1) {
        log.error("The value specified using `--gibbs_burnin` must be >= 1.")
        die()
    }

    if (!(params.gibbs_num_samples.toString().isInteger())) {
        log.error("Unknown `--gibbs_num_samples` entry: `${params.gibbs_num_samples}`")
        die()
    }

    if (params.gibbs_num_samples < 1) {
        log.error("The value specified using `--gibbs_num_samples` must be >= 1.")
        die()
    }

    if (!(params.gibbs_sampling_gap.toString().isInteger())) {
        log.error("Unknown `--gibbs_sampling_gap` entry: `${params.gibbs_sampling_gap}`")
        die()
    }

    if (params.gibbs_sampling_gap < 1) {
        log.error("The value specified using `--gibbs_sampling_gap` must be >= 1.")
        die()
    }


    // RSEM CalcCI options

    if (!(params.ci_credibility_level.toString().isFloat())) {
        log.error("Unknown `--ci_credibility_level` entry: `${params.ci_credibility_level}`")
        die()
    }

    if (!(params.ci_credibility_level >= 0 && params.ci_credibility_level <= 1)) {
        log.error("The value specified using `--ci_credibility_level` must exist in the range 0 to 1.")
        die()
    }

    if (!(params.ci_num_samples_per_count_vector.toString().isInteger())) {
        log.error("Unknown `--ci_num_samples_per_count_vector` entry: `${params.ci_num_samples_per_count_vector}`")
        die()
    }

    if (params.ci_num_samples_per_count_vector < 1) {
        log.error("The value specified using `--ci_num_samples_per_count_vector` must be >= 1.")
        die()
    }


    // MultiQC options

    if (!params.multiqc_config) {
        log.error("A configuration file for MultiQC is required. Please use the `--multiqc_config` option.")
        die()
    }

    if (file(params.multiqc_config).getExtension() != "yaml") {
        log.error("MultiQC config file `${params.multiqc_config}` must be a YAML file with the '.yaml' extension.")
        die()
    }


    // report options

    if (!params.execution_report.toString().endsWith('.html')) {
        log.error("The filename specified using `--execution_report` must end with '.html'")
        die()
    }

    if (!params.trace_report.toString().endsWith('.txt')) {
        log.error("The filename specified using `--trace_report` must end with '.txt'")
        die()
    }

    if (!params.timeline_report.toString().endsWith('.html')) {
        log.error("The filename specified using `--timeline_report` must end with '.html'")
        die()
    }

    def flowchart_extns = ['.dot', '.html', '.pdf', '.png', '.svg']

    if (!(flowchart_extns.any { params.flowchart.toString().endsWith(it) })) {
        log.error("The filename specified using `--flowchart` must end with one of ${flowchart_extns.join(", ")}")
        die()
    }
}

def validate_input_csv() {

    def csv = new CsvSplitter().target(file(params.csv)).options(header:true)
    def rows = csv.list()

    def fastq_columns = params.paired_end ? ["fastq1", "fastq2"] : ["fastq"]
    def required_columns = ["sample"] + fastq_columns

    required_columns.each { col ->

        if (!csv.columnsHeader.contains(col)) {
            log.error("Readgroup input file `${params.csv}` does not contain a '${col}' column. Exiting.")
            exit 1
        }
    }

    def readgroups = [:].withDefault { [] }
    def fastq_files = [:].withDefault { [] }

    def valid_file_extns = ["fastq", "fq", "fastq.gz", "fq.gz"]

    log.info("Validating ${rows.size()} entries...")

    rows.indexed(1).each { idx, row ->

        log.info("Validating entry ${idx}...")

        if (!row.sample) {
            log.error("Entry ${idx}: Invalid 'sample' value. Exiting.")
            exit 1
        }

        if (row.readgroup in readgroups[row.sample]) {
            log.error("Entry ${idx}: Sample `${row.sample}` must have unique readgroup entries. Exiting.")
            exit 1
        }

        readgroups[row.sample].add(row.readgroup)

        fastq_columns.each { fastq ->

            if (!row."${fastq}") {
                log.error("Entry ${idx}: Invalid '${fastq}' value. Exiting.")
                exit 1
            }

            fq = file(row."${fastq}")
            fq_extn = get_fastq_extn(fq)

            if (!valid_file_extns.contains(fq_extn)) {
                log.error("Entry ${idx}: ${fq} does not have a valid file extension. Must be one of: ${valid_file_extns.join(", ")}")
                exit 1
            }

            if (fq.getName() in fastq_files[row.sample]) {
                log.error("Entry ${idx}: File '${fq.getName()}' cannot be used more than once. Exiting.")
                exit 1
            }

            fastq_files[row.sample].add(fq.getName())
        }
    }

    log.info("Done")
}

def get_fastq_extn(fq) {
    def extn
    switch (fq) {
        case { it.name.endsWith('.fastq') }:
            extn = 'fastq'
            break
        case { it.name.endsWith('.fq') }:
            extn = 'fq'
            break
        case { it.name.endsWith('.fastq.gz') }:
            extn = 'fastq.gz'
            break
        case { it.name.endsWith('.fq.gz') }:
            extn = 'fq.gz'
            break
        default:
            extn = ''
            break
    }
    extn
}
