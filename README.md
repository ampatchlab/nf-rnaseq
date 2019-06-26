# ampatchlab/nf-rnaseq

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.1-brightgreen.svg)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

RNA-Seq nextflow pipeline

## Usage

```
Usage:
    nextflow run ampatchlab/nf-rnaseq [options]


Required arguments:

    --csv FILE
        Comma-separated list of sample and readgroup inputs

    --genome STR
        Reference genome name [Either: GRCh38, GRCm38; Default: null]


Optional arguments:

    --paired_end
        Expect entries for 'fastq1' and 'fastq2' in the input CSV

    --rgid_sep STR
        The separator used to create unique input readgroup IDs [Default: .]

    --adapters STR
        The adapters to trim [Either: TruSeq, null; Default: TruSeq]


Reference genome options:

    --fasta FILE
        Override the reference genome FASTA with FILE [Default: null]

    --gtf FILE
        Override the reference genome GTF with FILE [Default: null]


STAR genome generate options:

    --star_genome_chr_bin_n_bits INT
        Size of the bins for genome storage [Default: 18]

    --star_genome_sa_index_n_bases INT
        Length (bases) of the SA pre-indexing string [Default: 14]


Cutadapt options:

    --cutadapt_r1_adapter STR
        Sequence of the R1 adapter [Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]

    --cutadapt_r2_adapter STR
        Sequence of the R2 adapter [Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT]

    --cutadapt_base_qual_cutoff [INT,]INT
        Trim low-quality bases from each read [Default: 20]

    --cutadapt_min_read_length INT[:INT]
        Discard reads shorter than INT [Default: 20]


RSEM mparams:

    --rsem_fragment_length_min INT
        Minimum read/insert length allowed [default: 1]

    --rsem_fragment_length_max INT
        Maximum read/insert length allowed [default: 1000]

    --rsem_fragment_length_mean INT (single-end data only)
        The mean of the fragment length distribution, which is assumed to
        be Gaussian. A value of -1, disables the use of the fragment length
        distribution [Default: -1]

    --rsem_fragment_length_sd INT (single-end data only)
        The standard deviation of the fragment length distribution, which
        is assumed to be Gaussian. A value of 0, assumes that all fragments
        are of the same length [Default: 0]

    --rsem_estimate_rspd
        Set this option if you want to estimate the read start position
        distribution (RSPD) from data. Otherwise, RSEM will use a uniform
        RSPD

    --rsem_num_rspd_bins INT
        Number of bins in the RSPD. Only relevant when '--rsem_estimate_rspd'
        is specified. The default value is recommended. [Default: 20]

    --rsem_seed_length INT
        Seed length [Default: 25]


RSEM Gibbs options:

    --gibbs_burnin INT
        The number of burn-in rounds for RSEM's Gibbs sampler. Each round
        passes over the entire data set once. If RSEM can use multiple
        threads, multiple Gibbs samplers will start at the same time and
        all samplers share the same burn-in number [Default: 200]

    --gibbs_num_samples INT
        The total number of count vectors RSEM will collect from its Gibbs
        samplers [Default: 1000]

    --gibbs_sampling_gap INT
        The number of rounds between two succinct count vectors RSEM collects.
        If the count vector after round N is collected, the count vector after
        round N + INT will also be collected [Default: 1]


RSEM CalcCI options:

    --ci_credibility_level FLOAT
        The credibility level for credibility intervals [Default: 0.95]

    --ci_num_samples_per_count_vector INT
        The number of read generating probability vectors sampled per sampled
        count vector. The crebility intervals are calculated by first sampling
        P(C | D) and then sampling P(Theta | C) for each sampled count vector.
        This option controls how many Theta vectors are sampled per sampled
        count vector [Default: 50]


MultiQC options:

    --multiqc_config FILE
        MultiQC YAML config file [Default: [:]/assets/multiqc_config.yaml]


Output options:

    --inputs DIR
        Path where the input FASTQ files should be localized [Default: ./inputs]

    --outdir DIR
        Path where the results will be saved [Default: ./results]

    --refdir DIR
        Path where the reference index files will be saved [Default: ./reference]


Report options

    --execution_report STR
        Name of the Nextflow execution report to generate [Default: ./reports/execution_report.html]

    --trace_report STR
        Name of the Nextflow trace report to generate [Default: ./reports/trace_report.txt]

    --timeline_report STR
        Name of the Nextflow timeline report to generate [Default: ./reports/timeline_report.html]

    --flowchart STR
        Name of the Nextflow flowchart to generate [Default: ./reports/flowchart.png]


Standard options:

    --help
        Show this message and exit

    --version
        Show the pipeline version and exit
```

## Inputs

For paired-end data, the input CSV must have the following required columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq1: Absolute path of the 'R1' FASTQ file (required)
 * fastq2: Absolute path of the 'R2' FASTQ file (required)

For single-end data, the CSV must have the following columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq: Absolute path of the FASTQ file (required)

If a particular sample has multiple FASTQ files (or pairs of FASTQ files), then these may
be specified on additional lines with a unique readgroup identifier. All readgroups belonging
to a particular sample will be aligned and merged using STAR.

## Notes

### GRCh37

Input FASTA files must be sorted by karyotypic order or RNA-SeQC v1.1.8 (GATK) will complain.
Unfortunately, the current GRCh37 assembly available from GENCODE is sorted lexicographically.
This appears to be an artifact of the liftover process. To use the GRCh37 reference, please
re-sort the FASTA file using the instructions below. It is hoped that this issue is resolved
in later releases.

1. Obtain the curent GRCh37 FASTA and GTF files:
```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gtf.gz
```

2. Split the GRCh37 FASTA into separate chromosomes/contigs:
```
$ zcat GRCh37.primary_assembly.genome.fa.gz | awk '/^>/ { f=substr($1,2) ".fa" } { print > f }'
```

3. Replace the existing assembly with an ordered set of chromosomes/contigs:
```
$ cat chr{{1..22},X,Y,M}.fa GL*.fa | gzip > GRCh37.primary_assembly.genome.fa.gz
```

4. Clean up:
```
$ rm *.fa
```

### RSEM transcript plots

We use the `rsem-gen-transcript-plots` Rscript to create PDF output of transcripts plots
for a list of gene IDs for a particular sample. The list of gene IDs, must be Ensembl
gene IDs, one per line. In the output, read-depth contributed by unique reads is shown in
black while reads with multiple alignments are shown in red.

1. Navigate to the `readdepth` RSEM results directory:
```
$ cd ./results/RSEM/samples/${sample}/readdepth
```

2. Copy across the raw isoforms results file:
```
$ cp "$(dirname "$(readlink -e ../counts/${sample}.isoform.results.tsv)")/${sample}.isoforms.results" ./
```

3. Create the plots:
```
rsem-gen-transcript-plots "${sample}" gene_ids.txt 0 2 1 output.pdf
```
