# Stage 2 - clustering. Replaces the two QIIME runs and the mock / identity tables.
#
# q2-vsearch cluster-features-de-novo only uses the table to give each sequence an
# abundance and then calls `vsearch --cluster_size`. We carry the sample and abundance
# in the read headers instead, call vsearch directly and rebuild the table from the
# .uc membership files.

POOLED = config["subsetted_fasta"] if START == "subsetted_fasta" else f"{WORK}/pooled_reads.fasta"
ID97_ABUNDANCE = config["cluster97_order"]


if config.get("dereplicate", True):

    rule derep:
        input:
            POOLED,
        output:
            fasta=f"{WORK}/derep.fasta",
            uc=f"{WORK}/derep.uc",
        conda:
            "../envs/vsearch.yaml"
        log:
            f"{WORK}/logs/derep.log",
        shell:
            "vsearch --derep_fulllength {input} --output {output.fasta} --uc {output.uc} "
            "--sizeout --fasta_width 0 --minseqlength 1 &> {log}"

else:
    # Exact legacy behaviour (validation only, slow on large jobs): no dereplication, so every
    # read has abundance 1 exactly as in the old mock table. Each read is its own "seed".

    rule derep:
        input:
            POOLED,
        output:
            fasta=f"{WORK}/derep.fasta",
            uc=f"{WORK}/derep.uc",
        shell:
            r"""awk -v uc={output.uc} '{{sub(/\r$/,"")}} /^>/{{l=substr($0,2); print ">" l ";size=1"; """
            r"""printf("S\t0\t0\t*\t*\t*\t*\t*\t%s\t*\n", l) > uc; next}} {{print}}' {input} > {output.fasta}"""


rule cluster_100:
    input:
        f"{WORK}/derep.fasta",
    output:
        centroids=f"{WORK}/cluster_100_centroids.fasta",
        uc=f"{WORK}/cluster_100.uc",
    params:
        id=config["cluster_id_100"],
    threads: config["threads"]
    conda:
        "../envs/vsearch.yaml"
    log:
        f"{WORK}/logs/cluster_100.log",
    shell:
        "vsearch --cluster_size {input} --id {params.id} --sizein --xsize --qmask none "
        "--centroids {output.centroids} --uc {output.uc} --fasta_width 0 --minseqlength 1 "
        "--threads {threads} &> {log}"


rule otu_table_100:
    """Long-format 100% OTU table (OTU, sample, count) from the two .uc files."""
    input:
        derep_uc=f"{WORK}/derep.uc",
        cluster_uc=f"{WORK}/cluster_100.uc",
        centroids=f"{WORK}/cluster_100_centroids.fasta",
    output:
        long=f"{WORK}/otu_table_100_long.tsv.gz",
        seqs=f"{WORK}/otu_100_sequences.tsv",
    conda:
        "../envs/r.yaml"
    log:
        f"{WORK}/logs/otu_table_100.log",
    script:
        "../scripts/otu_table_100.R"


rule filter_100:
    input:
        long=f"{WORK}/otu_table_100_long.tsv.gz",
        seqs=f"{WORK}/otu_100_sequences.tsv",
    output:
        table=f"{OUT}/{NAME}_otu_table_100_greater_than_2.csv",
        counts=f"{OUT}/{NAME}_otu_table_100_filtered_counts.csv",
        fasta=f"{OUT}/{NAME}_rep_sequences_100_greater_than_2.fasta",
        fasta97=f"{WORK}/for_cluster_97.fasta",
    params:
        min_reads=config["min_total_reads_100"],
        mode=ID97_ABUNDANCE,
    conda:
        "../envs/r.yaml"
    log:
        f"{WORK}/logs/filter_100.log",
    script:
        "../scripts/filter_100.R"


rule cluster_97:
    input:
        f"{WORK}/for_cluster_97.fasta",
    output:
        centroids=f"{WORK}/cluster_97_centroids.fasta",
        uc=f"{WORK}/cluster_97.uc",
    params:
        id=config["cluster_id_97"],
    threads: config["threads"]
    conda:
        "../envs/vsearch.yaml"
    log:
        f"{WORK}/logs/cluster_97.log",
    shell:
        "vsearch --cluster_size {input} --id {params.id} --sizein --xsize --qmask none "
        "--centroids {output.centroids} --uc {output.uc} --fasta_width 0 --minseqlength 1 "
        "--threads {threads} &> {log}"


rule assign_97:
    input:
        table=f"{OUT}/{NAME}_otu_table_100_greater_than_2.csv",
        seqs=f"{WORK}/otu_100_sequences.tsv",
        uc=f"{WORK}/cluster_97.uc",
    output:
        xlsx=f"{OUT}/{NAME}_otu_table_100_with_reps.xlsx",
        fasta=f"{OUT}/{NAME}_final_sequences_to_subtype.fasta",
        tsv=f"{OUT}/{NAME}_final_sequences_to_subtype.tsv",
    conda:
        "../envs/r.yaml"
    log:
        f"{WORK}/logs/assign_97.log",
    script:
        "../scripts/assign_97.R"
