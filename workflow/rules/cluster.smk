# Stage 2 - clustering. Replaces the two QIIME runs and the mock / identity tables.
#
# q2-vsearch cluster-features-de-novo only uses the table to give each sequence an
# abundance and then calls `vsearch --cluster_size`. We carry the sample and abundance
# in the read headers instead, call vsearch directly and rebuild the table from the
# .uc membership files.

POOLED = config["subsetted_fasta"] if START == "subsetted_fasta" else f"{WORK}/pooled_reads.fasta"
ID97_ABUNDANCE = config["cluster97_order"]
PRECHIM = f"{WORK}/before_chimera_removal"


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
        table=f"{PRECHIM}/otu_table_100_greater_than_2.csv",
        counts=f"{PRECHIM}/otu_table_100_filtered_counts.csv",
        fasta=f"{PRECHIM}/rep_sequences_100_greater_than_2.fasta",
        fasta97=f"{PRECHIM}/for_cluster_97.fasta",
    params:
        min_reads=config["min_total_reads_100"],
        mode=ID97_ABUNDANCE,
    conda:
        "../envs/r.yaml"
    log:
        f"{WORK}/logs/filter_100.log",
    script:
        "../scripts/filter_100.R"


rule chimera_check:
    """De novo chimera check on the 100% OTUs that passed the read filter, with their total read
    counts as abundance (UCHIME3). Removal is done by the next rule."""
    input:
        table=f"{PRECHIM}/otu_table_100_greater_than_2.csv",
        fasta=f"{PRECHIM}/rep_sequences_100_greater_than_2.fasta",
    output:
        chimeras=f"{S3}/chimeras.fasta",
        sized=f"{S3}/otu_100_sized.fasta",
    conda:
        "../envs/vsearch.yaml"
    log:
        f"{S3}/logs/uchime.log",
    shell:
        r"""
        awk -F, 'NR==FNR{{if(FNR>1)t[$1]=$3;next}} /^>/{{id=substr($0,2);sub(/;.*/,"",id);print ">"id";size="t[id];next}}{{print}}' \
            {input.table} {input.fasta} > {output.sized}
        vsearch --uchime3_denovo {output.sized} --chimeras {output.chimeras} --nonchimeras /dev/null &> {log}
        """


rule remove_chimeras:
    """Drop the OTUs called chimeric (config remove_chimeras) from the table, the read counts and
    the sequences that go on to 97% clustering. With remove_chimeras false nothing is dropped and
    the calls only appear as the possible_chimera flag."""
    input:
        table=f"{PRECHIM}/otu_table_100_greater_than_2.csv",
        counts=f"{PRECHIM}/otu_table_100_filtered_counts.csv",
        fasta=f"{PRECHIM}/rep_sequences_100_greater_than_2.fasta",
        chimeras=f"{S3}/chimeras.fasta",
    output:
        table=f"{OUT}/{NAME}_otu_table_100_greater_than_2.csv",
        counts=f"{OUT}/{NAME}_otu_table_100_filtered_counts.csv",
        fasta=f"{OUT}/{NAME}_rep_sequences_100_greater_than_2.fasta",
        fasta97=f"{WORK}/for_cluster_97.fasta",
    params:
        drop_chimeras=config.get("remove_chimeras", True),
        mode=ID97_ABUNDANCE,
    conda:
        "../envs/r.yaml"
    log:
        f"{WORK}/logs/remove_chimeras.log",
    script:
        "../scripts/remove_chimeras.R"


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
