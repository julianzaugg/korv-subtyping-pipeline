# Stage 3 - subtyping preparation. Produces suggestions, evidence and figures for the expert.
# It never assigns a subtype: the expert completes the assignments file, then Stage 4 runs.

REFS = config["references"]
LIB = config["subtype_library"]
RULES = config["nonfunctional_rules"]
S3 = f"{OUT}/subtyping"
FINAL = f"{OUT}/{NAME}_final_sequences_to_subtype"


rule translate:
    input:
        fasta=f"{FINAL}.fasta",
        refs=REFS,
    output:
        faa=f"{FINAL}.faa",
        tsv=f"{S3}/translation.tsv",
        for_aln=f"{S3}/proteins_for_alignment.faa",
    threads: config["threads"]
    conda:
        "../envs/r.yaml"
    log:
        f"{S3}/logs/translate.log",
    script:
        "../scripts/translate.R"


rule align:
    """Same alignment as the manual workflow: mafft L-INS-i style with gap open 10, extension 1."""
    input:
        f"{S3}/proteins_for_alignment.faa",
    output:
        f"{FINAL}_with_references.aln",
    threads: config["threads"]
    conda:
        "../envs/mafft.yaml"
    log:
        f"{S3}/logs/mafft.log",
    shell:
        "mafft --localpair --reorder --op 10 --ep 1 --maxiterate 1000 --thread {threads} {input} > {output} 2> {log}"


rule tree:
    input:
        f"{FINAL}_with_references.aln",
    output:
        treefile=f"{FINAL}_with_references.treefile",
        iqtree=f"{FINAL}_with_references.iqtree",
    params:
        prefix=lambda wildcards, output: output.treefile.removesuffix(".treefile"),
    threads: config["threads"]
    conda:
        "../envs/iqtree.yaml"
    log:
        f"{S3}/logs/iqtree.log",
    shell:
        "iqtree2 -s {input} -m MFP -mset LG --threads-max {threads} -nt AUTO --prefix {params.prefix} -redo &> {log}"


rule chimera_check:
    input:
        tsv=f"{FINAL}.tsv",
    output:
        chimeras=f"{S3}/chimeras.fasta",
        sized=f"{S3}/final_sized.fasta",
    conda:
        "../envs/vsearch.yaml"
    log:
        f"{S3}/logs/uchime.log",
    shell:
        r"""
        awk -F'\t' 'NR>1{{print ">"$1";size="$4"\n"$5}}' {input.tsv} > {output.sized}
        vsearch --uchime3_denovo {output.sized} --chimeras {output.chimeras} --nonchimeras /dev/null &> {log}
        """


rule library_fasta:
    input:
        LIB,
    output:
        f"{S3}/library.fasta",
    conda:
        "../envs/vsearch.yaml"
    log:
        f"{S3}/logs/library_fasta.log",
    shell:
        r"""awk -F'\t' 'NR>1 && !seen[$1]++{{print ">"$1"\n"$2}}' {input} > {output} 2> {log}"""


rule library_search:
    input:
        query=f"{FINAL}.fasta",
        db=f"{S3}/library.fasta",
    output:
        f"{S3}/library_hits.tsv",
    threads: config["threads"]
    conda:
        "../envs/vsearch.yaml"
    log:
        f"{S3}/logs/library_search.log",
    shell:
        "vsearch --usearch_global {input.query} --db {input.db} --id 0.8 --strand both "
        "--maxaccepts 0 --maxrejects 0 --qmask none --dbmask none --threads {threads} "
        "--userout {output} --userfields query+target+id+qcov+tcov &> {log} || true; touch {output}"


rule suggestions:
    input:
        final=f"{FINAL}.tsv",
        translation=f"{S3}/translation.tsv",
        aln=f"{FINAL}_with_references.aln",
        library=LIB,
        hits=f"{S3}/library_hits.tsv",
        chimeras=f"{S3}/chimeras.fasta",
        collapsed=f"{OUT}/{NAME}_otu_table_100_with_reps.xlsx",
        rules=RULES,
        refs=REFS,
    output:
        xlsx=f"{OUT}/{NAME}_subtype_suggestions.xlsx",
        evidence=f"{S3}/evidence.tsv",
    conda:
        "../envs/r.yaml"
    log:
        f"{S3}/logs/suggestions.log",
    script:
        "../scripts/suggestions.R"


rule review_figures:
    input:
        treefile=f"{FINAL}_with_references.treefile",
        aln=f"{FINAL}_with_references.aln",
        evidence=f"{S3}/evidence.tsv",
        rules=RULES,
    output:
        tree_cladogram=f"{OUT}/review_figures/{NAME}_tree_cladogram.pdf",
        tree_distance=f"{OUT}/review_figures/{NAME}_tree_distance.pdf",
        heatmap=f"{OUT}/review_figures/{NAME}_reference_identity_heatmap.pdf",
        variable=f"{OUT}/review_figures/{NAME}_variable_region_alignment.pdf",
        alignment=f"{OUT}/review_figures/{NAME}_alignment.pdf",
    conda:
        "../envs/r.yaml"
    log:
        f"{S3}/logs/review_figures.log",
    script:
        "../scripts/review_figures.R"


rule review_page:
    """One self-contained HTML page: tree and alignment with linked selection and subtype assignment."""
    input:
        treefile=f"{FINAL}_with_references.treefile",
        aln=f"{FINAL}_with_references.aln",
        evidence=f"{S3}/evidence.tsv",
        suggestions=f"{OUT}/{NAME}_subtype_suggestions.xlsx",
        collapsed=f"{OUT}/{NAME}_otu_table_100_with_reps.xlsx",
        rules=RULES,
    output:
        html=f"{OUT}/{NAME}_subtype_review.html",
    params:
        analysis_name=NAME,
    conda:
        "../envs/r.yaml"
    log:
        f"{S3}/logs/review_page.log",
    script:
        "../scripts/review_page.R"
