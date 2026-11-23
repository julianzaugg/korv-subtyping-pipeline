# Stage 1 - preprocessing, per sample (replaces the CLC steps).
# merge overlapping pairs -> trim ENV primers (both required) -> quality/length/N filter
# -> subsample -> relabel "<internal_id>-<read number>".
# Only used when start_from = fastq.

PRE = f"{WORK}/samples"
SUBSAMPLE = config["subsample_depth"]


rule merge_pairs:
    input:
        unpack(raw_reads),
    output:
        temp(f"{PRE}/{{sid}}/merged.fastq"),
    params:
        minovlen=config["merge_min_overlap"],
        maxdiffs=config.get("merge_max_diffs", 10),
    threads: 1
    conda:
        "../envs/vsearch.yaml"
    log:
        f"{WORK}/logs/{{sid}}.merge.log",
    shell:
        "vsearch --fastq_mergepairs {input.r1} --reverse {input.r2} --fastqout {output} "
        "--fastq_minovlen {params.minovlen} --fastq_maxdiffs {params.maxdiffs} "
        "--threads {threads} &> {log}"


rule trim_primers:
    """Both primers must be present (CLC 'Discard when not found'); then length, N and
    quality filtering. Primers are removed first so end-quality trimming cannot cut
    into them and cause a read to be discarded."""
    input:
        f"{PRE}/{{sid}}/merged.fastq",
    output:
        primer_trimmed=temp(f"{PRE}/{{sid}}/primer_trimmed.fastq"),
        trimmed=temp(f"{PRE}/{{sid}}/trimmed.fastq"),
    params:
        fwd=config["primer_fwd"],
        rev_rc=config["primer_rev_rc"],
        err=config["primer_error_rate"],
        qual=config["quality_cutoff"],
        minlen=config["min_len"],
        maxlen=config["max_len"],
        maxn=config["max_n"],
    threads: 1
    conda:
        "../envs/cutadapt.yaml"
    log:
        f"{WORK}/logs/{{sid}}.cutadapt.log",
    shell:
        "(cutadapt -g {params.fwd}...{params.rev_rc} -e {params.err} --discard-untrimmed "
        "-j {threads} -o {output.primer_trimmed} {input} && "
        "cutadapt -q {params.qual},{params.qual} -m {params.minlen} -M {params.maxlen} "
        "--max-n {params.maxn} -j {threads} -o {output.trimmed} {output.primer_trimmed}) &> {log}"


rule subsample_relabel:
    input:
        f"{PRE}/{{sid}}/trimmed.fastq",
    output:
        f"{PRE}/{{sid}}/reads.fasta",
    params:
        depth="" if SUBSAMPLE is None else int(SUBSAMPLE),
        seed=config["subsample_seed"],
    conda:
        "../envs/cutadapt.yaml"
    log:
        f"{WORK}/logs/{{sid}}.subsample.log",
    shell:
        r"""
        exec 2> {log}
        total=$(seqkit stats -T {input} | awk 'NR==2{{print $4}}')
        if [ -z "{params.depth}" ] || [ "$total" -le "{params.depth}" ]; then
            seqkit fq2fa -w 0 {input}
        else
            seqkit sample -n {params.depth} -s {params.seed} -2 {input} | seqkit fq2fa -w 0
        fi | awk -v s={wildcards.sid} 'NR%2==1{{printf(">%s-%07d\n", s, ++n)}} NR%2==0{{print toupper($0)}}' > {output}
        """


rule read_stats:
    """One row per sample, same columns as the legacy *_data.xlsx read statistics."""
    input:
        raw=lambda wc: BY_INTERNAL[wc.sid]["R1"],
        merged=f"{PRE}/{{sid}}/merged.fastq",
        trimmed=f"{PRE}/{{sid}}/trimmed.fastq",
        reads=f"{PRE}/{{sid}}/reads.fasta",
    output:
        f"{PRE}/{{sid}}/read_stats.tsv",
    params:
        sample_id=lambda wc: BY_INTERNAL[wc.sid]["sample_id"],
    conda:
        "../envs/cutadapt.yaml"
    log:
        f"{WORK}/logs/{{sid}}.read_stats.log",
    shell:
        r"""
        exec 2> {log}
        st() {{ seqkit stats -T "$1" | awk 'NR==2{{print $4"\t"$7}}'; }}
        read raw_n raw_len < <(st {input.raw})
        read mrg_n mrg_len < <(st {input.merged})
        read trm_n trm_len < <(st {input.trimmed})
        read sub_n sub_len < <(st {input.reads})
        awk -v id="{params.sample_id}" -v raw=$raw_n -v mn=$mrg_n -v ml=$mrg_len -v tn=$trm_n -v tl=$trm_len -v sn=$sub_n \
          'BEGIN{{OFS="\t"; pct=(mn>0)?100*tn/mn:0;
                  print "sample_id","Raw read pairs","Number of reads","Avg.length","Number of reads after trim","Percentage trimmed (%)","Avg.length after trim","Number of reads subsampled";
                  print id,raw,mn,ml,tn,sprintf("%.2f",pct),tl,sn}}' > {output}
        """


rule pool_reads:
    input:
        expand(f"{PRE}/{{sid}}/reads.fasta", sid=INTERNAL),
    output:
        f"{WORK}/pooled_reads.fasta",
    conda:
        "../envs/cutadapt.yaml"
    log:
        f"{WORK}/logs/pool_reads.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule sample_read_stats:
    input:
        expand(f"{PRE}/{{sid}}/read_stats.tsv", sid=INTERNAL),
    output:
        f"{OUT}/{NAME}_sample_read_stats.tsv",
    conda:
        "../envs/cutadapt.yaml"
    log:
        f"{WORK}/logs/sample_read_stats.log",
    shell:
        "awk 'FNR==1 && NR!=1{{next}} {{print}}' {input} > {output} 2> {log}"
