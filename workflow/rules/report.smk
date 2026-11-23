# Stage 4 - final report. Only runs once the expert-completed assignments file is configured.

rule sample_map:
    """sample sheet with internal ids, so reports show the original sample_id plus metadata"""
    output:
        f"{WORK}/sample_map.tsv",
    params:
        rows=SAMPLE_ROWS,
    conda:
        "../envs/cutadapt.yaml"
    log:
        f"{WORK}/logs/sample_map.log",
    script:
        "../scripts/sample_map.py"


rule final_report:
    input:
        unpack(report_inputs),
    output:
        f"{OUT}/{NAME}_final_result.xlsx",
    params:
        extra_categories=config.get("extra_categories", ["Endogenous"]),
    conda:
        "../envs/r.yaml"
    log:
        f"{WORK}/logs/final_report.log",
    script:
        "../scripts/report.R"
