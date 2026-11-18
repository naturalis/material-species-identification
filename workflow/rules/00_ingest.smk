# Step 0: ingest. The sheet itself comes from workflow/scripts/make_sample_sheet.py.


rule ingest_check:
    input:
        sheet=config["sample_sheet"],
    output:
        report=f"{RESULTS}/00_ingest/ingest_report.tsv",
    log:
        f"{LOGS}/00_ingest/ingest_check.log",
    script:
        "../scripts/ingest_check.py"


# Off by default: see the note on ingest.md5 in config.yaml.
rule md5:
    input:
        r1=lambda wc: SHEET.loc[wc.library, "R1"],
        r2=lambda wc: SHEET.loc[wc.library, "R2"],
    output:
        f"{RESULTS}/00_ingest/md5/{{library}}.md5",
    log:
        f"{LOGS}/00_ingest/md5/{{library}}.log",
    shell:
        "md5sum {input.r1} {input.r2} > {output} 2> {log}"


rule md5_manifest:
    input:
        expand(f"{RESULTS}/00_ingest/md5/{{library}}.md5", library=LIBRARIES),
    output:
        f"{RESULTS}/00_ingest/md5_manifest.tsv",
    shell:
        "cat {input} > {output}"
