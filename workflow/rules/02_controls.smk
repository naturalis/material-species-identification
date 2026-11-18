# Step 2: control characterisation and the bleed test. The read-level half lives here. The
# taxonomic half needs the step 3 classifier and is the stub control_composition.


rule control_profile:
    input:
        lambda wc: merged_fastqs(wc.control),
    output:
        profile=f"{RESULTS}/02_controls/profiles/{{control}}.profile.tsv",
        top=f"{RESULTS}/02_controls/profiles/{{control}}.top_sequences.tsv",
    params:
        id=lambda wc: SHEET.loc[wc.control, "id"],
        motifs=" ".join(config["preprocess"]["adapter_motifs"]),
        cfg=config["controls"],
    log:
        f"{LOGS}/02_controls/{{control}}.profile.log",
    resources:
        mem_mb=config["mem_mb"]["control_profile"],
    shell:
        """
        seqkit seq -s -w 0 {input} 2> {log} \
            | python {SCRIPTS}/control_profile.py --id {params.id} \
                --top {params.cfg[top_sequences]} \
                --low-complexity-fraction {params.cfg[low_complexity_fraction]} \
                --adapter-motifs {params.motifs} \
                --out-profile {output.profile} --out-top {output.top} 2>> {log}
        """


rule blank_profile:
    input:
        readstats=f"{RESULTS}/01_preprocess/readstats.tsv",
        profiles=expand(f"{RESULTS}/02_controls/profiles/{{control}}.profile.tsv", control=CONTROLS),
        top=expand(
            f"{RESULTS}/02_controls/profiles/{{control}}.top_sequences.tsv", control=CONTROLS
        ),
    output:
        profile=f"{RESULTS}/02_controls/blank_profile.tsv",
        top=f"{RESULTS}/02_controls/control_top_sequences.tsv",
    log:
        f"{LOGS}/02_controls/blank_profile.log",
    shell:
        """
        python {SCRIPTS}/blank_profile.py --readstats {input.readstats} \
            --profiles {input.profiles} --top {input.top} \
            --out-profile {output.profile} --out-top {output.top} > {log} 2>&1
        """


def control_arguments(wildcards):
    return " ".join(
        f"--control {SHEET.loc[c, 'id']}={','.join(merged_fastqs(c))}" for c in CONTROLS
    )


# One job per sample source. The sources are compared one at a time, never pooled.
rule bleed_scan:
    input:
        source=lambda wc: merged_fastqs(wc.source),
        controls=[path for control in CONTROLS for path in merged_fastqs(control)],
    output:
        table=f"{RESULTS}/02_controls/bleed/{{source}}.tsv",
        top=f"{RESULTS}/02_controls/bleed/{{source}}.top_shared.tsv",
    params:
        id=lambda wc: SHEET.loc[wc.source, "id"],
        controls=control_arguments,
        motifs=" ".join(config["preprocess"]["adapter_motifs"]),
        cfg=config["bleed"],
    log:
        f"{LOGS}/02_controls/{{source}}.bleed_scan.log",
    threads: config["threads"]["bleed_scan"]
    resources:
        mem_mb=config["mem_mb"]["bleed_scan"],
    shell:
        """
        python {SCRIPTS}/bleed_scan.py --source-id {params.id} --source {input.source} \
            {params.controls} \
            --min-length {params.cfg[min_length]} --probe-reads {params.cfg[probe_reads]} \
            --exclude-motifs {params.motifs} \
            --min-self-hits {params.cfg[min_self_hits]} \
            --min-expected-shared {params.cfg[min_expected_shared]} \
            --flag-fraction {params.cfg[flag_fraction]} \
            --top-shared {params.cfg[top_shared]} --seed {params.cfg[seed]} \
            --out {output.table} --out-top {output.top} > {log} 2>&1
        """


rule bleed_test:
    input:
        tables=expand(f"{RESULTS}/02_controls/bleed/{{source}}.tsv", source=SAMPLES),
        top=expand(f"{RESULTS}/02_controls/bleed/{{source}}.top_shared.tsv", source=SAMPLES),
    output:
        table=f"{RESULTS}/02_controls/bleed_test.tsv",
        top=f"{RESULTS}/02_controls/bleed_top_shared.tsv",
    shell:
        """
        awk 'FNR > 1 || NR == 1' {input.tables} > {output.table}
        awk 'FNR > 1 || NR == 1' {input.top} > {output.top}
        """
