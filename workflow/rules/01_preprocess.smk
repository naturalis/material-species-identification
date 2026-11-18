# Step 1: adapter trimming, pair merging, read survival and fragment lengths. All libraries.


rule adapterremoval:
    input:
        r1=lambda wc: SHEET.loc[wc.library, "R1"],
        r2=lambda wc: SHEET.loc[wc.library, "R2"],
        ingest=f"{RESULTS}/00_ingest/ingest_report.tsv",
    output:
        collapsed=f"{RESULTS}/01_preprocess/{{library}}/{{library}}.collapsed.gz",
        truncated=f"{RESULTS}/01_preprocess/{{library}}/{{library}}.collapsed.truncated.gz",
        pair1=f"{RESULTS}/01_preprocess/{{library}}/{{library}}.pair1.truncated.gz",
        pair2=f"{RESULTS}/01_preprocess/{{library}}/{{library}}.pair2.truncated.gz",
        singleton=f"{RESULTS}/01_preprocess/{{library}}/{{library}}.singleton.truncated.gz",
        discarded=f"{RESULTS}/01_preprocess/{{library}}/{{library}}.discarded.gz",
        settings=f"{RESULTS}/01_preprocess/{{library}}/{{library}}.settings",
    params:
        basename=f"{RESULTS}/01_preprocess/{{library}}/{{library}}",
        cfg=config["preprocess"],
    log:
        f"{LOGS}/01_preprocess/{{library}}.adapterremoval.log",
    threads: config["threads"]["adapterremoval"]
    resources:
        mem_mb=config["mem_mb"]["default"],
    shell:
        """
        AdapterRemoval \
            --file1 {input.r1} --file2 {input.r2} \
            --basename {params.basename} --gzip --threads {threads} \
            --adapter1 {params.cfg[adapter1]} --adapter2 {params.cfg[adapter2]} \
            --collapse --collapse-deterministic \
            --minlength {params.cfg[min_length]} \
            --minalignmentlength {params.cfg[min_overlap]} \
            --mm {params.cfg[mismatch_rate]} \
            --qualitymax {params.cfg[quality_max]} \
            --seed {params.cfg[seed]} \
            > {log} 2>&1
        """


rule fragment_lengths:
    input:
        lambda wc: merged_fastqs(wc.library),
    output:
        f"{RESULTS}/01_preprocess/{{library}}/{{library}}.lengths.tsv",
    log:
        f"{LOGS}/01_preprocess/{{library}}.lengths.log",
    threads: config["threads"]["seqkit"]
    resources:
        mem_mb=config["mem_mb"]["default"],
    shell:
        """
        (printf 'length\\treads\\n'
         seqkit fx2tab -n -l -j {threads} {input} 2> {log} \
             | awk -F '\\t' '{{ n[$2]++ }} END {{ for (l in n) print l "\\t" n[l] }}' \
             | sort -n) > {output}
        """


rule readstats:
    input:
        sheet=config["sample_sheet"],
        ingest=f"{RESULTS}/00_ingest/ingest_report.tsv",
        settings=expand(
            f"{RESULTS}/01_preprocess/{{library}}/{{library}}.settings", library=LIBRARIES
        ),
        lengths=expand(
            f"{RESULTS}/01_preprocess/{{library}}/{{library}}.lengths.tsv", library=LIBRARIES
        ),
    output:
        stats=f"{RESULTS}/01_preprocess/readstats.tsv",
        lengths=f"{RESULTS}/01_preprocess/length_distributions.tsv",
    params:
        directory=f"{RESULTS}/01_preprocess",
        cfg=config["preprocess"],
    log:
        f"{LOGS}/01_preprocess/readstats.log",
    shell:
        """
        python {SCRIPTS}/readstats.py \
            --sheet {input.sheet} --ingest {input.ingest} --directory {params.directory} \
            --merge-rate-min {params.cfg[merge_rate_min]} \
            --merge-rate-max {params.cfg[merge_rate_max]} \
            --out-stats {output.stats} --out-lengths {output.lengths} > {log} 2>&1
        """


rule length_plot:
    input:
        f"{RESULTS}/01_preprocess/length_distributions.tsv",
    output:
        f"{RESULTS}/01_preprocess/length_distributions.png",
    params:
        cfg=config["preprocess"],
    log:
        f"{LOGS}/01_preprocess/length_plot.log",
    shell:
        """
        python {SCRIPTS}/plot_lengths.py --lengths {input} \
            --read-length {params.cfg[read_length]} --min-overlap {params.cfg[min_overlap]} \
            --out {output} > {log} 2>&1
        """
