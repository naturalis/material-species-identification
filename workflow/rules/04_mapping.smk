# Step 4: competitive mapping against the panel. Evidence for Gate 2, which stays manual.
#
# Input is step 3's exact-duplicate-collapsed FASTA: identical sequences map identically,
# and the count stays in the read name. The combined-panel mapping is only a prefilter, since
# a read that maps to any one panel genome also maps somewhere in the combined panel. The
# mapped reads are then mapped to every panel genome on its own, and assign_panel compares
# edit distances across genomes (§6: per-reference MAPQ, margin of 2 mismatches).

INDEX = f"{RESULTS}/04_mapping/index"


rule panel_index:
    input:
        config["references"]["panel_fasta"],
    output:
        marker=f"{INDEX}/index.done",
    params:
        dir=INDEX,
    log:
        f"{LOGS}/04_mapping/panel_index.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        rm -rf {params.dir}/panel.fasta* {params.dir}/single
        mkdir -p {params.dir}/single
        cp {input} {params.dir}/panel.fasta
        bwa index {params.dir}/panel.fasta > {log} 2>&1
        samtools faidx {params.dir}/panel.fasta
        for accession in $(cut -f1 {params.dir}/panel.fasta.fai); do
            samtools faidx {params.dir}/panel.fasta "$accession" > {params.dir}/single/"$accession".fasta
            bwa index {params.dir}/single/"$accession".fasta >> {log} 2>&1
        done
        touch {output.marker}
        """


rule map_panel:
    input:
        fasta=rules.screen_unique.output.fasta,
        index=rules.panel_index.output.marker,
        gate1=f"{RESULTS}/03_screen/gate1_evidence.tsv",
    output:
        bam=f"{RESULTS}/04_mapping/{{library}}/{{library}}.panel.bam",
        mapped=f"{RESULTS}/04_mapping/{{library}}/{{library}}.mapped.fasta.gz",
    params:
        aln=config["mapping"]["bwa_aln"],
        panel=f"{INDEX}/panel.fasta",
        sai=f"{RESULTS}/04_mapping/{{library}}/{{library}}.panel.sai",
    threads: config["threads"]["bwa"]
    log:
        f"{LOGS}/04_mapping/{{library}}.map_panel.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bwa aln -t {threads} {params.aln} {params.panel} {input.fasta} > {params.sai} 2> {log}
        bwa samse {params.panel} {params.sai} {input.fasta} 2>> {log} \\
            | samtools view -u -F 4 - \\
            | samtools sort -@ 2 -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        samtools fasta {output.bam} 2>> {log} | gzip > {output.mapped}
        rm {params.sai}
        """


rule assign_panel:
    input:
        mapped=rules.map_panel.output.mapped,
        index=rules.panel_index.output.marker,
    output:
        assigned=f"{RESULTS}/04_mapping/{{library}}/{{library}}.assigned.tsv.gz",
        metrics=f"{RESULTS}/04_mapping/{{library}}/{{library}}.metrics.tsv",
        ambiguous=f"{RESULTS}/04_mapping/{{library}}/{{library}}.ambiguous.tsv",
    params:
        aln=config["mapping"]["bwa_aln"],
        panel=f"{INDEX}/panel.fasta",
        single=f"{INDEX}/single",
        tmp=f"{RESULTS}/04_mapping/{{library}}/single",
        margin=config["mapping"]["panel_margin_mismatches"],
        min_mapq=config["mapping"]["min_mapq"],
    threads: config["threads"]["bwa"]
    log:
        f"{LOGS}/04_mapping/{{library}}.assign_panel.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        rm -rf {params.tmp}
        mkdir -p {params.tmp}
        for reference in {params.single}/*.fasta; do
            accession=$(basename "$reference" .fasta)
            bwa aln -t {threads} {params.aln} "$reference" {input.mapped} > {params.tmp}/"$accession".sai 2>> {log}
            bwa samse "$reference" {params.tmp}/"$accession".sai {input.mapped} 2>> {log} \\
                | samtools view -b -F 4 -o {params.tmp}/"$accession".bam -
            rm {params.tmp}/"$accession".sai
        done
        python {SCRIPTS}/assign_panel.py --bams {params.tmp}/*.bam --panel {params.panel} \\
            --library {wildcards.library} \\
            --margin {params.margin} --min-mapq {params.min_mapq} \\
            --assigned {output.assigned} --metrics {output.metrics} \\
            --ambiguous {output.ambiguous} >> {log} 2>&1
        rm -rf {params.tmp}
        """


rule per_reference_metrics:
    input:
        metrics=expand(rules.assign_panel.output.metrics, library=LIBRARIES),
        ambiguous=expand(rules.assign_panel.output.ambiguous, library=LIBRARIES),
        assigned=expand(rules.assign_panel.output.assigned, library=LIBRARIES),
    output:
        metrics=f"{RESULTS}/04_mapping/per_reference_metrics.tsv",
        ambiguous=f"{RESULTS}/04_mapping/ambiguity_bin.tsv",
        tracks=f"{RESULTS}/04_mapping/coverage_tracks.png",
    params:
        sheet=config["sample_sheet"],
        blank=config["gates"]["gate1_blank"],
        min_breadth=config["gates"]["gate2_min_breadth"],
        tracks=config["mapping"]["coverage_track_references"],
    log:
        f"{LOGS}/04_mapping/per_reference_metrics.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/panel_metrics.py --sheet {params.sheet} --metrics {input.metrics} \\
            --ambiguous {input.ambiguous} --min-breadth {params.min_breadth} \\
            --out {output.metrics} --ambiguous-out {output.ambiguous} > {log} 2>&1
        python {SCRIPTS}/plot_coverage_tracks.py --sheet {params.sheet} --blank {params.blank} \\
            --assigned {input.assigned} --metrics {output.metrics} --top {params.tracks} \\
            --out {output.tracks} >> {log} 2>&1
        """
