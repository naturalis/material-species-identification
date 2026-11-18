# Step 5: what is left of authentication after USER. Every number is reported next to the blank.
#
# Whether USER was full or partial is not known (§10.7). Full USER removes deamination
# almost everywhere, so a flat damage profile is uninformative. Partial USER leaves it on the
# terminal base, so a raised first position in one taxon and not another would separate old
# from fresh DNA. The profile is computed either way.


rule damage_profile:
    input:
        bam=rules.map_panel.output.bam,
        assigned=rules.assign_panel.output.assigned,
        index=rules.panel_index.output.marker,
    output:
        damage=f"{RESULTS}/05_auth/{{library}}/{{library}}.damage.tsv",
        lengths=f"{RESULTS}/05_auth/{{library}}/{{library}}.lengths.tsv",
    params:
        panel=f"{INDEX}/panel.fasta",
        positions=config["authentication"]["damage_positions"],
    log:
        f"{LOGS}/05_auth/{{library}}.damage_profile.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/damage_profile.py --bam {input.bam} --assigned {input.assigned} \\
            --panel {params.panel} --library {wildcards.library} --positions {params.positions} \\
            --damage {output.damage} --lengths {output.lengths} > {log} 2>&1
        """


rule complexity:
    input:
        rules.screen_unique.output.fasta,
    output:
        histogram=f"{RESULTS}/05_auth/{{library}}/{{library}}.duplicate_histogram.tsv",
        curve=f"{RESULTS}/05_auth/{{library}}/{{library}}.lc_extrap.tsv",
    params:
        step=config["authentication"]["preseq_step"],
    log:
        f"{LOGS}/05_auth/{{library}}.complexity.log",
    conda:
        "../envs/preseq.yaml"
    shell:
        """
        # Copies per distinct sequence, from the u<n>_<count> names, as preseq's -H histogram.
        gzip -dc {input} | awk '/^>/ {{ sub(/^>u[0-9]+_/, ""); n[$0]++ }}
            END {{ for (c in n) print c "\\t" n[c] }}' | sort -k1,1n > {output.histogram}
        # A saturated or tiny library can defeat the extrapolation; keep the header and say so.
        preseq lc_extrap -H {output.histogram} -s {params.step} -o {output.curve} > {log} 2>&1 \\
            || printf 'TOTAL_READS\\tEXPECTED_DISTINCT\\tLOWER_0.95CI\\tUPPER_0.95CI\\n' > {output.curve}
        """


rule auth_summary:
    input:
        damage=expand(rules.damage_profile.output.damage, library=LIBRARIES),
        lengths=expand(rules.damage_profile.output.lengths, library=LIBRARIES),
        curves=expand(rules.complexity.output.curve, library=LIBRARIES),
        histograms=expand(rules.complexity.output.histogram, library=LIBRARIES),
        library_lengths=f"{RESULTS}/01_preprocess/length_distributions.tsv",
    output:
        damage=f"{RESULTS}/05_auth/damage_summary.tsv",
        damage_plot=f"{RESULTS}/05_auth/damage_by_taxon.png",
        lengths=f"{RESULTS}/05_auth/length_by_taxon.tsv",
        lengths_plot=f"{RESULTS}/05_auth/length_by_taxon.png",
        complexity=f"{RESULTS}/05_auth/complexity.tsv",
    params:
        sheet=config["sample_sheet"],
        blank=config["gates"]["gate1_blank"],
        taxa=config["authentication"]["taxa"],
        interior=config["authentication"]["interior_positions"],
    log:
        f"{LOGS}/05_auth/auth_summary.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/auth_summary.py --sheet {params.sheet} --blank {params.blank} \\
            --taxa {params.taxa:q} --interior {params.interior} \\
            --damage {input.damage} --lengths {input.lengths} \\
            --curves {input.curves} --histograms {input.histograms} \\
            --library-lengths {input.library_lengths} \\
            --damage-out {output.damage} --damage-plot {output.damage_plot} \\
            --lengths-out {output.lengths} --lengths-plot {output.lengths_plot} \\
            --complexity-out {output.complexity} > {log} 2>&1
        """


# Added because step 5's control point is "every number
# reported alongside its blank equivalent", and a bare count next to a zero is not a
# comparison until the depth behind the zero is stated. See section 10.10.
rule background_test:
    input:
        metrics=f"{RESULTS}/04_mapping/per_reference_metrics.tsv",
        complexity=rules.auth_summary.output.complexity,
    output:
        table=f"{RESULTS}/05_auth/background_test.tsv",
        figure=f"{RESULTS}/05_auth/background_test.png",
    params:
        blank=config["gates"]["gate1_blank"],
    log:
        f"{LOGS}/05_auth/background_test.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/background_test.py --metrics {input.metrics} \\
            --complexity {input.complexity} --blank {params.blank} \\
            --out {output.table} --figure {output.figure} > {log} 2>&1
        """


# Needs a consensus of the endogenous human mitogenome, so it follows step 6.
rule contammix:
    input:
        bams=expand(f"{RESULTS}/04_mapping/{{library}}/{{library}}.panel.bam", library=LIBRARIES),
        metrics=f"{RESULTS}/04_mapping/per_reference_metrics.tsv",
    output:
        f"{RESULTS}/05_auth/contammix.tsv",
    shell:
        STUB
