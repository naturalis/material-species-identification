# Step 3: taxonomic screen, blastn plus an LCA in python. Evidence for
# Gate 1. Identical reads are searched once and carry their count in the query name.


rule screen_unique:
    input:
        lambda wc: merged_fastqs(wc.library),
    output:
        fasta=f"{RESULTS}/03_screen/{{library}}/{{library}}.unique.fasta.gz",
        stats=f"{RESULTS}/03_screen/{{library}}/{{library}}.unique.tsv",
    params:
        tmp=f"{RESULTS}/03_screen/{{library}}/tmp",
    threads: config["threads"]["screen_unique"]
    resources:
        mem_mb=config["mem_mb"]["screen_unique"],
    log:
        f"{LOGS}/03_screen/{{library}}.unique.log",
    shell:
        """
        mkdir -p {params.tmp}
        seqkit seq -s -w 0 {input} 2> {log} \\
            | LC_ALL=C sort -S {resources.mem_mb}M --parallel {threads} -T {params.tmp} \\
            | uniq -c \\
            | awk -v stats={output.stats} 'BEGIN {{ OFS = "\\t" }}
                {{ reads += $1; print ">u" NR "_" $1; print $2 }}
                END {{ print "reads", "distinct" > stats; print reads, NR > stats }}' \\
            | gzip > {output.fasta}
        rm -rf {params.tmp}
        """


rule screen_blastn:
    input:
        fasta=rules.screen_unique.output.fasta,
        db=rules.refseq_blastdb.output.marker,
    output:
        f"{RESULTS}/03_screen/{{library}}/{{library}}.blast.tsv.gz",
    params:
        db=config["references"]["refseq_mito_blastdb"],
        task=config["screen"]["blastn_task"],
        word_size=config["screen"]["blastn_word_size"],
        evalue=config["screen"]["blastn_evalue"],
        max_target_seqs=config["screen"]["blastn_max_target_seqs"],
    threads: config["threads"]["blastn"]
    log:
        f"{LOGS}/03_screen/{{library}}.blastn.log",
    conda:
        "../envs/blast.yaml"
    shell:
        """
        gzip -dc {input.fasta} \\
            | blastn -db {params.db} -task {params.task} -word_size {params.word_size} \\
                -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} -max_hsps 1 \\
                -num_threads {threads} \\
                -outfmt "6 qseqid saccver pident length bitscore qlen" 2> {log} \\
            | gzip > {output}
        """


rule screen_lca:
    input:
        rules.screen_blastn.output,
    output:
        f"{RESULTS}/03_screen/{{library}}/{{library}}.lca.tsv",
    params:
        taxonomy=config["references"]["refseq_mito_taxonomy"],
        identity=config["screen"]["strict_identity"],
        coverage=config["screen"]["strict_coverage"],
        artefacts=" ".join(config["screen"]["artefact_records"]),
    log:
        f"{LOGS}/03_screen/{{library}}.lca.log",
    shell:
        """
        python {SCRIPTS}/screen_lca.py --hits {input} --taxonomy {params.taxonomy} \\
            --strict-identity {params.identity} --strict-coverage {params.coverage} \\
            --artefact-records {params.artefacts} --out {output} > {log} 2>&1
        """


rule taxonomic_profile:
    input:
        lca=expand(f"{RESULTS}/03_screen/{{library}}/{{library}}.lca.tsv", library=LIBRARIES),
        stats=expand(f"{RESULTS}/03_screen/{{library}}/{{library}}.unique.tsv", library=LIBRARIES),
    output:
        f"{RESULTS}/03_screen/taxonomic_profile.tsv",
    params:
        sheet=config["sample_sheet"],
        results=f"{RESULTS}/03_screen",
    log:
        f"{LOGS}/03_screen/taxonomic_profile.log",
    shell:
        """
        python {SCRIPTS}/taxonomic_profile.py --sheet {params.sheet} \\
            --results {params.results} --out {output} > {log} 2>&1
        """


rule krona:
    input:
        f"{RESULTS}/03_screen/taxonomic_profile.tsv",
    output:
        f"{RESULTS}/03_screen/krona.html",
    shell:
        STUB


# The taxonomic half of step 2: blank composition, the revised Gate 1 threshold, and the
# composition comparison of the low-yield controls against the sample libraries.
rule control_composition:
    input:
        f"{RESULTS}/03_screen/taxonomic_profile.tsv",
    output:
        f"{RESULTS}/02_controls/control_composition.tsv",
    shell:
        STUB


rule gate1_evidence:
    input:
        f"{RESULTS}/03_screen/taxonomic_profile.tsv",
    output:
        f"{RESULTS}/03_screen/gate1_evidence.tsv",
    params:
        blank=config["gates"]["gate1_blank"],
        min_reads=config["gates"]["gate1_min_reads"],
        fold=config["gates"]["gate1_fold_over_blank"],
        zero_limit=config["gates"]["gate1_zero_count_limit"],
    log:
        f"{LOGS}/03_screen/gate1_evidence.log",
    shell:
        """
        python {SCRIPTS}/gate1_evidence.py --profile {input} --blank {params.blank} \\
            --min-reads {params.min_reads} --fold {params.fold} \\
            --zero-limit {params.zero_limit} --out {output} > {log} 2>&1
        """
