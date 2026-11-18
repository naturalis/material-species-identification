# Reference acquisition. `snakemake make_references --sdm conda`. Not part of `all`.


rule refseq_mito_fetch:
    output:
        fasta=config["references"]["refseq_mito_fasta"],
        taxonomy=config["references"]["refseq_mito_taxonomy"],
        release=f"{REFSEQ_DIR}/RELEASE_NUMBER",
        manifest=f"{REFSEQ_DIR}/manifest.part.tsv",
    params:
        cfg=config["references"],
        outdir=REFSEQ_DIR,
    log:
        f"{LOGS}/references/refseq_mito_fetch.log",
    shell:
        """
        PYTHONPATH={SCRIPTS} python {SCRIPTS}/fetch_refseq_mito.py \
            --url {params.cfg[refseq_mito_url]} --release-url {params.cfg[refseq_release_url]} \
            --outdir {params.outdir} --fasta {output.fasta} --taxonomy {output.taxonomy} \
            --release-file {output.release} --manifest {output.manifest} > {log} 2>&1
        """


rule refseq_blastdb:
    input:
        config["references"]["refseq_mito_fasta"],
    output:
        marker=f"{config['references']['refseq_mito_blastdb']}.njs",
        version=f"{Path(config['references']['refseq_mito_blastdb']).parent.parent}/makeblastdb.version",
    params:
        db=config["references"]["refseq_mito_blastdb"],
    log:
        f"{LOGS}/references/refseq_blastdb.log",
    conda:
        "../envs/blast.yaml"
    shell:
        """
        makeblastdb -version > {output.version}
        makeblastdb -in {input} -dbtype nucl -parse_seqids -out {params.db} > {log} 2>&1
        """


rule panel_fetch:
    input:
        request=config["references"]["panel_csv"],
        taxonomy=config["references"]["refseq_mito_taxonomy"],
    output:
        fasta=config["references"]["panel_fasta"],
        table=config["references"]["panel_tsv"],
        manifest=f"{PANEL_DIR}/manifest.part.tsv",
    params:
        cfg=config["references"],
        outdir=PANEL_DIR,
    log:
        f"{LOGS}/references/panel_fetch.log",
    shell:
        """
        PYTHONPATH={SCRIPTS} python {SCRIPTS}/fetch_panel.py --request {input.request} \
            --efetch-url {params.cfg[efetch_url]} --delay {params.cfg[efetch_delay]} \
            --outdir {params.outdir} --refseq-taxonomy {input.taxonomy} \
            --panel-fasta {output.fasta} --panel-tsv {output.table} \
            --manifest {output.manifest} > {log} 2>&1
        """


rule make_references:
    input:
        refseq=f"{REFSEQ_DIR}/manifest.part.tsv",
        panel=f"{PANEL_DIR}/manifest.part.tsv",
        cds=f"{PANEL_DIR}/cds/manifest.part.tsv",
        lineage=expand("resources/lineage/{lineage}/manifest.part.tsv",
                       lineage=[s.replace(" ", "_") for s in config["barcode"].get("lineage", {})]),
        marker=f"{config['references']['refseq_mito_blastdb']}.njs",
        version=f"{Path(config['references']['refseq_mito_blastdb']).parent.parent}/makeblastdb.version",
    output:
        "resources/MANIFEST.tsv",
    params:
        cfg=config["references"],
    log:
        f"{LOGS}/references/merge_manifest.log",
    shell:
        """
        PYTHONPATH={SCRIPTS} python {SCRIPTS}/merge_manifest.py \
            --parts {input.refseq} {input.panel} {input.cds} {input.lineage} --blastdb {params.cfg[refseq_mito_blastdb]} \
            --blast-version {input.version} --records {params.cfg[refseq_mito_fasta]} \
            --out {output} > {log} 2>&1
        """
