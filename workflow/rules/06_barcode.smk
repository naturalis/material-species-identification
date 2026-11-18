# Step 6: consensus sequences of the target taxa, and what they match.
#
# Targets are panel species named in config barcode.taxa (§10.11). Their reads are the ones
# step 4 assigned to that genome alone, after position deduplication.
#   Route 1: MitoGeneExtractor builds the barcode gene from the reads, guided by the protein.
#   Route 2: de novo assembly (GetOrganelle): not run, see §10.11.
#   Route 3: mapping consensus on the panel genome, masked below consensus_min_depth (§6).
# Both barcode sequences are searched against the curated and the uncurated BOLD libraries
# (§3b), locally. The human consensus gets a haplogroup and a mixture check. A taxon with a
# lineage set in config is placed against complete mitogenomes of its species and family,
# and simulated reads from its closest lineage measure how mapping recovery depends on length.

BARCODE = config["barcode"]
BARCODE_DIR = f"{RESULTS}/06_barcode"
LINEAGE_DIR = "resources/lineage"
CDS_DIR = f"{PANEL_DIR}/cds"


def slug(name):
    return name.replace(" ", "_")


TARGETS = {slug(name): name for name in BARCODE["taxa"]}
LINEAGE_TARGETS = {slug(name): name for name in BARCODE.get("lineage", {})}
HUMAN = slug(BARCODE["haplogroup_taxon"])
SAMPLE_IDS = SHEET.loc[SAMPLES, "id"].tolist()

with open(config["references"]["panel_csv"], newline="") as _fh:
    PANEL_ACCESSION = {
        row["species"]: row["accession"]
        for row in __import__("csv").DictReader(_fh)
        if row["use"] == "panel"
    }


def accession(wildcards):
    return PANEL_ACCESSION[TARGETS[wildcards.taxon]]


def single_reference(wildcards):
    return f"{INDEX}/single/{accession(wildcards)}.fasta"


wildcard_constraints:
    taxon="|".join(TARGETS),
    lineage="|".join(LINEAGE_TARGETS) or "none",


rule cds_fetch:
    output:
        table=f"{CDS_DIR}/barcode_gene.tsv",
        proteins=[f"{CDS_DIR}/{PANEL_ACCESSION[t]}.{BARCODE['gene'][0]}.faa" for t in TARGETS.values()],
        manifest=f"{CDS_DIR}/manifest.part.tsv",
    params:
        accessions=[PANEL_ACCESSION[t] for t in TARGETS.values()],
        genes=BARCODE["gene"],
        cfg=config["references"],
        outdir=CDS_DIR,
    log:
        f"{LOGS}/references/cds_fetch.log",
    shell:
        """
        PYTHONPATH={SCRIPTS} python {SCRIPTS}/fetch_cds.py --accessions {params.accessions} \\
            --gene-names {params.genes} --efetch-url {params.cfg[efetch_url]} \\
            --delay {params.cfg[efetch_delay]} --outdir {params.outdir} --table {output.table} \\
            --manifest {output.manifest} > {log} 2>&1
        """


rule lineage_fetch:
    input:
        fasta=config["references"]["refseq_mito_fasta"],
        taxonomy=config["references"]["refseq_mito_taxonomy"],
    output:
        fasta=f"{LINEAGE_DIR}/{{lineage}}/genomes.fasta",
        table=f"{LINEAGE_DIR}/{{lineage}}/genomes.tsv",
        genbank=f"{LINEAGE_DIR}/{{lineage}}/genomes.gb.gz",
        manifest=f"{LINEAGE_DIR}/{{lineage}}/manifest.part.tsv",
    params:
        query=lambda wc: BARCODE["lineage"][LINEAGE_TARGETS[wc.lineage]]["query"],
        family=lambda wc: BARCODE["lineage"][LINEAGE_TARGETS[wc.lineage]]["family"],
        delay=config["references"]["efetch_delay"],
    log:
        f"{LOGS}/references/lineage_fetch.{{lineage}}.log",
    shell:
        """
        PYTHONPATH={SCRIPTS} python {SCRIPTS}/fetch_lineage.py --query {params.query:q} \\
            --family {params.family} --refseq-fasta {input.fasta} \\
            --refseq-taxonomy {input.taxonomy} --delay {params.delay} \\
            --fasta {output.fasta} --table {output.table} --genbank {output.genbank} \\
            --manifest {output.manifest} > {log} 2>&1
        """


rule partition_taxa:
    input:
        assigned=rules.assign_panel.output.assigned,
        mapped=rules.map_panel.output.mapped,
    output:
        [f"{BARCODE_DIR}/reads/{{library}}.{taxon}.fasta.gz" for taxon in TARGETS],
    params:
        taxa=list(TARGETS.values()),
    log:
        f"{LOGS}/06_barcode/{{library}}.partition_taxa.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/partition_taxa.py --assigned {input.assigned} --mapped {input.mapped} \\
            --taxa {params.taxa:q} --out {output} > {log} 2>&1
        """


rule taxon_map:
    input:
        reads=f"{BARCODE_DIR}/reads/{{library}}.{{taxon}}.fasta.gz",
        index=rules.panel_index.output.marker,
    output:
        bam=f"{BARCODE_DIR}/{{taxon}}/bam/{{library}}.bam",
    params:
        reference=single_reference,
        aln=config["mapping"]["bwa_aln"],
        min_mapq=config["mapping"]["min_mapq"],
        id=lambda wc: SHEET.loc[wc.library, "id"],
    log:
        f"{LOGS}/06_barcode/{{taxon}}.{{library}}.taxon_map.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bwa aln {params.aln} {params.reference} {input.reads} > {output.bam}.sai 2> {log}
        bwa samse -r '@RG\\tID:{params.id}\\tSM:{params.id}' {params.reference} {output.bam}.sai \\
            {input.reads} 2>> {log} \\
            | samtools view -u -F 4 -q {params.min_mapq} - \\
            | samtools sort -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        rm {output.bam}.sai
        """


rule consensus:
    input:
        bams=expand(rules.taxon_map.output.bam, library=LIBRARIES, allow_missing=True),
        index=rules.panel_index.output.marker,
    output:
        pooled_bam=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.samples.bam",
        consensus=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.consensus.fasta",
        refcoords=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.consensus_refcoords.fasta",
        low_depth=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.consensus_refcoords.min1.fasta",
        libraries=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.libraries.fasta",
        libraries_refcoords=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.libraries_refcoords.fasta",
    params:
        samples=expand(rules.taxon_map.output.bam, library=SAMPLES, allow_missing=True),
        ids=[SHEET.loc[library, "id"] for library in LIBRARIES],
        depth=BARCODE["consensus_min_depth"],
        fraction=BARCODE["consensus_call_fraction"],
    log:
        f"{LOGS}/06_barcode/{{taxon}}.consensus.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        # simple mode: a base needs depth {params.depth} and {params.fraction} of reads, else N.
        # Insertions shown in the consensus proper; refcoords keeps reference positions.
        samtools merge -f -o {output.pooled_bam} {params.samples} 2> {log}
        samtools index {output.pooled_bam}
        call() {{
            samtools consensus -a -m simple -d "$1" -c {params.fraction} -H 1.01 \\
                --show-ins "$2" --show-del "$3" "$4" 2>> {log} | sed "1s/^>.*/>$5/"
        }}
        call {params.depth} yes no {output.pooled_bam} '{wildcards.taxon} samples pooled' > {output.consensus}
        call {params.depth} no yes {output.pooled_bam} '{wildcards.taxon} samples pooled' > {output.refcoords}
        call 1 no yes {output.pooled_bam} '{wildcards.taxon} samples pooled, depth 1' > {output.low_depth}
        : > {output.libraries}
        : > {output.libraries_refcoords}
        set -- {params.ids}
        for bam in {input.bams}; do
            call {params.depth} yes no "$bam" "$1" >> {output.libraries}
            call {params.depth} no yes "$bam" "$1" >> {output.libraries_refcoords}
            shift
        done
        """


rule consensus_stats:
    input:
        bams=expand(rules.taxon_map.output.bam, library=LIBRARIES, allow_missing=True),
        refcoords=rules.consensus.output.refcoords,
        libraries=rules.consensus.output.libraries_refcoords,
    output:
        sites=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.site_counts.tsv.gz",
        summary=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.consensus_summary.tsv",
        mixture=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.mixture.tsv",
    params:
        reference=single_reference,
        ids=[SHEET.loc[library, "id"] for library in LIBRARIES],
        depth=BARCODE["consensus_min_depth"],
        mixture_depth=BARCODE["mixture_depth"],
        minor=BARCODE["minor_fraction"],
    log:
        f"{LOGS}/06_barcode/{{taxon}}.consensus_stats.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/consensus_stats.py --reference {params.reference} --bams {input.bams} \\
            --libraries {params.ids} --samples {SAMPLE_IDS} --pooled {input.refcoords} \\
            --per-library {input.libraries} --min-depth {params.depth} \\
            --mixture-depth {params.mixture_depth} --minor-fraction {params.minor} \\
            --site-counts {output.sites} --summary {output.summary} \\
            --mixture {output.mixture} > {log} 2>&1
        """


rule coi_mge:
    input:
        reads=expand(f"{BARCODE_DIR}/reads/{{library}}.{{taxon}}.fasta.gz", library=SAMPLES,
                     allow_missing=True),
        cds=rules.cds_fetch.output.table,
    output:
        consensus=f"{BARCODE_DIR}/{{taxon}}/mge/{{taxon}}.mge_consensus.fasta",
        alignment=f"{BARCODE_DIR}/{{taxon}}/mge/{{taxon}}.mge_alignment.fasta",
    params:
        protein=lambda wc: f"{CDS_DIR}/{accession(wc)}.{BARCODE['gene'][0]}.faa",
        dir=f"{BARCODE_DIR}/{{taxon}}/mge",
        code=BARCODE["genetic_code"],
        depth=BARCODE["consensus_min_depth"],
        threshold=BARCODE["consensus_call_fraction"],
    log:
        f"{LOGS}/06_barcode/{{taxon}}.coi_mge.log",
    conda:
        "../envs/mitogeneextractor.yaml"
    shell:
        """
        rm -rf {params.dir}/work
        mkdir -p {params.dir}/work
        zcat {input.reads} > {params.dir}/work/reads.fasta
        MitoGeneExtractor -d {params.dir}/work/reads.fasta -p {params.protein} \\
            -o {params.dir}/work/aln_ -c {params.dir}/work/cons_ -V {params.dir}/work \\
            -C {params.code} -t {params.threshold} \\
            --minSeqCoverageInAlignment_total {params.depth} > {log} 2>&1
        cat {params.dir}/work/cons_*.fas > {output.consensus}
        cat {params.dir}/work/aln_*.fas > {output.alignment}
        rm -rf {params.dir}/work
        """


rule coi_routes:
    input:
        refcoords=rules.consensus.output.refcoords,
        low_depth=rules.consensus.output.low_depth,
        mge=rules.coi_mge.output.consensus,
        cds=rules.cds_fetch.output.table,
    output:
        fasta=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.{BARCODE['gene'][0]}.fasta",
        table=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.{BARCODE['gene'][0]}_routes.tsv",
    params:
        accession=accession,
    log:
        f"{LOGS}/06_barcode/{{taxon}}.coi_routes.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/coi_routes.py --taxon {wildcards.taxon} --accession {params.accession} \\
            --cds {input.cds} --refcoords {input.refcoords} --low-depth {input.low_depth} \\
            --mge {input.mge} --fasta {output.fasta} --table {output.table} > {log} 2>&1
        """


rule bold_search:
    input:
        rules.coi_routes.output.fasta,
    output:
        f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.bold_hits.tsv.gz",
    params:
        curated=BARCODE["bold_curated_db"],
        uncurated=BARCODE["bold_uncurated_db"],
        max_targets=BARCODE["bold_max_target_seqs"],
        query=f"{BARCODE_DIR}/{{taxon}}/{{taxon}}.bold_query.fasta",
    threads: config["threads"]["blastn"]
    log:
        f"{LOGS}/06_barcode/{{taxon}}.bold_search.log",
    conda:
        "../envs/blast.yaml"
    shell:
        """
        # Unknown bases as N, deletions removed. Searched locally; nothing leaves the machine.
        awk '/^>/ {{ print; next }} {{ gsub(/[*]/, ""); gsub(/-/, "N"); print }}' {input} > {params.query}
        : > {output}.tmp
        for db in curated:{params.curated} uncurated:{params.uncurated}; do
            blastn -task megablast -query {params.query} -db "${{db#*:}}" \\
                -max_target_seqs {params.max_targets} -num_threads {threads} \\
                -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qseq sseq" \\
                2>> {log} | sed "s/^/${{db%%:*}}\\t/" >> {output}.tmp
        done
        gzip -c {output}.tmp > {output}
        rm {output}.tmp {params.query}
        """


rule haplogroup:
    input:
        consensus=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.consensus.fasta",
        refcoords=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.consensus_refcoords.fasta",
        libraries=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.libraries.fasta",
        summary=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.consensus_summary.tsv",
        sites=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.site_counts.tsv.gz",
        index=rules.panel_index.output.marker,
    output:
        haplogroups=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.haplogroups.tsv",
        contributors=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.contributors.tsv",
        sites=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.mixed_sites.tsv",
        haplotypes=f"{BARCODE_DIR}/{HUMAN}/{HUMAN}.contributor_haplotypes.fasta",
    params:
        reference=f"{INDEX}/single/{PANEL_ACCESSION[BARCODE['haplogroup_taxon']]}.fasta",
        tree=BARCODE["haplogroup_tree"],
        depth=BARCODE["mixture_depth"],
        fraction=BARCODE["phase_fraction"],
        noise=BARCODE["minor_fraction"],
        work=f"{BARCODE_DIR}/{HUMAN}/haplogrep_work",
    log:
        f"{LOGS}/06_barcode/haplogroup.log",
    conda:
        "../envs/haplogrep.yaml"
    shell:
        """
        # haplogrep3 runs locally against the PhyloTree build it ships with.
        python {SCRIPTS}/haplogroups.py --reference {params.reference} \\
            --pooled {input.refcoords} --consensus {input.consensus} \\
            --libraries {input.libraries} --summary {input.summary} --site-counts {input.sites} \\
            --tree {params.tree} --depth {params.depth} --fraction {params.fraction} \\
            --noise {params.noise} --work {params.work} --haplogroups {output.haplogroups} \\
            --contributors {output.contributors} --sites {output.sites} \\
            --haplotypes {output.haplotypes} > {log} 2>&1
        # haplogrep3 unpacks a bwa JNI library into the working directory; it is not an output.
        rm -rf {params.work} jbwa-data
        """


rule lineage_align:
    input:
        rules.lineage_fetch.output.fasta,
    output:
        f"{BARCODE_DIR}/{{lineage}}/lineage/genomes.aln.fasta",
    threads: config["threads"]["bwa"]
    log:
        f"{LOGS}/06_barcode/{{lineage}}.lineage_align.log",
    conda:
        "../envs/lineage.yaml"
    shell:
        "mafft --thread {threads} --auto {input} > {output} 2> {log}"


rule lineage_place:
    input:
        reads=expand(f"{BARCODE_DIR}/reads/{{library}}.{{lineage}}.fasta.gz", library=SAMPLES,
                     allow_missing=True),
        genomes=rules.lineage_fetch.output.fasta,
        table=rules.lineage_fetch.output.table,
        alignment=rules.lineage_align.output,
    output:
        clusters=f"{BARCODE_DIR}/{{lineage}}/lineage/clusters.tsv",
        genomes=f"{BARCODE_DIR}/{{lineage}}/lineage/genome_scores.tsv",
        reads=f"{BARCODE_DIR}/{{lineage}}/lineage/read_support.tsv",
        summary=f"{BARCODE_DIR}/{{lineage}}/lineage/cluster_support.tsv",
    params:
        aln=config["mapping"]["bwa_aln"],
        min_mapq=config["mapping"]["min_mapq"],
        species=lambda wc: LINEAGE_TARGETS[wc.lineage],
        reference=lambda wc: PANEL_ACCESSION[LINEAGE_TARGETS[wc.lineage]],
        cut=lambda wc: BARCODE["lineage"][LINEAGE_TARGETS[wc.lineage]]["cluster_distance"],
        work=f"{BARCODE_DIR}/{{lineage}}/lineage/work",
    threads: config["threads"]["bwa"]
    log:
        f"{LOGS}/06_barcode/{{lineage}}.lineage_place.log",
    conda:
        "../envs/lineage.yaml"
    shell:
        """
        python {SCRIPTS}/lineage_place.py --reads {input.reads} --genomes {input.genomes} \\
            --table {input.table} --alignment {input.alignment} --species {params.species:q} \\
            --reference {params.reference} --cut {params.cut} --aln-params {params.aln:q} \\
            --min-mapq {params.min_mapq} --threads {threads} --work {params.work} \\
            --clusters {output.clusters} --genome-scores {output.genomes} \\
            --read-support {output.reads} --summary {output.summary} > {log} 2>&1
        rm -rf {params.work}
        """


rule recovery_by_length:
    input:
        genomes=rules.lineage_fetch.output.fasta,
        scores=rules.lineage_place.output.genomes,
        clusters=rules.lineage_place.output.clusters,
        index=rules.panel_index.output.marker,
        lengths=expand(rules.damage_profile.output.lengths, library=SAMPLES),
    output:
        table=f"{BARCODE_DIR}/{{lineage}}/recovery_by_length.tsv",
        predicted=f"{BARCODE_DIR}/{{lineage}}/recovery_predicted_lengths.tsv",
    params:
        species=lambda wc: LINEAGE_TARGETS[wc.lineage],
        reference=lambda wc: PANEL_ACCESSION[LINEAGE_TARGETS[wc.lineage]],
        panel=f"{INDEX}/panel.fasta",
        single=f"{INDEX}/single",
        aln=config["mapping"]["bwa_aln"],
        margin=config["mapping"]["panel_margin_mismatches"],
        min_mapq=config["mapping"]["min_mapq"],
        sim=BARCODE["recovery_simulation"],
        work=f"{BARCODE_DIR}/{{lineage}}/recovery_work",
    threads: config["threads"]["bwa"]
    log:
        f"{LOGS}/06_barcode/{{lineage}}.recovery_by_length.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/simulate_recovery.py --genomes {input.genomes} \\
            --genome-scores {input.scores} --clusters {input.clusters} \\
            --species {params.species:q} --reference {params.reference} \\
            --panel {params.panel} --single {params.single} --aln-params {params.aln:q} \\
            --margin {params.margin} --min-mapq {params.min_mapq} --threads {threads} \\
            --min-length {params.sim[min_length]} --max-length {params.sim[max_length]} \\
            --step {params.sim[step]} --per-length {params.sim[per_length]} \\
            --seed {params.sim[seed]} --observed-lengths {input.lengths} \\
            --work {params.work} \\
            --out {output.table} --predicted {output.predicted} > {log} 2>&1
        rm -rf {params.work}
        """


rule coi_quality:
    input:
        expand(rules.coi_routes.output.fasta, taxon=TARGETS),
    output:
        f"{BARCODE_DIR}/coi_quality.tsv",
    params:
        code=BARCODE["genetic_code"],
        bold=BARCODE["bold_barcode_standard"],
        validator=BARCODE["validator_min_length"],
    log:
        f"{LOGS}/06_barcode/coi_quality.log",
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python {SCRIPTS}/coi_quality.py --fasta {input} --table {params.code} \\
            --bold-min-length {params.bold} --validator-min-length {params.validator} \\
            --out {output} > {log} 2>&1
        """


rule gene_coverage_figure:
    input:
        sites=expand(rules.consensus_stats.output.sites, taxon=TARGETS),
        cds=rules.cds_fetch.output.table,
    output:
        f"{BARCODE_DIR}/gene_coverage.png",
    params:
        taxa=list(TARGETS),
        accessions=[PANEL_ACCESSION[t] for t in TARGETS.values()],
        depth=BARCODE["consensus_min_depth"],
    log:
        f"{LOGS}/06_barcode/gene_coverage_figure.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/plot_gene_coverage.py --site-counts {input.sites} \\
            --taxa {params.taxa} --accessions {params.accessions} --cds {input.cds} \\
            --min-depth {params.depth} --out {output} > {log} 2>&1
        """


rule barcode_summary:
    input:
        summaries=expand(rules.consensus_stats.output.summary, taxon=TARGETS),
        mixtures=expand(rules.consensus_stats.output.mixture, taxon=TARGETS),
        routes=expand(rules.coi_routes.output.table, taxon=TARGETS),
        hits=expand(rules.bold_search.output, taxon=TARGETS),
        haplogroups=rules.haplogroup.output.haplogroups,
        lineage=expand(rules.lineage_place.output.summary, lineage=LINEAGE_TARGETS),
        genome_scores=expand(rules.lineage_place.output.genomes, lineage=LINEAGE_TARGETS),
        recovery=expand(rules.recovery_by_length.output.table, lineage=LINEAGE_TARGETS),
        predicted=expand(rules.recovery_by_length.output.predicted, lineage=LINEAGE_TARGETS),
    output:
        bold=f"{BARCODE_DIR}/bold_summary.tsv",
        figure=f"{BARCODE_DIR}/step6_summary.png",
    params:
        taxa=list(TARGETS),
    log:
        f"{LOGS}/06_barcode/barcode_summary.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        python {SCRIPTS}/barcode_summary.py --taxa {params.taxa} --hits {input.hits} \\
            --summaries {input.summaries} --lineage {input.lineage} \\
            --genome-scores {input.genome_scores} --recovery {input.recovery} \\
            --predicted {input.predicted} --bold-out {output.bold} \\
            --figure {output.figure} > {log} 2>&1
        """
