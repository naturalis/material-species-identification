# Step 7: phylogenetic placement. Stub.


rule placement:
    input:
        expand(rules.coi_routes.output.fasta, taxon=TARGETS),
        expand(rules.consensus.output.consensus, taxon=TARGETS),
    output:
        tree=f"{RESULTS}/07_placement/tree.treefile",
        sites=f"{RESULTS}/07_placement/diagnostic_sites.tsv",
    shell:
        STUB
