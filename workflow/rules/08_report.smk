# Step 8: report. Stub.


rule report:
    input:
        f"{RESULTS}/05_auth/contammix.tsv",
        f"{RESULTS}/07_placement/tree.treefile",
    output:
        f"{RESULTS}/report.md",
    shell:
        STUB
