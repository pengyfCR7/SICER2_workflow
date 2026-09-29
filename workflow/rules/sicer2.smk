rule pool_treatment_bed:
    input:
        lambda wildcards: SAMPLE_INFO["group_treatments"][wildcards.group]
    output:
        bed=temp(str(WORK_DIR / "pooled/{group}.treatment.bed")),
        stats=temp(str(WORK_DIR / "pooled/{group}.treatment.stats.txt"))
    log:
        str(LOG_DIR / "pool/{group}.treatment.log")
    params:
        samtools=COMMANDS["samtools"],
        bedtools=COMMANDS["bedtools"],
        script=str(POOL_SCRIPT),
        chrom_sizes=CHROM_SIZES
    shell:
        """
        mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
        python3 {params.script:q} \
          --samtools {params.samtools:q} \
          --bedtools {params.bedtools:q} \
          --chrom-sizes {params.chrom_sizes:q} \
          --output {output.bed:q} \
          --stats {output.stats:q} \
          --bam {input:q} > {log:q} 2>&1
        """


rule pool_control_bed:
    input:
        lambda wildcards: SAMPLE_INFO["group_controls"][wildcards.group]
    output:
        bed=temp(str(WORK_DIR / "pooled/{group}.control.bed")),
        stats=temp(str(WORK_DIR / "pooled/{group}.control.stats.txt"))
    log:
        str(LOG_DIR / "pool/{group}.control.log")
    params:
        samtools=COMMANDS["samtools"],
        bedtools=COMMANDS["bedtools"],
        script=str(POOL_SCRIPT),
        chrom_sizes=CHROM_SIZES
    run:
        if wildcards.group not in GROUPS_WITH_CONTROL:
            raise ValueError(f"Group {wildcards.group!r} has no control")
        shell(
            "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}); "
            "python3 {params.script:q} "
            "--samtools {params.samtools:q} --bedtools {params.bedtools:q} "
            "--chrom-sizes {params.chrom_sizes:q} --output {output.bed:q} "
            "--stats {output.stats:q} --bam {input:q} > {log:q} 2>&1"
        )


def call_inputs(wildcards):
    inputs = {"treatment": pooled_treatment(wildcards.group)}
    if wildcards.group in GROUPS_WITH_CONTROL:
        inputs["control"] = pooled_control(wildcards.group)
    return inputs


rule sicer2_call:
    input:
        unpack(call_inputs)
    output:
        islands=str(SICER_DIR / "{group}/{group}.islands.bed"),
        summary=str(SICER_DIR / "{group}/{group}.islands.summary.txt"),
        raw=directory(str(SICER_DIR / "{group}/raw"))
    threads:
        THREADS
    log:
        str(LOG_DIR / "sicer2_call/{group}.log")
    params:
        python=COMMANDS["sicer2_python"],
        script=str(RUN_SCRIPT),
        chrom_sizes=CHROM_SIZES,
        control=lambda wildcards, input: f"--control_file {input.control}" if hasattr(input, "control") else "",
        rt=int(CALL_CFG["redundancy_threshold"]),
        window=int(CALL_CFG["window_size"]),
        fragment=int(CALL_CFG["fragment_size"]),
        egf=float(CALL_CFG["effective_genome_fraction"]),
        gap=int(CALL_CFG["gap_size"]),
        fdr=float(CALL_CFG["false_discovery_rate"]),
        extra=CALL_CFG.get("extra", "")
    shell:
        """
        mkdir -p {output.raw:q} $(dirname {log:q})
        {params.python:q} {params.script:q} call \
          --chrom-sizes {params.chrom_sizes:q} \
          --treatment_file {input.treatment:q} {params.control} \
          --output-directory {output.raw:q} \
          --redundancy_threshold {params.rt} \
          --window_size {params.window} \
          --fragment_size {params.fragment} \
          --effective_genome_fraction {params.egf} \
          --gap_size {params.gap} \
          --false_discovery_rate {params.fdr} \
          --cpu {threads} \
          --stable-islands {output.islands:q} \
          --stable-summary {output.summary:q} \
          {params.extra} > {log:q} 2>&1
        """


def diff_inputs(wildcards):
    test, reference = contrast_groups(wildcards.contrast)
    if (test in GROUPS_WITH_CONTROL) != (reference in GROUPS_WITH_CONTROL):
        raise ValueError(f"Contrast {wildcards.contrast} mixes groups with and without control")
    inputs = {
        "test": pooled_treatment(test),
        "reference": pooled_treatment(reference),
    }
    if wildcards.contrast in CONTRASTS_WITH_CONTROL:
        inputs["test_control"] = pooled_control(test)
        inputs["reference_control"] = pooled_control(reference)
    return inputs


rule sicer2_diff:
    input:
        unpack(diff_inputs)
    output:
        all=str(DIFF_DIR / "{contrast}/{contrast}.all_islands.txt"),
        increased_fdr=str(DIFF_DIR / "{contrast}/{contrast}.increased.fdr.txt"),
        decreased_fdr=str(DIFF_DIR / "{contrast}/{contrast}.decreased.fdr.txt"),
        increased_filtered=str(DIFF_DIR / "{contrast}/{contrast}.increased.filtered.txt"),
        decreased_filtered=str(DIFF_DIR / "{contrast}/{contrast}.decreased.filtered.txt"),
        raw=directory(str(DIFF_DIR / "{contrast}/raw"))
    threads:
        THREADS
    log:
        str(LOG_DIR / "sicer2_diff/{contrast}.log")
    params:
        python=COMMANDS["sicer2_python"],
        script=str(RUN_SCRIPT),
        chrom_sizes=CHROM_SIZES,
        controls=lambda wildcards, input: (
            f"--control_file {input.test_control} {input.reference_control}"
            if hasattr(input, "test_control") else ""
        ),
        rt=int(CALL_CFG["redundancy_threshold"]),
        window=int(CALL_CFG["window_size"]),
        fragment=int(CALL_CFG["fragment_size"]),
        egf=float(CALL_CFG["effective_genome_fraction"]),
        gap=int(CALL_CFG["gap_size"]),
        fdr=float(CALL_CFG["false_discovery_rate"]),
        fdr_df=float(DIFF_CFG["false_discovery_rate_df"]),
        min_fc=float(DIFF_CFG["min_fold_change"]),
        extra=DIFF_CFG.get("extra", "")
    shell:
        """
        mkdir -p {output.raw:q} $(dirname {log:q})
        {params.python:q} {params.script:q} diff \
          --chrom-sizes {params.chrom_sizes:q} \
          --treatment_file {input.test:q} {input.reference:q} {params.controls} \
          --output-directory {output.raw:q} \
          --redundancy_threshold {params.rt} \
          --window_size {params.window} \
          --fragment_size {params.fragment} \
          --effective_genome_fraction {params.egf} \
          --gap_size {params.gap} \
          --false_discovery_rate {params.fdr} \
          --false_discovery_rate_df {params.fdr_df} \
          --min-fold-change {params.min_fc} \
          --cpu {threads} \
          --all-islands {output.all:q} \
          --increased-fdr {output.increased_fdr:q} \
          --decreased-fdr {output.decreased_fdr:q} \
          --increased-filtered {output.increased_filtered:q} \
          --decreased-filtered {output.decreased_filtered:q} \
          {params.extra} > {log:q} 2>&1
        """


rule run_manifest:
    input:
        calls=expand(str(SICER_DIR / "{group}/{group}.islands.summary.txt"), group=GROUPS),
        diffs=expand(str(DIFF_DIR / "{contrast}/{contrast}.all_islands.txt"), contrast=CONTRASTS),
        spec=str(MANIFEST_SPEC)
    output:
        str(RESULT_DIR / "run_manifest.txt")
    params:
        python=COMMANDS["sicer2_python"],
        script=str(MANIFEST_SCRIPT)
    shell:
        "{params.python:q} {params.script:q} --spec {input.spec:q} --output {output:q}"
