rule pool_treatment_bed:
    input:
        lambda wc: SAMPLES["group_treatments"][wc.group]
    output:
        temp(str(WORK_DIR / "pooled/{group}.treatment.bed"))
    log:
        str(LOG_DIR / "pool/{group}.treatment.log")
    shell:
        """
        mkdir -p $(dirname {output:q}) $(dirname {log:q})
        python3 {POOL_SCRIPT:q} \
          --samtools {COMMANDS[samtools]:q} \
          --bedtools {COMMANDS[bedtools]:q} \
          --chrom-sizes {CHROM_SIZES:q} \
          --output {output:q} \
          --bam {input:q} > {log:q} 2>&1
        """


rule pool_control_bed:
    input:
        lambda wc: SAMPLES["group_controls"][wc.group]
    output:
        temp(str(WORK_DIR / "pooled/{group}.control.bed"))
    log:
        str(LOG_DIR / "pool/{group}.control.log")
    run:
        if not input:
            raise ValueError(f"Group {wildcards.group} has no control")
        shell(
            "mkdir -p $(dirname {output:q}) $(dirname {log:q}); "
            "python3 {POOL_SCRIPT:q} "
            "--samtools {COMMANDS[samtools]:q} --bedtools {COMMANDS[bedtools]:q} "
            "--chrom-sizes {CHROM_SIZES:q} --output {output:q} "
            "--bam {input:q} > {log:q} 2>&1"
        )


def call_inputs(wc):
    files = {"treatment": treatment_bed(wc.group)}
    if SAMPLES["group_controls"][wc.group]:
        files["control"] = control_bed(wc.group)
    return files


rule sicer:
    input:
        unpack(call_inputs)
    output:
        directory(str(RESULT_DIR / "sicer2/{group}"))
    threads:
        THREADS
    log:
        str(LOG_DIR / "sicer/{group}.log")
    params:
        control=lambda wc, input: f"--control_file {input.control}" if hasattr(input, "control") else ""
    shell:
        """
        mkdir -p $(dirname {output:q}) $(dirname {log:q})
        {COMMANDS[sicer2_python]:q} {RUN_SCRIPT:q} call \
          --chrom-sizes {CHROM_SIZES:q} \
          --treatment_file {input.treatment:q} {params.control} \
          --output_directory {output:q} \
          --redundancy_threshold {CALL_CFG[redundancy_threshold]} \
          --window_size {CALL_CFG[window_size]} \
          --fragment_size {CALL_CFG[fragment_size]} \
          --effective_genome_fraction {CALL_CFG[effective_genome_fraction]} \
          --gap_size {CALL_CFG[gap_size]} \
          --false_discovery_rate {CALL_CFG[false_discovery_rate]} \
          --cpu {threads} {CALL_CFG[extra]} > {log:q} 2>&1
        """


def diff_inputs(wc):
    test = CONTRASTS[wc.contrast]["test"]
    reference = CONTRASTS[wc.contrast]["reference"]
    test_has_control = bool(SAMPLES["group_controls"][test])
    reference_has_control = bool(SAMPLES["group_controls"][reference])
    if test_has_control != reference_has_control:
        raise ValueError(f"Contrast {wc.contrast} mixes groups with and without control")

    files = {"test": treatment_bed(test), "reference": treatment_bed(reference)}
    if test_has_control:
        files["test_control"] = control_bed(test)
        files["reference_control"] = control_bed(reference)
    return files


rule sicer_df:
    input:
        unpack(diff_inputs)
    output:
        directory(str(RESULT_DIR / "sicer2_diff/{contrast}"))
    threads:
        THREADS
    log:
        str(LOG_DIR / "sicer_df/{contrast}.log")
    params:
        controls=lambda wc, input: (
            f"--control_file {input.test_control} {input.reference_control}"
            if hasattr(input, "test_control") else ""
        )
    shell:
        """
        mkdir -p $(dirname {output:q}) $(dirname {log:q})
        {COMMANDS[sicer2_python]:q} {RUN_SCRIPT:q} diff \
          --chrom-sizes {CHROM_SIZES:q} \
          --treatment_file {input.test:q} {input.reference:q} {params.controls} \
          --output_directory {output:q} \
          --redundancy_threshold {CALL_CFG[redundancy_threshold]} \
          --window_size {CALL_CFG[window_size]} \
          --fragment_size {CALL_CFG[fragment_size]} \
          --effective_genome_fraction {CALL_CFG[effective_genome_fraction]} \
          --gap_size {CALL_CFG[gap_size]} \
          --false_discovery_rate {CALL_CFG[false_discovery_rate]} \
          --false_discovery_rate_df {DIFF_CFG[false_discovery_rate_df]} \
          --cpu {threads} {DIFF_CFG[extra]} > {log:q} 2>&1
        """
