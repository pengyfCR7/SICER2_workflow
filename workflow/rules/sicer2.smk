rule pool_treatment_bed:
    input:
        lambda wc: GROUP_TREATMENTS[wc.group]
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
        lambda wc: GROUP_CONTROLS[wc.group]
    output:
        temp(str(WORK_DIR / "pooled/{group}.control.bed"))
    log:
        str(LOG_DIR / "pool/{group}.control.log")
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


rule sicer:
    input:
        unpack(lambda wc: CALL_INPUTS[wc.group])
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


rule sicer_df:
    input:
        unpack(lambda wc: DIFF_INPUTS[wc.contrast])
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
