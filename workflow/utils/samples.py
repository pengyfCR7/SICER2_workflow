import pandas as pd


def load_samples(sample_table):
    samples = pd.read_csv(sample_table, sep="\t", dtype=str).fillna("")
    columns = ["group", "sample", "treatment_bam", "control_bam", "layout"]

    missing = set(columns) - set(samples.columns)
    if missing:
        raise ValueError(f"Missing columns in sample table: {missing}")

    samples = samples[columns].apply(lambda column: column.str.strip())
    samples["layout"] = samples["layout"].str.upper()

    if samples["sample"].duplicated().any():
        raise ValueError("Duplicate sample names")
    if not samples["layout"].eq("SE").all():
        raise ValueError("SICER2_workflow supports only layout=SE")
    if samples[["group", "sample", "treatment_bam"]].eq("").any().any():
        raise ValueError("group, sample and treatment_bam are required")

    for group, rows in samples.groupby("group", sort=False):
        if rows["control_bam"].eq("").nunique() > 1:
            raise ValueError(f"Group {group} mixes samples with and without control")

    return samples
