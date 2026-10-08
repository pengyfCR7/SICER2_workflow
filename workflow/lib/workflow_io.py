from pathlib import Path

import pandas as pd


def load_samples(path):
    """Read and check the replicate table."""
    path = Path(path).resolve()
    columns = ["group", "sample", "treatment_bam", "control_bam", "layout"]
    samples = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    missing = [column for column in columns if column not in samples]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")

    samples = samples[columns].apply(lambda column: column.str.strip())
    if samples.empty:
        raise ValueError(f"No samples found in {path}")
    if samples[["group", "sample", "treatment_bam"]].eq("").any().any():
        raise ValueError("group, sample and treatment_bam are required")
    if samples["sample"].duplicated().any():
        raise ValueError("sample names must be unique")
    if not samples["layout"].str.upper().eq("SE").all():
        raise ValueError("SICER2_workflow v1 supports only layout=SE")
    samples["layout"] = "SE"

    # Resolve paths once, so repeated controls and symlinks become the same path.
    for column in ["treatment_bam", "control_bam"]:
        resolved = []
        for value in samples[column]:
            if not value:
                resolved.append("")
                continue
            bam = Path(value).expanduser()
            if not bam.is_absolute():
                bam = path.parent / bam
            resolved.append(str(bam.resolve(strict=True)))
        samples[column] = resolved

    for group, rows in samples.groupby("group", sort=False):
        has_control = rows["control_bam"].ne("")
        if has_control.any() and not has_control.all():
            raise ValueError(f"Group {group} mixes samples with and without control")

    return samples


def load_contrasts(path, groups):
    """Read contrasts and confirm that their groups exist."""
    path = Path(path).resolve()
    columns = ["contrast", "test", "reference"]
    contrasts = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    missing = [column for column in columns if column not in contrasts]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")

    contrasts = contrasts[columns].apply(lambda column: column.str.strip())
    if contrasts["contrast"].duplicated().any():
        raise ValueError("contrast names must be unique")
    if contrasts[["contrast", "test", "reference"]].eq("").any().any():
        raise ValueError("contrast, test and reference are required")

    unknown = set(contrasts["test"]) | set(contrasts["reference"])
    unknown -= set(groups)
    if unknown:
        raise ValueError(f"Unknown groups in contrasts: {sorted(unknown)}")

    return contrasts
