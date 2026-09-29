import csv
from pathlib import Path


def read_table(path, required_columns):
    path = Path(path).resolve()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required_columns) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]


def real_file(path, table_dir):
    path = Path(path).expanduser()
    if not path.is_absolute():
        path = table_dir / path
    return str(path.resolve(strict=True))


def unique_files(paths):
    unique = []
    seen = set()
    for path in paths:
        stat = Path(path).stat()
        identity = (stat.st_dev, stat.st_ino)
        if identity not in seen:
            seen.add(identity)
            unique.append(path)
    return unique


def load_samples(path):
    table = Path(path).resolve()
    rows = read_table(
        table,
        ["group", "sample", "treatment_bam", "control_bam", "layout"],
    )
    if not rows:
        raise ValueError(f"No samples found in {table}")
    if len({row["sample"] for row in rows}) != len(rows):
        raise ValueError("sample names must be unique")

    for row in rows:
        if not row["group"] or not row["sample"] or not row["treatment_bam"]:
            raise ValueError("group, sample and treatment_bam are required")
        if row["layout"].upper() != "SE":
            raise ValueError("SICER2_workflow v1 supports only layout=SE")
        row["treatment_bam"] = real_file(row["treatment_bam"], table.parent)
        if row["control_bam"]:
            row["control_bam"] = real_file(row["control_bam"], table.parent)

    groups = list(dict.fromkeys(row["group"] for row in rows))
    treatments = {}
    controls = {}
    for group in groups:
        members = [row for row in rows if row["group"] == group]
        treatments[group] = [row["treatment_bam"] for row in members]
        group_controls = [row["control_bam"] for row in members]
        if any(group_controls) and not all(group_controls):
            raise ValueError(f"Group {group} mixes samples with and without control")
        controls[group] = unique_files(group_controls) if all(group_controls) else []

    return {
        "groups": groups,
        "group_treatments": treatments,
        "group_controls": controls,
    }


def load_contrasts(path, groups):
    rows = read_table(path, ["contrast", "test", "reference"])
    contrasts = {}
    for row in rows:
        name = row["contrast"]
        if name in contrasts:
            raise ValueError(f"Duplicate contrast: {name}")
        if row["test"] not in groups or row["reference"] not in groups:
            raise ValueError(f"Unknown group in contrast: {name}")
        contrasts[name] = row
    return contrasts
