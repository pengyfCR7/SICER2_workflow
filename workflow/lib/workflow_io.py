import csv
import os
import re
from pathlib import Path


NAME_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+$")
SAMPLE_COLUMNS = {"group", "sample", "treatment_bam", "control_bam", "layout"}
CONTRAST_COLUMNS = {"contrast", "test", "reference"}
TAIR10_NUCLEAR_LENGTHS = {
    "Chr1": 30427671,
    "Chr2": 19698289,
    "Chr3": 23459830,
    "Chr4": 18585056,
    "Chr5": 26975502,
}


def _check_name(value, field, line_number):
    if not value or not NAME_PATTERN.fullmatch(value):
        raise ValueError(
            f"Invalid {field} {value!r} at line {line_number}; "
            "use letters, numbers, '.', '_' or '-'"
        )


def _file_identity(path):
    resolved = Path(path).resolve(strict=True)
    stat = resolved.stat()
    return (stat.st_dev, stat.st_ino), str(resolved)


def _unique_files(paths):
    unique = []
    seen = set()
    for path in paths:
        identity, resolved = _file_identity(path)
        if identity not in seen:
            seen.add(identity)
            unique.append(resolved)
    return unique


def load_samples(path, require_files=True):
    table = Path(path).resolve()
    with table.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = SAMPLE_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing sample columns: {sorted(missing)}")

        rows = []
        seen_samples = set()
        treatment_identities = {}
        for line_number, row in enumerate(reader, start=2):
            if not any((value or "").strip() for value in row.values()):
                continue
            group = (row.get("group") or "").strip()
            sample = (row.get("sample") or "").strip()
            treatment = (row.get("treatment_bam") or "").strip()
            control = (row.get("control_bam") or "").strip()
            layout = (row.get("layout") or "").strip().upper()
            _check_name(group, "group", line_number)
            _check_name(sample, "sample", line_number)
            if sample in seen_samples:
                raise ValueError(f"Duplicate sample {sample!r} at line {line_number}")
            seen_samples.add(sample)
            if layout != "SE":
                raise ValueError(
                    f"SICER2_workflow v1 supports only layout=SE; "
                    f"sample {sample!r} uses {layout!r}"
                )
            if not treatment:
                raise ValueError(f"Missing treatment_bam at line {line_number}")

            def resolve_input(value):
                candidate = Path(value).expanduser()
                if not candidate.is_absolute():
                    candidate = table.parent / candidate
                if require_files:
                    return _file_identity(candidate)[1]
                return str(candidate.resolve(strict=False))

            treatment = resolve_input(treatment)
            control = resolve_input(control) if control else ""
            if require_files:
                identity, _ = _file_identity(treatment)
                if identity in treatment_identities:
                    other = treatment_identities[identity]
                    raise ValueError(
                        f"Treatment BAM is reused by samples {other!r} and {sample!r}: "
                        f"{treatment}"
                    )
                treatment_identities[identity] = sample
            rows.append(
                {
                    "group": group,
                    "sample": sample,
                    "treatment_bam": treatment,
                    "control_bam": control,
                    "layout": layout,
                }
            )

    if not rows:
        raise ValueError(f"No samples found in {table}")

    groups = list(dict.fromkeys(row["group"] for row in rows))
    group_rows = {group: [row for row in rows if row["group"] == group] for group in groups}
    group_treatments = {}
    group_controls = {}
    for group, members in group_rows.items():
        controls = [row["control_bam"] for row in members]
        if any(controls) and not all(controls):
            raise ValueError(f"Group {group!r} mixes rows with and without control_bam")
        group_treatments[group] = [row["treatment_bam"] for row in members]
        group_controls[group] = _unique_files(controls) if controls and controls[0] else []

    return {
        "table": str(table),
        "rows": rows,
        "groups": groups,
        "group_rows": group_rows,
        "group_treatments": group_treatments,
        "group_controls": group_controls,
    }


def load_contrasts(path, groups):
    table = Path(path).resolve()
    if not table.exists():
        raise FileNotFoundError(f"Contrast table not found: {table}")
    with table.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = CONTRAST_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing contrast columns: {sorted(missing)}")
        contrasts = []
        seen = set()
        for line_number, row in enumerate(reader, start=2):
            if not any((value or "").strip() for value in row.values()):
                continue
            name = (row.get("contrast") or "").strip()
            test = (row.get("test") or "").strip()
            reference = (row.get("reference") or "").strip()
            _check_name(name, "contrast", line_number)
            if name in seen:
                raise ValueError(f"Duplicate contrast {name!r} at line {line_number}")
            seen.add(name)
            if test not in groups or reference not in groups:
                raise ValueError(
                    f"Contrast {name!r} references unknown group(s): {test!r}, {reference!r}"
                )
            if test == reference:
                raise ValueError(f"Contrast {name!r} compares a group with itself")
            contrasts.append({"contrast": name, "test": test, "reference": reference})
    return {row["contrast"]: row for row in contrasts}


def validate_parameters(call_cfg, diff_cfg):
    redundancy = int(call_cfg["redundancy_threshold"])
    window = int(call_cfg["window_size"])
    fragment = int(call_cfg["fragment_size"])
    egf = float(call_cfg["effective_genome_fraction"])
    gap = int(call_cfg["gap_size"])
    call_fdr = float(call_cfg["false_discovery_rate"])
    diff_fdr = float(diff_cfg["false_discovery_rate_df"])
    min_fc = float(diff_cfg["min_fold_change"])
    if redundancy < 1:
        raise ValueError("redundancy_threshold must be >= 1")
    if window <= 0 or gap < 0 or gap % window != 0:
        raise ValueError("window_size must be positive and gap_size must be a non-negative multiple")
    if fragment <= 0:
        raise ValueError("fragment_size must be positive for SE data")
    if not 0 < egf <= 1:
        raise ValueError("effective_genome_fraction must be in (0, 1]")
    if not 0 < call_fdr <= 1:
        raise ValueError("false_discovery_rate must be in (0, 1]")
    if not 0 < diff_fdr <= 1:
        raise ValueError("false_discovery_rate_df must be in (0, 1]")
    if min_fc < 1:
        raise ValueError("min_fold_change must be >= 1")


def validate_contrast_controls(contrasts, group_controls):
    for name, contrast in contrasts.items():
        test_has = bool(group_controls[contrast["test"]])
        reference_has = bool(group_controls[contrast["reference"]])
        if test_has != reference_has:
            raise ValueError(
                f"Contrast {name!r} mixes a group with control and a group without control"
            )


def validate_tair10_nuclear_chrom_sizes(path):
    observed = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[1].isdigit():
                raise ValueError(f"Invalid chrom sizes line {line_number}: {line.rstrip()}")
            chrom, length = fields[0], int(fields[1])
            if chrom in observed:
                raise ValueError(f"Duplicate chromosome in chrom sizes: {chrom}")
            observed[chrom] = length
    if observed != TAIR10_NUCLEAR_LENGTHS:
        raise ValueError(
            "chrom_sizes must contain exactly TAIR10 Chr1-Chr5 with assembled lengths; "
            f"observed {observed}"
        )
