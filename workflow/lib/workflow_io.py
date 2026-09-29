import csv
import re
from pathlib import Path


NAME_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+$")
SAMPLE_COLUMNS = {"group", "sample", "treatment_bam", "control_bam", "layout"}
CONTRAST_COLUMNS = {"contrast", "test", "reference"}


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


def load_samples(path):
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
                return _file_identity(candidate)[1]

            treatment = resolve_input(treatment)
            control = resolve_input(control) if control else ""
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
        "groups": groups,
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
