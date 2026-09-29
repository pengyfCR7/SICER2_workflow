#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
from pathlib import Path


def line_count(path, has_header=True):
    with open(path, encoding="utf-8") as handle:
        count = sum(1 for line in handle if line.strip())
    return max(count - int(has_header), 0)


def package_version(python):
    code = "from importlib.metadata import version; print(version('SICER'))"
    return subprocess.check_output([python, "-c", code], text=True).strip()


def file_record(path):
    resolved = Path(path).resolve(strict=True)
    stat = resolved.stat()
    return f"{stat.st_dev}:{stat.st_ino}:{resolved}"


def result_counts(files):
    return {
        name: line_count(path, has_header=(name != "islands"))
        for name, path in files.items()
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    with open(args.spec, encoding="utf-8") as handle:
        spec = json.load(handle)

    fields = [
        "record_type", "id", "test", "reference", "treatment_files", "control_files",
        "file_identities", "parameters", "sicer2_module", "sicer_package",
        "result_files", "result_counts",
    ]
    version = package_version(spec["sicer2_python"])
    rows = []
    for group in spec["groups"]:
        rows.append(
            {
                "record_type": "group_call",
                "id": group["group"],
                "test": "",
                "reference": "",
                "treatment_files": ";".join(group["treatment_files"]),
                "control_files": ";".join(group["control_files"]),
                "file_identities": ";".join(
                    file_record(path)
                    for path in group["treatment_files"] + group["control_files"]
                ),
                "parameters": json.dumps(spec["call_parameters"], sort_keys=True),
                "sicer2_module": spec["sicer2_module"],
                "sicer_package": version,
                "result_files": json.dumps(group["result_files"], sort_keys=True),
                "result_counts": json.dumps(result_counts(group["result_files"]), sort_keys=True),
            }
        )
    for contrast in spec["contrasts"]:
        rows.append(
            {
                "record_type": "differential",
                "id": contrast["contrast"],
                "test": contrast["test"],
                "reference": contrast["reference"],
                "treatment_files": (
                    "test=" + ";".join(contrast["test_treatment_files"]) +
                    "|reference=" + ";".join(contrast["reference_treatment_files"])
                ),
                "control_files": (
                    "test=" + ";".join(contrast["test_control_files"]) +
                    "|reference=" + ";".join(contrast["reference_control_files"])
                ),
                "file_identities": ";".join(
                    file_record(path) for path in (
                        contrast["test_treatment_files"] +
                        contrast["reference_treatment_files"] +
                        contrast["test_control_files"] +
                        contrast["reference_control_files"]
                    )
                ),
                "parameters": json.dumps(
                    {
                        "call": spec["call_parameters"],
                        "differential": spec["diff_parameters"],
                    },
                    sort_keys=True,
                ),
                "sicer2_module": spec["sicer2_module"],
                "sicer_package": version,
                "result_files": json.dumps(contrast["result_files"], sort_keys=True),
                "result_counts": json.dumps(result_counts(contrast["result_files"]), sort_keys=True),
            }
        )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
