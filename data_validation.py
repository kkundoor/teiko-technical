import csv
from pathlib import Path


COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

CELL_COLUMNS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

SUBJECT_FIELDS = [
    "project",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
]


class DataValidationError(ValueError):
    pass


def _required_text(row, column, line_number):
    value = row.get(column)

    if value is None or value.strip() == "":
        raise DataValidationError(
            f"line {line_number}: {column} is required"
        )

    return value.strip()


def _integer(row, column, line_number):
    raw = row.get(column)

    if raw is None or raw.strip() == "":
        raise DataValidationError(
            f"line {line_number}: {column} is required"
        )

    try:
        return int(raw)
    except ValueError as exc:
        raise DataValidationError(
            f"line {line_number}: {column} must be an integer"
        ) from exc


def load_validated_rows(path):
    path = Path(path)

    if not path.exists():
        raise DataValidationError(f"input file not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []

        if len(header) != len(set(header)):
            raise DataValidationError(
                "input contains duplicate column names"
            )

        expected = set(COLUMNS)
        actual = set(header)

        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected)

        if missing or unexpected:
            details = []

            if missing:
                details.append(
                    f"missing columns: {', '.join(missing)}"
                )

            if unexpected:
                details.append(
                    f"unexpected columns: {', '.join(unexpected)}"
                )

            raise DataValidationError("; ".join(details))

        rows = []
        seen_samples = set()
        subject_metadata = {}

        for line_number, source_row in enumerate(reader, start=2):
            if None in source_row:
                raise DataValidationError(
                    f"line {line_number}: too many values"
                )

            row = {
                "project": _required_text(
                    source_row, "project", line_number
                ),
                "subject": _required_text(
                    source_row, "subject", line_number
                ),
                "condition": _required_text(
                    source_row, "condition", line_number
                ),
                "age": _integer(
                    source_row, "age", line_number
                ),
                "sex": _required_text(
                    source_row, "sex", line_number
                ),
                "treatment": _required_text(
                    source_row, "treatment", line_number
                ),
                "response": (
                    source_row["response"].strip()
                    if source_row["response"]
                    and source_row["response"].strip()
                    else None
                ),
                "sample": _required_text(
                    source_row, "sample", line_number
                ),
                "sample_type": _required_text(
                    source_row, "sample_type", line_number
                ),
                "time_from_treatment_start": _integer(
                    source_row,
                    "time_from_treatment_start",
                    line_number,
                ),
            }

            if row["age"] < 0:
                raise DataValidationError(
                    f"line {line_number}: age cannot be negative"
                )

            for column in CELL_COLUMNS:
                count = _integer(
                    source_row, column, line_number
                )

                if count < 0:
                    raise DataValidationError(
                        f"line {line_number}: "
                        f"{column} cannot be negative"
                    )

                row[column] = count

            if sum(row[column] for column in CELL_COLUMNS) == 0:
                raise DataValidationError(
                    f"line {line_number}: "
                    "cell-count total must be greater than zero"
                )

            sample = row["sample"]

            if sample in seen_samples:
                raise DataValidationError(
                    f"line {line_number}: "
                    f"duplicate sample id {sample}"
                )

            seen_samples.add(sample)

            subject = row["subject"]
            metadata = tuple(row[field] for field in SUBJECT_FIELDS)

            if subject in subject_metadata:
                if subject_metadata[subject] != metadata:
                    raise DataValidationError(
                        f"line {line_number}: "
                        f"inconsistent metadata for subject {subject}"
                    )
            else:
                subject_metadata[subject] = metadata

            rows.append(row)

    if not rows:
        raise DataValidationError("input contains no data rows")

    return rows
